package io.github.yok.bdg.core.solver;

import com.google.common.base.Preconditions;
import io.github.yok.bdg.core.state.MeanFieldState;
import java.util.Locale;
import lombok.RequiredArgsConstructor;
import lombok.extern.slf4j.Slf4j;

/**
 * 固定した Zeeman 場 h のもとで、総粒子数 N が目標値に一致するように 平均化学ポテンシャル μ を二分法で探索します。
 *
 * <p>
 * まず初期μから「必要な方向（片側）」にのみ探索して、目標を挟む区間（挟み込み区間）を作ります。 その後、その区間内で二分法により μ を詰めます。
 * </p>
 */
@Slf4j
@RequiredArgsConstructor
public final class ChemicalPotentialBisectionOptimizer implements ChemicalPotentialOptimizer {

    /**
     * 内側ループ（SCF）です。
     */
    private final SelfConsistencySolver scfSolver;

    /**
     * μ探索の最大反復回数です。
     */
    private final int maxMuIterations;

    /**
     * 粒子数の許容誤差です。
     */
    private final double particleTolerance;

    /**
     * μ探索の探索下限です（安全側の丸めに使います）。
     */
    private final double muMin = -100.0;

    /**
     * μ探索の探索上限です（安全側の丸めに使います）。
     */
    private final double muMax = 100.0;

    /**
     * 片側探索での初期移動幅です。
     */
    private final double initialSearchWidth = 0.1;

    /**
     * 片側探索で幅を広げる回数上限です。
     */
    private final int maxWidthExpansions = 20;

    /**
     * SCFが未収束になった点を避けるため、収束側へ寄せて「収束する点」を探す最大試行回数です。
     */
    private final int maxBackoffStepsWhenNotConverged = 30;

    /**
     * 固定した h のもとで、総粒子数 N が目標値に一致するように μ を探索します。
     *
     * @param initialState 初期状態（SCFの初期値）です
     * @param totalParticleNumber 目標の総粒子数 N (=N1+N2) です
     * @param zeemanField Zeeman 場 h（有効磁場、h=(mu1-mu2)/2）です
     * @return 最適化結果です
     * @throws NullPointerException 初期状態が null の場合
     * @throws IllegalArgumentException 総粒子数が正でない、または許容誤差が正でない場合
     * @throws IllegalStateException 初期点でSCFが収束しない、または挟み込み区間が作れない場合
     */
    @Override
    public OptimizationResult optimize(MeanFieldState initialState, int totalParticleNumber,
            double zeemanField) {
        Preconditions.checkNotNull(initialState, "初期状態が null です。");
        Preconditions.checkArgument(totalParticleNumber > 0, "総粒子数Nは正の値である必要があります。N=%s",
                totalParticleNumber);
        Preconditions.checkArgument(particleTolerance > 0.0, "粒子数の許容誤差は正の値である必要があります。tol=%s",
                particleTolerance);

        final double targetN = totalParticleNumber;
        final double h = zeemanField;

        // 初期中心 μ0（warm start：initialState の μ平均）
        final double mu0 =
                clampMuToSearchRange(0.5 * (initialState.getMu1() + initialState.getMu2()));

        log.info("μ探索を開始します。目標N={}、h={}、初期μ={}、探索範囲=[{}, {}]、許容誤差={}、最大反復回数={}", fmt5(targetN),
                fmt5(h), fmt5(mu0), fmt5(muMin), fmt5(muMax), fmt5(particleTolerance),
                maxMuIterations);

        // 初期点を評価（ここが収束しないと warm start の起点が作れない）
        Evaluation start = evaluateAtMuWithFixedZeemanField(initialState, mu0, h);
        if (!start.converged) {
            throw new IllegalStateException("初期μでSCFが収束しませんでした。mu0=" + mu0 + ", h=" + h);
        }

        // 初期点が既に許容誤差内なら終了
        if (start.absErr(targetN) < particleTolerance) {
            log.info("初期μで目標に一致しました。h={}、μ={}、N={}、誤差={}", fmt5(h), fmt5(start.muAvg),
                    fmt5(start.totalN), fmt5(start.absErr(targetN)));
            return new OptimizationResult(start.mu1, start.mu2, start.state, 0);
        }

        // 片側探索で挟み込み区間（収束点2つ、符号反転）を作る
        EnclosingInterval interval = findEnclosingIntervalByOneSidedSearch(start, targetN, h);

        Evaluation lower = interval.lower;
        Evaluation upper = interval.upper;

        Evaluation best = chooseBetterByParticleError(lower, upper, targetN);

        // 粒子数 N が変わらない plateau を検出して無駄計算を止める
        double prevN = Double.NaN;
        final double plateauEps = 0.01;

        // 挟み込み区間内で二分法
        for (int iter = 1; iter <= maxMuIterations; iter++) {
            final double muMid = 0.5 * (lower.muAvg + upper.muAvg);

            // μmid に近い側の収束済み状態を seed にする
            final Evaluation near =
                    (Math.abs(muMid - lower.muAvg) <= Math.abs(upper.muAvg - muMid)) ? lower
                            : upper;
            final Evaluation far = (near == lower) ? upper : lower;

            Evaluation mid = evaluateAtMuWithFixedZeemanField(near.state, muMid, h);

            logMuIteration(iter, targetN, h, lower, upper, mid, near);

            // lower/upper は収束点として維持する（未収束点を採用しない）
            if (!mid.converged) {
                EnclosingInterval repaired = shrinkIntervalKeepingConvergedEndpoints(lower, upper,
                        near, far, muMid, targetN, h);
                lower = repaired.lower;
                upper = repaired.upper;
                best = updateBestIfBetter(best, lower, targetN);
                best = updateBestIfBetter(best, upper, targetN);
                continue;
            }

            best = updateBestIfBetter(best, mid, targetN);

            // μを動かしてもNが変わらない（plateau）場合の打ち切り
            if (Double.isFinite(prevN) && Math.abs(mid.totalN - prevN) < plateauEps) {
                log.warn("μを更新しても N が変化しないため探索を打ち切ります。反復回数={}、h={}、μ={}、N={}、誤差={}、最良誤差={}", iter,
                        fmt5(h), fmt5(mid.muAvg), fmt5(mid.totalN), fmt5(mid.absErr(targetN)),
                        fmt5(best.absErr(targetN)));
                return new OptimizationResult(best.mu1, best.mu2, best.state, iter);
            }
            prevN = mid.totalN;

            if (mid.absErr(targetN) <= particleTolerance) {
                log.info("μ探索を終了します（許容誤差内）。反復回数={}、h={}、μ={}、N={}、誤差={}", iter, fmt5(h),
                        fmt5(mid.muAvg), fmt5(mid.totalN), fmt5(mid.absErr(targetN)));
                return new OptimizationResult(mid.mu1, mid.mu2, mid.state, iter);
            }

            // 二分更新：符号で区間を更新（lower/upperは収束点）
            final double fLower = lower.f(targetN);
            final double fMid = mid.f(targetN);

            if (fLower * fMid <= 0.0) {
                upper = mid;
            } else {
                lower = mid;
            }

            best = updateBestIfBetter(best, lower, targetN);
            best = updateBestIfBetter(best, upper, targetN);
        }

        log.info("μ探索を終了します（上限到達）。h={}、最良μ={}、最良N={}、最良誤差={}", fmt5(h), fmt5(best.muAvg),
                fmt5(best.totalN), fmt5(best.absErr(targetN)));

        return new OptimizationResult(best.mu1, best.mu2, best.state, maxMuIterations);
    }

    /**
     * 初期点から「必要な方向（片側）」へ探索し、f(μ)=N(μ;h)-N_target が符号反転する 収束点2つ（挟み込み区間）を作ります。
     *
     * @param start 初期点（収束済み）です
     * @param targetN 目標の総粒子数です
     * @param h 固定のZeeman場です
     * @return 挟み込み区間（lower/upperはいずれも収束点）です
     * @throws IllegalStateException 区間を作れない場合
     */
    private EnclosingInterval findEnclosingIntervalByOneSidedSearch(Evaluation start,
            double targetN, double h) {
        final double f0 = start.f(targetN);

        // f0<0 なら Nが不足 → μを増やす方向
        // f0>0 なら Nが過剰 → μを減らす方向
        final int direction = (f0 < 0.0) ? +1 : -1;
        final String dirLabel = (direction > 0) ? "μを増やす方向（Nを増やす側）" : "μを減らす方向（Nを減らす側）";

        log.info("挟み込み区間の探索を開始します。初期N={}、目標N={}、差分f={}、探索方向={}", fmt5(start.totalN), fmt5(targetN),
                fmt5(f0), dirLabel);

        // 直前の収束点（secant用）
        Evaluation prev = null;
        Evaluation current = start;
        double width = initialSearchWidth;

        for (int i = 0; i < maxWidthExpansions; i++) {
            double muTryRaw =
                    proposeNextMuForOneSidedSearch(prev, current, targetN, direction, width);
            double muTry = clampMuToSearchRange(muTryRaw);

            // 丸めで進まなくなったら、これ以上外へ行けない
            if (muTry == current.muAvg) {
                throw new IllegalStateException(
                        "探索範囲の端に到達しましたが、挟み込み区間を作れませんでした。" + " direction=" + dirLabel + ", mu="
                                + current.muAvg + ", range=[" + muMin + "," + muMax + "]");
            }

            Evaluation tried = evaluateAtMuWithFixedZeemanField(current.state, muTry, h);

            logOneSidedTry(i + 1, targetN, h, current, muTry, tried, dirLabel);

            if (!tried.converged) {
                Evaluation recovered = findConvergedPointBetween(current, muTry, targetN, h);
                if (recovered == null) {
                    log.warn("片側探索で未収束が続き、収束点を前進できませんでした。h={}、現在μ={}、試行μ={}", fmt5(h),
                            fmt5(current.muAvg), fmt5(muTry));
                    width *= 2.0;
                    continue;
                }
                // 収束点を少しでも前進できたら、それを基準に続行（prev を更新）
                prev = current;
                current = recovered;
                width *= 2.0;
                continue;
            }

            // 収束したので、符号反転するなら挟み込み確定
            double fCur = current.f(targetN);
            double fTry = tried.f(targetN);

            if (fCur * fTry <= 0.0) {
                EnclosingInterval interval = orderAsLowerUpper(current, tried);
                log.info("挟み込み区間を確定しました。h={}、下側μ={} (N={})、上側μ={} (N={})", fmt5(h),
                        fmt5(interval.lower.muAvg), fmt5(interval.lower.totalN),
                        fmt5(interval.upper.muAvg), fmt5(interval.upper.totalN));
                return interval;
            }

            // 符号が変わらないなら、前進して続行（prev を更新）
            prev = current;
            current = tried;
            width *= 2.0;
        }

        throw new IllegalStateException("片側探索で挟み込み区間を作れませんでした。h=" + h + ", targetN=" + targetN
                + ", startMu=" + start.muAvg + ", range=[" + muMin + "," + muMax + "]");
    }

    /**
     * current（収束点）と、未収束または外側の点 muTarget の間で、 収束する点を見つけて「少しでも外側へ前進」させます。
     *
     * <p>
     * 見つかった場合は muTarget に近い収束点を返します。 見つからない場合は null を返します。
     * </p>
     *
     * @param current 現在の収束点です
     * @param muTarget 外側の目標μです
     * @param targetN 目標の総粒子数です
     * @param h 固定のZeeman場です
     * @return 前進できた収束点、見つからない場合は null です
     */
    private Evaluation findConvergedPointBetween(Evaluation current, double muTarget,
            double targetN, double h) {
        double a = current.muAvg;
        double b = muTarget;

        // 方向が逆なら何もしない
        if (a == b) {
            return null;
        }

        Evaluation best = null;

        for (int step = 0; step < maxBackoffStepsWhenNotConverged; step++) {
            double muTry = 0.5 * (a + b);

            Evaluation ev = evaluateAtMuWithFixedZeemanField(current.state, muTry, h);

            log.debug("未収束回避の再探索：h={}、μ={}、SCF={}、N={}、誤差={}", fmt5(h), fmt5(muTry),
                    ev.converged ? "収束" : "未収束", fmt5(ev.totalN), fmt5(ev.absErr(targetN)));

            if (ev.converged) {
                best = ev;
                a = muTry; // さらに外側へ寄せて、もっと進めるか試す
            } else {
                b = muTry; // 収束側へ戻す
            }
        }

        return best;
    }

    /**
     * mid が未収束だった場合に、収束側へ寄せて区間を縮めつつ、 lower/upper が「収束点のまま」になるように挟み込み区間を修復します。
     *
     * @param lower 現在の下側（収束点）です
     * @param upper 現在の上側（収束点）です
     * @param near seedに使った側（収束側）です
     * @param far 反対側（収束点）です
     * @param muMid 未収束になったμです
     * @param targetN 目標の総粒子数です
     * @param h 固定のZeeman場です
     * @return 修復後の挟み込み区間です（lower/upperはいずれも収束点）
     */
    private EnclosingInterval shrinkIntervalKeepingConvergedEndpoints(Evaluation lower,
            Evaluation upper, Evaluation near, Evaluation far, double muMid, double targetN,
            double h) {

        double a = near.muAvg;
        double b = muMid;

        Evaluation recovered = null;

        for (int step = 0; step < maxBackoffStepsWhenNotConverged; step++) {
            double muTry = 0.5 * (a + b);

            Evaluation ev = evaluateAtMuWithFixedZeemanField(near.state, muTry, h);

            log.debug("未収束回避の再試行：h={}、μ={}、SCF={}、N={}、誤差={}", fmt5(h), fmt5(muTry),
                    ev.converged ? "収束" : "未収束", fmt5(ev.totalN), fmt5(ev.absErr(targetN)));

            if (ev.converged) {
                recovered = ev;
                break;
            }
            b = muTry;
        }

        if (recovered == null) {
            log.warn("未収束回避に失敗しました。h={}、収束側μ={}、未収束μ={}。区間は維持して続行します。", fmt5(h), fmt5(near.muAvg),
                    fmt5(muMid));
            return new EnclosingInterval(lower, upper);
        }

        double fNear = near.f(targetN);
        double fFar = far.f(targetN);
        double fRec = recovered.f(targetN);

        boolean nearAndRec = (fNear * fRec <= 0.0);
        boolean recAndFar = (fRec * fFar <= 0.0);

        if (nearAndRec) {
            return orderAsLowerUpper(near, recovered);
        }
        if (recAndFar) {
            return orderAsLowerUpper(recovered, far);
        }

        log.warn("再取得した収束点で挟み込み区間を維持できませんでした。h={}。区間は維持して続行します。", fmt5(h));
        return new EnclosingInterval(lower, upper);
    }

    /**
     * 固定した h のもとで、指定した μ に対して SCF を実行し、 粒子数などの評価結果を返します。
     *
     * @param seed SCFの初期値（warm start）です
     * @param muAvg 平均化学ポテンシャル μ です
     * @param h 固定のZeeman場です
     * @return 評価結果です（収束可否、粒子数、状態など）
     */
    private Evaluation evaluateAtMuWithFixedZeemanField(MeanFieldState seed, double muAvg,
            double h) {
        Preconditions.checkNotNull(seed, "SCFのseedが null です。");

        double mu1 = muAvg + h;
        double mu2 = muAvg - h;

        MeanFieldState st = copyState(seed);
        st.setMu1(mu1);
        st.setMu2(mu2);

        SelfConsistencySolver.SolveResult sr = scfSolver.solve(st);

        double n1 = sumArray(st.getN1());
        double n2 = sumArray(st.getN2());
        double total = n1 + n2;

        return new Evaluation(muAvg, mu1, mu2, total, sr.isConverged(), st);
    }

    /**
     * μ探索（二分法）の各反復の状況を、日本語で分かるログとして出力します。
     *
     * @param iter 反復回数です
     * @param targetN 目標の総粒子数です
     * @param h 固定のZeeman場です
     * @param lower 下側の収束点です
     * @param upper 上側の収束点です
     * @param mid 中点の評価結果です
     * @param seedFrom seedに使った側（lowerまたはupper）です
     */
    private void logMuIteration(int iter, double targetN, double h, Evaluation lower,
            Evaluation upper, Evaluation mid, Evaluation seedFrom) {

        String seedSide = (seedFrom == lower) ? "下側（低μ）" : "上側（高μ）";

        if (mid.converged) {
            log.info("反復{}：h={}、区間=[{}, {}]、μ={}（seed={}）→ SCF=収束、N={}、誤差={}、μ1={}、μ2={}", iter,
                    fmt5(h), fmt5(lower.muAvg), fmt5(upper.muAvg), fmt5(mid.muAvg), seedSide,
                    fmt5(mid.totalN), fmt5(mid.absErr(targetN)), fmt5(mid.mu1), fmt5(mid.mu2));
        } else {
            log.info("反復{}：h={}、区間=[{}, {}]、μ={}（seed={}）→ SCF=未収束", iter, fmt5(h),
                    fmt5(lower.muAvg), fmt5(upper.muAvg), fmt5(mid.muAvg), seedSide);
        }
    }

    /**
     * 片側探索の試行ログを、日本語で分かる形に出力します。
     *
     * @param attempt 試行番号です
     * @param targetN 目標の総粒子数です
     * @param h 固定のZeeman場です
     * @param current 現在の収束点です
     * @param muTry 試行したμです
     * @param tried 試行結果です
     * @param dirLabel 探索方向の説明文です
     */
    private void logOneSidedTry(int attempt, double targetN, double h, Evaluation current,
            double muTry, Evaluation tried, String dirLabel) {

        if (tried.converged) {
            log.info("区間探索{}：方向={}、h={}、基準μ={} → 試行μ={}：SCF=収束、N={}、誤差={}", attempt, dirLabel,
                    fmt5(h), fmt5(current.muAvg), fmt5(muTry), fmt5(tried.totalN),
                    fmt5(tried.absErr(targetN)));
        } else {
            log.info("区間探索{}：方向={}、h={}、基準μ={} → 試行μ={}：SCF=未収束", attempt, dirLabel, fmt5(h),
                    fmt5(current.muAvg), fmt5(muTry));
        }
    }

    /**
     * 2点をμの大小で並べ替えて、下側/上側として返します。
     *
     * @param a 1点目です
     * @param b 2点目です
     * @return 下側/上側に整列した挟み込み区間です
     */
    private EnclosingInterval orderAsLowerUpper(Evaluation a, Evaluation b) {
        return (a.muAvg <= b.muAvg) ? new EnclosingInterval(a, b) : new EnclosingInterval(b, a);
    }

    /**
     * 現在の最良値より誤差が小さい場合に、最良値を更新します。
     *
     * @param best 現在の最良値です
     * @param candidate 比較対象です
     * @param targetN 目標の総粒子数です
     * @return 更新後の最良値です
     */
    private Evaluation updateBestIfBetter(Evaluation best, Evaluation candidate, double targetN) {
        if (candidate == null || !candidate.converged) {
            return best;
        }
        if (best == null) {
            return candidate;
        }
        return (candidate.absErr(targetN) < best.absErr(targetN)) ? candidate : best;
    }

    /**
     * 2点のうち、粒子数誤差が小さい方を返します（いずれも収束点であることを想定）。
     *
     * @param a 比較対象1です
     * @param b 比較対象2です
     * @param targetN 目標の総粒子数です
     * @return 誤差が小さい方です
     */
    private Evaluation chooseBetterByParticleError(Evaluation a, Evaluation b, double targetN) {
        if (a == null) {
            return b;
        }
        if (b == null) {
            return a;
        }
        return (a.absErr(targetN) <= b.absErr(targetN)) ? a : b;
    }

    /**
     * μを探索範囲に丸め込みます。
     *
     * @param mu μです
     * @return 範囲内に丸めたμです
     */
    private double clampMuToSearchRange(double mu) {
        if (mu < muMin) {
            return muMin;
        }
        if (mu > muMax) {
            return muMax;
        }
        return mu;
    }

    /**
     * 配列の総和を計算します。
     *
     * @param a 配列です
     * @return 総和です
     */
    private static double sumArray(double[] a) {
        double s = 0.0;
        for (double v : a) {
            s += v;
        }
        return s;
    }

    /**
     * 平均場状態をコピーします（Δ, n1, n2, μ1, μ2 を複製します）。
     *
     * @param src 元の状態です
     * @return コピーした状態です
     */
    private static MeanFieldState copyState(MeanFieldState src) {
        int n = src.getSiteCount();
        MeanFieldState dst = new MeanFieldState(n);
        System.arraycopy(src.getDelta(), 0, dst.getDelta(), 0, n);
        System.arraycopy(src.getN1(), 0, dst.getN1(), 0, n);
        System.arraycopy(src.getN2(), 0, dst.getN2(), 0, n);
        dst.setMu1(src.getMu1());
        dst.setMu2(src.getMu2());
        return dst;
    }

    /**
     * 数値を小数点以下5桁までの文字列に整形します。
     *
     * @param v 数値です
     * @return 整形文字列です
     */
    private static String fmt5(double v) {
        return String.format(Locale.ROOT, "%.5f", v);
    }

    /**
     * 収束点2つで構成される「挟み込み区間（下側/上側）」です。
     *
     * <p>
     * ここでの挟み込み区間とは、f(μ)=N(μ;h)-N_target の符号が 下側と上側で反転している区間を指します。
     * </p>
     */
    private static final class EnclosingInterval {

        /**
         * 下側（μが小さい側）の収束点です。
         */
        final Evaluation lower;

        /**
         * 上側（μが大きい側）の収束点です。
         */
        final Evaluation upper;

        /**
         * 挟み込み区間を作成します。
         *
         * @param lower 下側の収束点です
         * @param upper 上側の収束点です
         */
        EnclosingInterval(Evaluation lower, Evaluation upper) {
            this.lower = lower;
            this.upper = upper;
        }
    }

    /**
     * ある μ に対する評価結果です（SCFの収束可否と粒子数など）。
     */
    private static final class Evaluation {

        /**
         * 平均化学ポテンシャル μ です。
         */
        final double muAvg;

        /**
         * μ1（mu+h）です。
         */
        final double mu1;

        /**
         * μ2（mu-h）です。
         */
        final double mu2;

        /**
         * 総粒子数 N=N1+N2 です。
         */
        final double totalN;

        /**
         * SCFが収束したかどうかです。
         */
        final boolean converged;

        /**
         * 収束時（または最終時点）の平均場状態です。
         */
        final MeanFieldState state;

        /**
         * 評価結果を作成します。
         *
         * @param muAvg 平均化学ポテンシャルです
         * @param mu1 μ1です
         * @param mu2 μ2です
         * @param totalN 総粒子数です
         * @param converged SCF収束フラグです
         * @param state 状態です
         */
        Evaluation(double muAvg, double mu1, double mu2, double totalN, boolean converged,
                MeanFieldState state) {
            this.muAvg = muAvg;
            this.mu1 = mu1;
            this.mu2 = mu2;
            this.totalN = totalN;
            this.converged = converged;
            this.state = state;
        }

        /**
         * f(μ)=N(μ;h)-N_target を返します。
         *
         * @param targetN 目標の総粒子数です
         * @return 差分です
         */
        double f(double targetN) {
            return totalN - targetN;
        }

        /**
         * |N-N_target| を返します。
         *
         * @param targetN 目標の総粒子数です
         * @return 絶対誤差です
         */
        double absErr(double targetN) {
            return Math.abs(totalN - targetN);
        }
    }

    /**
     * 片側探索で次に試す μ を提案します。
     *
     * <p>
     * prev がある場合は secant（割線）で targetN に当てに行きます。 ただし、片側探索の前提（direction
     * 方向に進む）を壊す提案は捨て、通常ステップにフォールバックします。 また、飛びすぎ防止のため current からの移動量を ±width に制限します。
     * </p>
     *
     * @param prev 直前の収束点（初回は null）
     * @param current 現在の収束点
     * @param targetN 目標粒子数
     * @param direction 探索方向（+1: μを増やす, -1: μを減らす）
     * @param width 現在の探索幅（正の値）
     * @return 次に試す μ（丸め前）
     */
    private static double proposeNextMuForOneSidedSearch(Evaluation prev, Evaluation current,
            double targetN, int direction, double width) {

        double muTryRaw;

        if (prev != null) {
            double fPrev = prev.f(targetN);
            double fCur = current.f(targetN);
            double denom = (fCur - fPrev);

            if (Math.abs(denom) > 1e-12) {
                muTryRaw = current.muAvg + (0.0 - fCur) * (current.muAvg - prev.muAvg) / denom;

                // 一方向探索の前提を壊さない（逆方向の提案は捨てる）
                if ((direction > 0 && muTryRaw <= current.muAvg)
                        || (direction < 0 && muTryRaw >= current.muAvg)) {
                    muTryRaw = current.muAvg + direction * width;
                }
            } else {
                muTryRaw = current.muAvg + direction * width;
            }

            // secantが大きく跳ねるのを抑える（±width以内に丸める）
            double step = muTryRaw - current.muAvg;
            if (step > width) {
                step = width;
            }
            if (step < -width) {
                step = -width;
            }
            muTryRaw = current.muAvg + step;

        } else {
            // 初回は従来どおり
            muTryRaw = current.muAvg + direction * width;
        }

        return muTryRaw;
    }
}
