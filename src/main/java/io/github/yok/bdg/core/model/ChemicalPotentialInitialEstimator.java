package io.github.yok.bdg.core.model;

import static com.google.common.base.Preconditions.*;
import io.github.yok.bdg.core.lattice.SquareLattice2D;

/**
 * 固定粒子数（総粒子数 N）に整合する化学ポテンシャル {@code mu-initial}（平均 μ）の初期値を、 BdG 対角化なしで見積もるクラスです。
 *
 * <p>
 * 外部ポテンシャル {@code V(i)} は {@link Potential#valueAt(int)} から取得し、 局所密度近似（LDA）により局所化学ポテンシャル
 * {@code muLoc(i)=mu - V(i)} を作ります。 2D 正方格子分散の k 和により一様系密度を評価し、サイト和で {@code N(mu)} を計算して、 二分法で
 * {@code N(mu)=Ntarget} を満たす {@code mu} を求めます。
 * </p>
 *
 * <h2>近似</h2>
 * <ul>
 * <li>既定（非相互作用 LDA）： {@code n_hom(muLoc,T)=2*mean_k f(eps_k - muLoc)}（係数 2 はスピン和）。</li>
 * <li>オプション（Δなし Hartree-LDA）： {@code muEff = muLoc - U*(n/2)} を用い、サイトごとに
 * {@code n = 2*mean_k f(eps_k - muEff)} を固定点反復で解きます（対角化は不要）。 {@code hartreeEnabled=true} の場合のみ有効です。
 * </li>
 * </ul>
 *
 * <h2>用途</h2>
 * <ul>
 * <li>SCF 前の seed（精密解ではない）</li>
 * <li>外側 μ 探索が「挟み込み」を作りやすい初期 μ を安価に与える</li>
 * </ul>
 *
 * <p>
 * 注意：本クラスの k 和は周期境界の波数 {@code kx=2π nx/Lx, ky=2π ny/Ly} を用いる想定です。
 * </p>
 */
public final class ChemicalPotentialInitialEstimator {

    /**
     * 2D 正方格子です。
     *
     * <p>
     * k 和（分散 {@code epsK}）の構築に格子サイズ {@code Lx, Ly} を使用し、 サイト和（{@code N(mu)=sum_i ...}）の計算でサイト総数
     * {@code Nsite} を使用します。
     * </p>
     */
    private final SquareLattice2D lattice;

    /**
     * 外部ポテンシャルです。
     *
     * <p>
     * サイト {@code i} に対して {@link Potential#valueAt(int)} が {@code V(i)} を返す想定です。
     * 推定ロジックはポテンシャルのパラメータ（例: {@code omegaHo}）を直接参照せず、 {@code valueAt} の結果のみを用います。
     * </p>
     */
    private final Potential potential;

    /**
     * ホッピング {@code t} です。
     *
     * <p>
     * 2D 正方格子（最近接）分散 {@code eps_k = -2t(cos kx + cos ky)} の係数です。 {@code t>0} を前提とします。
     * </p>
     */
    private final double t;

    /**
     * 相互作用 {@code U} です。
     *
     * <p>
     * {@code hartreeEnabled=true} の場合のみ参照します。
     * </p>
     */
    private final double u;

    /**
     * 温度 {@code T} です。
     *
     * <p>
     * フェルミ分布 {@code f(x)=1/(exp(x/T)+1)} の幅を与えます。 数値安定のため {@code x/T} が十分大きい/小さい領域では 0/1 に飽和させます。
     * </p>
     */
    private final double temperature;

    /**
     * 二分法の最大反復回数です。
     *
     * <p>
     * {@code N(mu)} は {@code mu} に対して単調増加であることを前提に二分法で解きます。 初期値推定用途のため、厳密解ではなく上限回数で打ち切る設計です。
     * </p>
     */
    private final int maxIterations;

    /**
     * 早期終了の許容誤差です。
     *
     * <p>
     * {@code |N(mu)-Ntarget| <= tolerance} を満たした時点で探索を終了します。
     * </p>
     */
    private final double tolerance;

    /**
     * 挟み込み区間作りの余裕（マージン）です。
     *
     * <p>
     * 初期の探索区間を {@code [ -4t - Vmax - bracketMargin, +4t + Vmax + bracketMargin ]}
     * のように広げ、目標粒子数を挟み込める確率を上げます。
     * </p>
     */
    private final double bracketMargin;

    /**
     * Δなし Hartree-LDA を有効にするかどうかです。
     */
    private final boolean hartreeEnabled;

    /**
     * Hartree 固定点反復の最大反復回数です。
     *
     * <p>
     * {@code hartreeEnabled=false} の場合は参照しません。
     * </p>
     */
    private final int hartreeMaxIterations;

    /**
     * Hartree 固定点反復の収束判定閾値です。
     *
     * <p>
     * {@code hartreeEnabled=false} の場合は参照しません。
     * </p>
     */
    private final double hartreeTolerance;

    /**
     * Hartree 固定点反復の緩和係数です（0より大きく1以下）。
     *
     * <p>
     * {@code hartreeEnabled=false} の場合は参照しません。
     * </p>
     */
    private final double hartreeMixing;

    /**
     * k 点分散 {@code eps_k} の前計算配列です。
     *
     * <p>
     * 2D 正方格子（周期境界）に対して {@code eps_k = -2t(cos kx + cos ky)} を {@code kx=2π nx/Lx},
     * {@code ky=2π ny/Ly} で評価し、配列（サイズ {@code Lx*Ly}）として保持します。
     * </p>
     */
    private final double[] epsK;

    /**
     * サイトポテンシャル {@code V(i)} の前計算配列です。
     *
     * <p>
     * 各サイト {@code i} に対して {@link Potential#valueAt(int)} を 1 回だけ呼び出し、 {@code V(i)}（サイズ
     * {@code Nsite}）として保持します。
     * </p>
     */
    private final double[] vSite;

    /**
     * 推定器を生成します（前計算を含む）。
     *
     * @param lattice 2D 正方格子（null 不可）
     * @param potential 外部ポテンシャル（null 不可）
     * @param t ホッピング（正）
     * @param u 相互作用 U（0 可）
     * @param temperature 温度（正）
     * @param maxIterations 二分法の最大反復回数（正）
     * @param tolerance 早期終了許容誤差（正）
     * @param bracketMargin 挟み込み余裕（0以上）
     * @param hartreeEnabled Δなし Hartree-LDA を有効にするかどうか
     * @param hartreeMaxIterations Hartree 固定点反復の最大反復回数（正。Hartree 無効時は未参照）
     * @param hartreeTolerance Hartree 固定点反復の収束判定閾値（正。Hartree 無効時は未参照）
     * @param hartreeMixing Hartree 固定点反復の緩和係数（0より大きく1以下。Hartree 無効時は未参照）
     * @throws IllegalArgumentException 引数が不正な場合
     */
    public ChemicalPotentialInitialEstimator(SquareLattice2D lattice, Potential potential, double t,
            double u, double temperature, int maxIterations, double tolerance, double bracketMargin,
            boolean hartreeEnabled, int hartreeMaxIterations, double hartreeTolerance,
            double hartreeMixing) {
        this.lattice = checkNotNull(lattice, "lattice は null 不可です");
        this.potential = checkNotNull(potential, "potential は null 不可です");

        checkArgument(t > 0.0, "t は 0 より大きい必要があります: %s", t);
        checkArgument(temperature > 0.0, "temperature は 0 より大きい必要があります: %s", temperature);
        checkArgument(maxIterations > 0, "maxIterations は 1 以上である必要があります: %s", maxIterations);
        checkArgument(tolerance > 0.0, "tolerance は 0 より大きい必要があります: %s", tolerance);
        checkArgument(bracketMargin >= 0.0, "bracketMargin は 0 以上である必要があります: %s", bracketMargin);

        this.t = t;
        this.u = u;
        this.temperature = temperature;
        this.maxIterations = maxIterations;
        this.tolerance = tolerance;
        this.bracketMargin = bracketMargin;

        // Hartree は「有効フラグ && U!=0」のときだけ動かす
        boolean enabled = hartreeEnabled && (u != 0.0);
        this.hartreeEnabled = enabled;

        if (enabled) {
            checkArgument(hartreeMaxIterations > 0, "hartreeMaxIterations は 1 以上である必要があります: %s",
                    hartreeMaxIterations);
            checkArgument(hartreeTolerance > 0.0, "hartreeTolerance は 0 より大きい必要があります: %s",
                    hartreeTolerance);
            checkArgument(hartreeMixing > 0.0 && hartreeMixing <= 1.0,
                    "hartreeMixing は (0,1] の範囲が必要です: %s", hartreeMixing);

            this.hartreeMaxIterations = hartreeMaxIterations;
            this.hartreeTolerance = hartreeTolerance;
            this.hartreeMixing = hartreeMixing;
        } else {
            // 無効時は参照しない値（呼び出し側にダミーを要求しないため、クラス内で固定）
            this.hartreeMaxIterations = 0;
            this.hartreeTolerance = 0.0;
            this.hartreeMixing = 0.0;
        }

        // 前計算
        this.epsK = precomputeEpsK(this.lattice, this.t);
        this.vSite = precomputeVSite(this.lattice, this.potential);
    }

    /**
     * 目標総粒子数（スピン和）{@code targetTotalN} を満たす {@code mu-initial} を推定します。
     *
     * @param targetTotalN 目標総粒子数（スピン和）
     * @return 推定された μ（平均 μ の seed）
     * @throws IllegalArgumentException {@code targetTotalN} が範囲外の場合
     * @throws IllegalStateException 挟み込み区間の生成に失敗した場合
     */
    public double estimateMuInitial(double targetTotalN) {
        checkArgument(targetTotalN > 0.0, "targetTotalN は 0 より大きい必要があります: %s", targetTotalN);

        int nSite = lattice.siteCount();
        double nUpperBound = 2.0 * nSite; // 各サイト 0..2（スピン和）
        checkArgument(targetTotalN <= nUpperBound, "targetTotalN が上限（2*サイト数）を超えています: %s > %s",
                targetTotalN, nUpperBound);

        double vMax = max(vSite);

        double muMin = -4.0 * t - vMax - bracketMargin;
        double muMax = +4.0 * t + vMax + bracketMargin;

        double nAtMin = totalN(muMin);
        double nAtMax = totalN(muMax);

        // 念のためレンジ拡張（通常は不要想定）
        int expand = 0;
        while ((targetTotalN < nAtMin || targetTotalN > nAtMax) && expand < 20) {
            double span = (muMax - muMin);
            muMin -= span;
            muMax += span;
            nAtMin = totalN(muMin);
            nAtMax = totalN(muMax);
            expand++;
        }

        if (targetTotalN < nAtMin || targetTotalN > nAtMax) {
            throw new IllegalStateException("目標粒子数を挟み込む探索区間を作れませんでした。" + " targetTotalN="
                    + targetTotalN + ", N(muMin)=" + nAtMin + ", N(muMax)=" + nAtMax + ", muMin="
                    + muMin + ", muMax=" + muMax);
        }

        double lo = muMin;
        double hi = muMax;

        for (int it = 0; it < maxIterations; it++) {
            double mid = 0.5 * (lo + hi);
            double nMid = totalN(mid);
            double err = nMid - targetTotalN;

            if (Math.abs(err) <= tolerance) {
                return mid;
            }
            if (err < 0.0) {
                lo = mid;
            } else {
                hi = mid;
            }
        }

        return 0.5 * (lo + hi);
    }

    /**
     * 指定した μ に対する総粒子数（スピン和）{@code N(mu)} を計算します。
     *
     * @param mu 平均化学ポテンシャル
     * @return 総粒子数（スピン和）
     */
    private double totalN(double mu) {
        double sum = 0.0;
        for (int i = 0; i < vSite.length; i++) {
            double muMinusV = mu - vSite[i]; // μ − V(i)
            sum += siteDensity(muMinusV);
        }
        return sum;
    }

    /**
     * サイト i の密度（スピン和、0〜2）を計算します。
     *
     * <p>
     * Hartree が無効な場合は非相互作用 LDA（k 和）です。 Hartree が有効な場合は Δなし Hartree-LDA を固定点反復で解きます。
     * </p>
     *
     * @param muMinusV {@code μ − V(i)} です
     * @return サイト密度（スピン和、0〜2）です
     */
    private double siteDensity(double muMinusV) {
        // 非相互作用密度を初期値にする
        double n = nHom(muMinusV);

        if (!hartreeEnabled) {
            return clamp02(n);
        }

        // muEff = (μ − V) − U*(n/2)
        for (int it = 0; it < hartreeMaxIterations; it++) {
            double muEff = muMinusV - u * (0.5 * n);
            double nNew = nHom(muEff);

            double nMix = (1.0 - hartreeMixing) * n + hartreeMixing * nNew;
            if (Math.abs(nMix - n) <= hartreeTolerance) {
                return clamp02(nMix);
            }
            n = nMix;
        }

        return clamp02(n);
    }

    /**
     * 一様 2D 正方格子のスピン和密度 {@code n_hom(mu,T)} を計算します（非相互作用）。
     *
     * @param muLoc 化学ポテンシャル（例：{@code μ − V} や {@code μ_eff}）
     * @return スピン和密度（0〜2）
     */
    private double nHom(double muLoc) {
        double acc = 0.0;
        for (double e : epsK) {
            acc += fermi(e - muLoc, temperature);
        }
        return 2.0 * (acc / epsK.length);
    }

    /**
     * フェルミ分布 {@code f(x)=1/(exp(x/T)+1)} を数値安定に評価します。
     *
     * @param x エネルギー差（例: {@code eps_k - muLoc}）
     * @param T 温度（正）
     * @return フェルミ分布値（0〜1）
     */
    private static double fermi(double x, double T) {
        double y = x / T;
        if (y > 50.0) {
            return 0.0;
        }
        if (y < -50.0) {
            return 1.0;
        }
        return 1.0 / (Math.exp(y) + 1.0);
    }

    /**
     * 値を 0.0〜2.0 の範囲にクランプします。
     *
     * @param v 値です
     * @return クランプ後の値です
     */
    private static double clamp02(double v) {
        if (v < 0.0) {
            return 0.0;
        }
        if (v > 2.0) {
            return 2.0;
        }
        return v;
    }

    /**
     * 2D 正方格子（周期境界）の分散 {@code eps_k} を前計算します。
     *
     * @param lattice 2D 正方格子
     * @param t ホッピング
     * @return 分散配列（サイズ: {@code Lx*Ly}）
     */
    private static double[] precomputeEpsK(SquareLattice2D lattice, double t) {
        int lx = lattice.lx();
        int ly = lattice.ly();

        double[] eps = new double[lx * ly];
        int idx = 0;

        for (int ny = 0; ny < ly; ny++) {
            double ky = 2.0 * Math.PI * ny / ly;
            double cy = Math.cos(ky);

            for (int nx = 0; nx < lx; nx++) {
                double kx = 2.0 * Math.PI * nx / lx;
                double cx = Math.cos(kx);

                eps[idx++] = -2.0 * t * (cx + cy);
            }
        }

        return eps;
    }

    /**
     * サイトごとの外部ポテンシャル {@code V(i)} を前計算します。
     *
     * @param lattice 2D 正方格子
     * @param potential 外部ポテンシャル
     * @return サイトポテンシャル配列（サイズ: {@code Nsite}）
     */
    private static double[] precomputeVSite(SquareLattice2D lattice, Potential potential) {
        int n = lattice.siteCount();
        double[] v = new double[n];
        for (int i = 0; i < n; i++) {
            v[i] = potential.valueAt(i);
        }
        return v;
    }

    /**
     * 配列の最大値を返します。
     *
     * @param a 対象配列（長さ 1 以上を想定）
     * @return 最大値
     */
    private static double max(double[] a) {
        double m = a[0];
        for (int i = 1; i < a.length; i++) {
            m = Math.max(m, a[i]);
        }
        return m;
    }
}
