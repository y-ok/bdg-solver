package io.github.yok.bdg.core.solver;

import io.github.yok.bdg.app.BdgProperties;
import io.github.yok.bdg.core.linearalgebra.EigenDecompositionBackend;
import io.github.yok.bdg.core.linearalgebra.EigenDecompositionBackend.EigenDecompositionResult;
import io.github.yok.bdg.core.model.BdgModel;
import io.github.yok.bdg.core.model.MeanFieldUpdater;
import io.github.yok.bdg.core.state.MeanFieldIncrements;
import io.github.yok.bdg.core.state.MeanFieldState;
import java.util.Locale;
import lombok.Getter;
import lombok.Value;
import lombok.extern.slf4j.Slf4j;
import org.ejml.data.DMatrixRMaj;

/**
 * BdG の自己無撞着（SCF）ループを実行するクラスです。
 *
 * <p>
 * 平均場状態 → BdG 行列構築 → 固有分解 → 更新式 → mixing → 収束判定、を反復します。
 * </p>
 */
@Getter
@Slf4j
public final class SelfConsistencySolver {

    /**
     * BdG 行列を構築するモデルです。
     */
    private final BdgModel model;

    /**
     * 固有分解を行うバックエンドです。
     */
    private final EigenDecompositionBackend eigenBackend;

    /**
     * 更新式（n1, n2, Δ）を計算するコンポーネントです。
     */
    private final MeanFieldUpdater updater;

    /**
     * 最大反復回数です。
     */
    private final int maxIterations;

    /**
     * 密度（n1, n2）に適用する mixing 係数です。
     */
    private final double densityAlpha;

    /**
     * 秩序パラメータ（Δ）に適用する mixing 係数です。
     */
    private final double orderParameterAlpha;

    /**
     * 密度（n1, n2）の収束判定に用いる許容誤差（残差の絶対値）です。
     */
    private final double densityToleranceAbs;

    /**
     * 秩序パラメータ（Δ）の収束判定に用いる許容誤差（残差の絶対値）です。
     */
    private final double orderParameterToleranceAbs;

    /**
     * 秩序パラメータ（Δ）の収束判定に用いる許容誤差（残差RMSの絶対値）です。
     */
    private final double orderParameterToleranceRmsAbs;

    /**
     * SCF ソルバを生成します。
     *
     * @param model BdG 行列を構築するモデルです（null 不可）
     * @param eigenBackend 固有分解バックエンドです（null 不可）
     * @param updater 更新式です（null 不可）
     * @param scf SCF 設定です（null 不可）
     * @throws IllegalArgumentException 引数が不正な場合に発生します
     */
    public SelfConsistencySolver(BdgModel model, EigenDecompositionBackend eigenBackend,
            MeanFieldUpdater updater, BdgProperties.Scf scf) {

        if (model == null) {
            throw new IllegalArgumentException("model は null 不可です");
        }
        if (eigenBackend == null) {
            throw new IllegalArgumentException("eigenBackend は null 不可です");
        }
        if (updater == null) {
            throw new IllegalArgumentException("updater は null 不可です");
        }
        if (scf == null) {
            throw new IllegalArgumentException("scf は null 不可です");
        }
        if (scf.getMixing() == null) {
            throw new IllegalArgumentException("scf.mixing は null 不可です");
        }
        if (scf.getTolerance() == null) {
            throw new IllegalArgumentException("scf.tolerance は null 不可です");
        }

        int maxIter = scf.getMaxIter();
        double dAlpha = scf.getMixing().getDensityAlpha();
        double opAlpha = scf.getMixing().getOrderParameterAlpha();
        double dTol = scf.getTolerance().getDensityAbs();
        double opTol = scf.getTolerance().getOrderParameterAbs();
        double opTolRms = scf.getTolerance().getOrderParameterRmsAbs();

        if (maxIter <= 0) {
            throw new IllegalArgumentException("scf.maxIter は 1 以上が必要です: " + maxIter);
        }
        if (!(dAlpha > 0.0 && dAlpha <= 1.0)) {
            throw new IllegalArgumentException("scf.mixing.densityAlpha は (0, 1] が必要です: " + dAlpha);
        }
        if (!(opAlpha > 0.0 && opAlpha <= 1.0)) {
            throw new IllegalArgumentException(
                    "scf.mixing.orderParameterAlpha は (0, 1] が必要です: " + opAlpha);
        }
        if (!(dTol > 0.0)) {
            throw new IllegalArgumentException(
                    "scf.tolerance.densityAbs は 0 より大きい必要があります: " + dTol);
        }
        if (!(opTol > 0.0)) {
            throw new IllegalArgumentException(
                    "scf.tolerance.orderParameterAbs は 0 より大きい必要があります: " + opTol);
        }
        if (!(opTolRms > 0.0)) {
            throw new IllegalArgumentException(
                    "scf.tolerance.orderParameterRmsAbs は 0 より大きい必要があります: " + opTolRms);
        }

        this.model = model;
        this.eigenBackend = eigenBackend;
        this.updater = updater;

        this.maxIterations = maxIter;
        this.densityAlpha = dAlpha;
        this.orderParameterAlpha = opAlpha;
        this.densityToleranceAbs = dTol;
        this.orderParameterToleranceAbs = opTol;
        this.orderParameterToleranceRmsAbs = opTolRms;
    }

    /**
     * SCF ループを実行して、収束した平均場状態を返します。
     *
     * @param state 初期状態です（配列長はモデルのサイト数と一致が必要です）
     * @return 実行結果です
     * @throws IllegalArgumentException state が null、またはサイズ不整合の場合に発生します
     * @throws IllegalStateException モデルが実対称を要求する前提と矛盾する場合に発生します
     */
    public SolveResult solve(MeanFieldState state) {
        if (state == null) {
            throw new IllegalArgumentException("state は null 不可です");
        }

        int siteCount = model.siteCount();
        ensureStateSizeMatchesModel(state, siteCount);

        if (!model.isRealSymmetric()) {
            throw new IllegalStateException("このソルバは実対称 BdG 行列を前提としています");
        }

        long t0 = System.nanoTime();
        double nTotal0 = totalParticleNumber(state);

        log.info(
                "SCFを開始します。最大反復={}、mixing係数（密度）={}、mixing係数（秩序パラメータ）={}、"
                        + "許容誤差（密度）={}、許容誤差（秩序パラメータmax）={}、許容誤差（秩序パラメータRMS）={}、N={}",
                maxIterations, fmt5(densityAlpha), fmt5(orderParameterAlpha),
                fmt5(densityToleranceAbs), fmt5(orderParameterToleranceAbs),
                fmt5(orderParameterToleranceRmsAbs), fmt5(nTotal0));

        boolean converged = false;
        double lastMaxResidual = Double.POSITIVE_INFINITY;
        int performedIterations = 0;

        for (int iter = 1; iter <= maxIterations; iter++) {
            performedIterations = iter;

            // 1) BdG 行列（2N×2N）
            DMatrixRMaj bdgMatrix = model.buildBdgMatrix(state);

            // 2) 固有分解
            EigenDecompositionResult eigen = eigenBackend.decomposeSymmetricAndSort(bdgMatrix);

            // 3) 更新候補（mixing 前）
            MeanFieldIncrements inc = updater.computeIncrements(state, eigen);

            // 4) 残差（mixing 前）
            Residual residual = measureResidual(state, inc);

            // 収束判定用に最大残差を取得
            double maxDensityResidual = residual.maxDensity();
            double maxOrderParameterResidual = residual.getMaxOrderParameter();
            double rmsOrderParameterResidual = residual.getRmsOrderParameter();

            // 収束判定用に最後の最大残差を保存
            lastMaxResidual = Math.max(maxDensityResidual,
                    Math.max(maxOrderParameterResidual, rmsOrderParameterResidual));

            // 5) mixing して state を更新（毎回必ず）
            applyMixing(state, inc, densityAlpha, orderParameterAlpha);

            // 6) 収束判定（mixing 前 residual で判定）
            boolean densityOk = maxDensityResidual < densityToleranceAbs;
            // Δは RMS を主判定、max は暴れ検知（偽収束排除のため両方要求）
            boolean opMaxOk = maxOrderParameterResidual < orderParameterToleranceAbs;
            boolean opRmsOk = rmsOrderParameterResidual < orderParameterToleranceRmsAbs;

            // Δは RMS を主判定、max は暴れ検知（偽収束排除のため両方要求）
            boolean orderParameterOk = opMaxOk && opRmsOk;

            double nTotal = totalParticleNumber(state);

            log.info(
                    "SCF反復 {} / {}：残差（絶対値）の最大（密度={}／許容={} [{}]、秩序パラメータmax={}／許容={} [{}]、"
                            + "秩序パラメータRMS={}／許容={} [{}]）、N={}",
                    iter, maxIterations, fmt5(maxDensityResidual), fmt5(densityToleranceAbs),
                    densityOk ? "OK" : "NG", fmt5(maxOrderParameterResidual),
                    fmt5(orderParameterToleranceAbs), opMaxOk ? "OK" : "NG",
                    fmt5(rmsOrderParameterResidual), fmt5(orderParameterToleranceRmsAbs),
                    opRmsOk ? "OK" : "NG", fmt5(nTotal));

            if (densityOk && orderParameterOk) {
                converged = true;
                // 収束到達ログ（INFO）
                log.info("SCFが収束条件を満たしました。反復回数={}（密度残差={}、秩序パラメータ残差 max={}、RMS={}）", iter,
                        fmt5(maxDensityResidual), fmt5(maxOrderParameterResidual),
                        fmt5(rmsOrderParameterResidual));
                break;
            }
        }

        long elapsedMs = (System.nanoTime() - t0) / 1_000_000L;

        double nTotalFinal = totalParticleNumber(state);

        if (converged) {
            log.info("SCFが収束しました。反復回数={}、所要時間={}ms（μ1={} μ2={}、N={}）", performedIterations,
                    elapsedMs, fmt5(state.getMu1()), fmt5(state.getMu2()), fmt5(nTotalFinal));
        } else {
            log.warn("SCFが未収束で終了しました。反復回数={}、所要時間={}ms（μ1={} μ2={}、N={}）", performedIterations,
                    elapsedMs, fmt5(state.getMu1()), fmt5(state.getMu2()), fmt5(nTotalFinal));
        }

        return new SolveResult(converged, performedIterations, lastMaxResidual);
    }

    /**
     * 実行結果を表すクラスです。
     */
    @Value
    public static class SolveResult {

        /**
         * 収束したかどうかです。
         */
        boolean converged;

        /**
         * 実行した反復回数です。
         */
        int iterations;

        /**
         * 最終反復での最大変化量です。
         */
        double lastMaxChange;
    }

    /**
     * state の配列長がモデルのサイト数と一致することを検証します。
     *
     * @param state 状態です
     * @param siteCount サイト数です
     * @throws IllegalArgumentException サイズ不整合の場合に発生します
     */
    private static void ensureStateSizeMatchesModel(MeanFieldState state, int siteCount) {
        if (state.getSiteCount() != siteCount) {
            throw new IllegalArgumentException("state.getSiteCount と model.getSiteCount が一致しません: "
                    + state.getSiteCount() + " vs " + siteCount);
        }
        if (state.getDelta().length != siteCount || state.getN1().length != siteCount
                || state.getN2().length != siteCount) {
            throw new IllegalArgumentException("state の配列長が siteCount と一致しません");
        }
    }

    /**
     * mixing を適用して state を更新します。
     *
     * 密度（n1, n2）と秩序パラメータ（Δ）で mixing 係数を分けます。 また、秩序パラメータは全体符号（Δ と -Δ の等価性）による揺れを抑えるため、
     * 旧値との内積で符号合わせしてから mixing します。
     *
     * @param state 更新対象です
     * @param inc 更新候補値です
     * @param densityAlpha 密度（n1, n2）に適用する mixing 係数です（0 より大きく 1 以下）
     * @param orderParameterAlpha 秩序パラメータ（Δ）に適用する mixing 係数です（0 より大きく 1 以下）
     */
    private static void applyMixing(MeanFieldState state, MeanFieldIncrements inc,
            double densityAlpha, double orderParameterAlpha) {

        // 密度 n1
        double[] n1Old = state.getN1();
        // 密度 n2
        double[] n2Old = state.getN2();
        // Δ（秩序パラメータ）
        double[] opOld = state.getDelta();

        double[] n1Calc = inc.getN1Calculated();
        double[] n2Calc = inc.getN2Calculated();
        double[] opCalc = inc.getDeltaCalculated();

        int siteCount = state.getSiteCount();

        // Δ と -Δ は等価になり得るため、全体符号の揺れを抑える（実数Δ前提）
        double sign = orderParameterSignFactor(opOld, opCalc);

        for (int i = 0; i < siteCount; i++) {
            n1Old[i] = densityAlpha * n1Calc[i] + (1.0 - densityAlpha) * n1Old[i];
            n2Old[i] = densityAlpha * n2Calc[i] + (1.0 - densityAlpha) * n2Old[i];

            double opCalcAligned = sign * opCalc[i];
            opOld[i] = orderParameterAlpha * opCalcAligned + (1.0 - orderParameterAlpha) * opOld[i];
        }
    }

    /**
     * 秩序パラメータ（Δ）の全体符号を、旧値との整合が取れるように決めます（実数Δ前提）。
     *
     * Δ と -Δ は物理的に等価になり得るため、反復で全体符号が反転すると |Δ_calc - Δ_old| が人工的に大きくなり、収束判定や mixing が不安定になります。
     *
     * @param oldValues 旧値（現在の state）
     * @param newValues 新値（更新候補）
     * @return newValues に掛ける符号（+1 または -1）
     */
    private static double orderParameterSignFactor(double[] oldValues, double[] newValues) {
        double dot = 0.0;
        for (int i = 0; i < oldValues.length; i++) {
            dot += oldValues[i] * newValues[i];
        }
        return (dot < 0.0) ? -1.0 : 1.0;
    }

    /**
     * 数値を小数点以下5桁までの文字列に整形します。
     *
     * @param v 数値です
     * @return 整形した文字列です
     */
    private static String fmt5(double v) {
        return String.format(Locale.ROOT, "%.5f", v);
    }

    /**
     * mixing 前（= 自己無撞着残差）の最大値を計測します。
     *
     * <p>
     * ここでの residual は「計算された更新候補値」と「現在の state」の不一致です。 mixing 係数 α に依存しない収束判定ができます。
     * </p>
     *
     * @param state 現在の平均場状態
     * @param inc 更新候補値（mixing 前）
     * @return 残差メトリクス
     */
    private static Residual measureResidual(MeanFieldState state, MeanFieldIncrements inc) {
        int siteCount = state.getSiteCount();

        double[] n1Old = state.getN1();
        double[] n2Old = state.getN2();
        double[] opOld = state.getDelta();

        double[] n1Calc = inc.getN1Calculated();
        double[] n2Calc = inc.getN2Calculated();
        double[] opCalc = inc.getDeltaCalculated();

        double maxN1 = 0.0;
        double maxN2 = 0.0;

        // Δ残差：max は全サイト、RMS は密度重み付き
        double maxOrderParameter = 0.0;

        double sumWeightedSq = 0.0; // Σ w_i * diff^2
        double sumWeight = 0.0; // Σ w_i

        // Δ と -Δ は等価になり得るため、全体符号の揺れを抑える（実数Δ前提）
        double sign = orderParameterSignFactor(opOld, opCalc);

        for (int i = 0; i < siteCount; i++) {
            // --- 密度（n1, n2）の残差（max） ---
            double dn1 = Math.abs(n1Calc[i] - n1Old[i]);
            double dn2 = Math.abs(n2Calc[i] - n2Old[i]);
            if (dn1 > maxN1) {
                maxN1 = dn1;
            }
            if (dn2 > maxN2) {
                maxN2 = dn2;
            }

            // --- 秩序パラメータ（Δ）の残差 ---
            // 全体符号を旧値に合わせてから差分を取る
            double opCalcAligned = sign * opCalc[i];
            double diff = opCalcAligned - opOld[i];

            // max |Δ_calc - Δ_old| は全サイトで取る
            double ad = Math.abs(diff);
            if (ad > maxOrderParameter) {
                maxOrderParameter = ad;
            }

            // RMS は密度重み付き（真空領域は自然に寄与が小さくなる）
            double w = n1Old[i] + n2Old[i]; // w_i = n_tot(i)
            if (w > 0.0) {
                sumWeightedSq += w * diff * diff;
                sumWeight += w;
            }
        }

        // 全サイトが真空に近いなどで sumWeight=0 の場合は 0 とする（この状況でΔ収束判定は別途 max が効く）
        double rmsOrderParameter = (sumWeight > 0.0) ? Math.sqrt(sumWeightedSq / sumWeight) : 0.0;

        return Residual.of(maxN1, maxN2, maxOrderParameter, rmsOrderParameter);
    }

    /**
     * 残差の内訳（mixing 前）です。
     *
     * residual は「更新候補（計算値）」と「現在の state（旧値）」の不一致を表します。<br>
     * mixing 係数に依存しないため、自己無撞着の収束判定に使えます。
     */
    @Value
    private static class Residual {

        /**
         * 各サイト i について {@code |n1_calc[i] - n1_old[i]|} を計算し、その最大値を取ります。
         */
        double maxN1;

        /**
         * 各サイト i について {@code |n2_calc[i] - n2_old[i]|} を計算し、その最大値を取ります。
         */
        double maxN2;

        /**
         * 各サイト i について {@code |OP_calc[i] - OP_old[i]|} を計算し、その最大値を取ります。<br>
         * ここで OP は秩序パラメータ（実装上はΔ）を表します。
         */
        double maxOrderParameter;

        /**
         * 各サイト i について {@code (OP_calc[i] - OP_old[i])^2} を計算し、その平均の平方根（RMS）を取ります。<br>
         * ここで OP は秩序パラメータ（実装上は Δ）を表します。
         *
         * <p>
         * 注意：これは OP（Δ）そのものの平均ではなく、「差分」のRMSです。<br>
         * そのため、FFLOのように空間的に符号反転する場合でも相殺でゼロになりにくく、収束判定として安定します。
         * </p>
         */
        double rmsOrderParameter;

        /**
         * 残差オブジェクトを生成します。
         *
         * @param maxN1 n1 の最大残差です
         * @param maxN2 n2 の最大残差です
         * @param maxOrderParameter 秩序パラメータの最大残差です
         * @param rmsOrderParameter 秩序パラメータの残差RMSです
         * @return 残差オブジェクトです
         */
        static Residual of(double maxN1, double maxN2, double maxOrderParameter,
                double rmsOrderParameter) {
            return new Residual(maxN1, maxN2, maxOrderParameter, rmsOrderParameter);
        }

        /**
         * 密度（n1, n2）の残差（絶対値）の最大値を返します。
         *
         * @return 密度残差の最大値です
         */
        double maxDensity() {
            return Math.max(maxN1, maxN2);
        }
    }

    /**
     * 現在の state から全粒子数 N (= Σ_i (n1_i + n2_i)) を計算します。
     *
     * @param state 平均場状態です
     * @return 全粒子数です
     */
    private static double totalParticleNumber(MeanFieldState state) {
        double sum = 0.0;
        double[] n1 = state.getN1();
        double[] n2 = state.getN2();
        for (int i = 0; i < state.getSiteCount(); i++) {
            sum += n1[i] + n2[i];
        }
        return sum;
    }
}
