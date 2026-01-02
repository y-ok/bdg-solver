package io.github.yok.bdg.core.state;

import io.github.yok.bdg.app.BdgProperties;
import io.github.yok.bdg.core.lattice.SquareLattice2D;
import io.github.yok.bdg.core.model.BdgModel;
import io.github.yok.bdg.core.model.ChemicalPotentialInitialEstimator;
import io.github.yok.bdg.core.model.ToroidalTrapPotential;

/**
 * トロイダルトラップを想定した平均場初期値生成です。
 *
 * <ul>
 * <li>密度: 外部ポテンシャル V(i) と μ(seed) から作る LDA 風（フェルミ分布形）</li>
 * <li>秩序パラメータ(Δ): 密度プロファイルに追従させて与える（= V(i) にのみ依存して自然に中心/リングが決まる）</li>
 * <li>μ1, μ2: 探索用 seed として muInitial を両方に設定（h は後段で反映する想定）</li>
 * </ul>
 */
public final class ToroidalLdaRingInitializer implements MeanFieldStateInitializer {

    /** 秩序パラメータ(Δ) の初期振幅の既定値です。 */
    private static final double DEFAULT_ORDER_PARAMETER0 = 0.1;

    /** 初期化用途として、温度の下限を設けます。 */
    private static final double MIN_T_EFF = 1e-3;

    private final BdgProperties properties;
    private final BdgModel model;
    private final ToroidalTrapPotential trapPotential;

    /**
     * 初期値生成器を作成します。
     *
     * @param properties 設定です
     * @param model BdG モデルです
     * @param trapPotential トラップポテンシャルです
     * @throws IllegalArgumentException 引数が null の場合に発生します
     */
    public ToroidalLdaRingInitializer(BdgProperties properties, BdgModel model,
            ToroidalTrapPotential trapPotential) {
        if (properties == null) {
            throw new IllegalArgumentException("properties は null 不可です");
        }
        if (model == null) {
            throw new IllegalArgumentException("model は null 不可です");
        }
        if (trapPotential == null) {
            throw new IllegalArgumentException("trapPotential は null 不可です");
        }
        this.properties = properties;
        this.model = model;
        this.trapPotential = trapPotential;
    }

    /**
     * 平均場初期状態を作成します。
     *
     * <p>
     * まず密度を LDA 風に作り、（N 指定があれば）全体スケールして N を合わせます。 その後、秩序パラメータは密度に追従させて与えます（ポテンシャル形状に自然に依存）。
     * </p>
     *
     * @return 初期状態です
     */
    @Override
    public MeanFieldState create() {
        int siteCount = model.siteCount();
        MeanFieldState state = new MeanFieldState(siteCount);

        // μ探索の初期値（seed）
        double muInitial = resolveMuInitial();
        state.setMu1(muInitial);
        state.setMu2(muInitial);

        // LDA 風密度生成に使う温度（下限あり）
        double tEff = Math.max(properties.getModel().getTemperature(), MIN_T_EFF);
        double beta = 1.0 / tEff;

        double[] n1 = state.getN1();
        double[] n2 = state.getN2();
        double[] op = state.getDelta();

        // 1) 密度（LDA風）
        for (int i = 0; i < siteCount; i++) {
            double v = trapPotential.valueAt(i);

            // n ≈ 1/(1+exp(beta*(V-μ))) を初期値として使用（0〜1 にクランプ）
            double nLda = fermiLikeOccupation(beta * (v - muInitial));
            n1[i] = clamp01(nLda);
            n2[i] = clamp01(nLda);
        }

        // 2) 総粒子数 N が指定されている場合は、初期密度を N に合わせて全体スケール
        int totalN = properties.getModel().getTotalParticleNumber();
        if (totalN > 0) {
            scaleDensityToTargetTotalN(n1, n2, totalN, siteCount);
        }

        // 3) 秩序パラメータ（Δ）は密度に追従させる（= V(i) にのみ依存）
        // ※ここでは両成分の平均密度を使う（初期偏極は入れない前提に揃える）
        for (int i = 0; i < siteCount; i++) {
            double nAvg = 0.5 * (n1[i] + n2[i]);
            op[i] = DEFAULT_ORDER_PARAMETER0 * nAvg;
        }

        return state;
    }

    /**
     * フェルミ分布形の占有率 1/(1+exp(z)) を返します（オーバーフロー対策あり）。
     *
     * @param z z=beta*(V-μ) です
     * @return 占有率です
     */
    private static double fermiLikeOccupation(double z) {
        if (z > 50.0) {
            return 0.0;
        }
        if (z < -50.0) {
            return 1.0;
        }
        return 1.0 / (1.0 + Math.exp(z));
    }

    /**
     * 初期密度 (n1+n2) の平均が、目標総粒子数 N に一致するように全体スケールします。
     *
     * @param n1 成分1密度です
     * @param n2 成分2密度です
     * @param totalN 目標総粒子数です
     * @param siteCount サイト数です
     */
    private static void scaleDensityToTargetTotalN(double[] n1, double[] n2, int totalN,
            int siteCount) {
        double targetTotalPerSite = totalN / (double) siteCount;

        double sum = 0.0;
        for (int i = 0; i < siteCount; i++) {
            sum += (n1[i] + n2[i]);
        }
        double currentTotalPerSite = sum / siteCount;

        if (currentTotalPerSite <= 0.0) {
            // 全ゼロでスケール不能 → 目標Nに一致する一様密度で埋める
            double nPerComponent = 0.5 * targetTotalPerSite;
            for (int i = 0; i < siteCount; i++) {
                n1[i] = clamp01(nPerComponent);
                n2[i] = clamp01(nPerComponent);
            }
            return;
        }

        double scale = targetTotalPerSite / currentTotalPerSite;
        for (int i = 0; i < siteCount; i++) {
            n1[i] = clamp01(n1[i] * scale);
            n2[i] = clamp01(n2[i] * scale);
        }
    }

    /**
     * 値を 0.0〜1.0 の範囲にクランプします。
     *
     * @param v 値です
     * @return クランプ後の値です
     */
    private static double clamp01(double v) {
        if (v < 0.0) {
            return 0.0;
        }
        if (v > 1.0) {
            return 1.0;
        }
        return v;
    }

    /**
     * mu-initial（平均 μ の seed）を決定します。
     *
     * <p>
     * 設定で自動推定が有効な場合は、対角化なしの推定器で mu-initial を見積もります。 無効な場合は設定値 {@code bdg.model.mu-initial}
     * をそのまま返します。
     * </p>
     *
     * @return mu-initial です
     */
    private double resolveMuInitial() {
        BdgProperties.Model modelProps = properties.getModel();
        BdgProperties.Model.MuInitialEstimation est = modelProps.getMuInitialEstimation();

        // 推定無効の場合、設定値を使用
        if (est == null || !est.isEnabled()) {
            return modelProps.getMuInitial();
        }

        double targetTotalN = modelProps.getTotalParticleNumber();

        // 推定に必要な格子は、ポテンシャルが保持しているものを使う
        SquareLattice2D lattice2D = trapPotential.getLattice2D();

        ChemicalPotentialInitialEstimator estimator = new ChemicalPotentialInitialEstimator(
                lattice2D, trapPotential, modelProps.getT(), modelProps.getU(),
                modelProps.getTemperature(), est.getMaxIterations(), est.getTolerance(),
                est.getBracketMargin(), est.isHartreeEnabled(), est.getHartreeMaxIterations(),
                est.getHartreeTolerance(), est.getHartreeMixing());

        return estimator.estimateMuInitial(targetTotalN);
    }
}
