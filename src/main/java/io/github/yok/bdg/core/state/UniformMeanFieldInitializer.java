package io.github.yok.bdg.core.state;

import io.github.yok.bdg.app.BdgProperties;
import io.github.yok.bdg.core.model.BdgModel;

/**
 * 一様な平均場初期値を生成します（全サイト同一）。
 *
 * <ul>
 * <li>n1, n2: 総粒子数 N とサイト数から一様密度を作り、成分1/2に半分ずつ割り当てます</li>
 * <li>Δ: 定数 delta0 を全サイトに設定します</li>
 * <li>μ1, μ2: 探索用の初期値 muInitial を両方に設定します（Zeeman場 h は後段で反映する想定）</li>
 * </ul>
 */
public final class UniformMeanFieldInitializer implements MeanFieldStateInitializer {

    /**
     * 初期ギャップの既定値です。
     */
    private static final double DEFAULT_DELTA0 = 0.1;

    /**
     * 総粒子数 N が未指定のときの n1 既定値（1サイトあたり）です。
     */
    private static final double DEFAULT_N1 = 0.5;

    /**
     * 総粒子数 N が未指定のときの n2 既定値（1サイトあたり）です。
     */
    private static final double DEFAULT_N2 = 0.5;

    private final BdgProperties properties;
    private final BdgModel model;

    /**
     * 初期値生成器を作成します。
     *
     * @param properties 設定です
     * @param model BdGモデルです
     * @throws IllegalArgumentException 引数が null の場合に発生します
     */
    public UniformMeanFieldInitializer(BdgProperties properties, BdgModel model) {
        if (properties == null) {
            throw new IllegalArgumentException("properties は null 不可です");
        }
        if (model == null) {
            throw new IllegalArgumentException("model は null 不可です");
        }
        this.properties = properties;
        this.model = model;
    }

    /**
     * 一様な平均場初期状態を作成します。
     *
     * @return 初期状態です
     */
    @Override
    public MeanFieldState create() {
        int siteCount = model.siteCount();
        MeanFieldState state = new MeanFieldState(siteCount);

        // μ探索の初期値（seed）
        double muInitial = properties.getModel().getMuInitial();
        state.setMu1(muInitial);
        state.setMu2(muInitial);

        // 総粒子数 N から一様密度を作る
        int totalN = properties.getModel().getTotalParticleNumber();

        double n1;
        double n2;

        if (totalN > 0) {
            // 1サイトあたりの総密度（成分1+成分2）
            double nTotalPerSite = totalN / (double) siteCount;

            // 成分1/2 に半分ずつ割り当て
            n1 = clamp01(0.5 * nTotalPerSite);
            n2 = clamp01(0.5 * nTotalPerSite);
        } else {
            // N を使わない場合のデフォルト（半充填）
            n1 = DEFAULT_N1;
            n2 = DEFAULT_N2;
        }

        double delta0 = DEFAULT_DELTA0;

        fill(state.getN1(), n1);
        fill(state.getN2(), n2);
        fill(state.getDelta(), delta0);

        return state;
    }

    /**
     * 配列を指定値で埋めます。
     *
     * @param a 配列です
     * @param v 設定値です
     */
    private static void fill(double[] a, double v) {
        for (int i = 0; i < a.length; i++) {
            a[i] = v;
        }
    }

    /**
     * 値を 0.0〜1.0 の範囲に丸め込みます。
     *
     * @param v 値です
     * @return 丸め込み後の値です
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
}
