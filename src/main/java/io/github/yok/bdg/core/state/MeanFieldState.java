package io.github.yok.bdg.core.state;

import lombok.Getter;
import lombok.Setter;

/**
 * 自己無撞着（SCF）ループの平均場状態を保持するクラスです。
 *
 * <p>
 * サイト数 N に対して Δ, n1, n2 を長さ N の配列で保持します。
 * </p>
 */
@Getter
public final class MeanFieldState {

    /**
     * サイト数です。
     */
    private final int siteCount;

    /**
     * 秩序パラメータ Δ_i の配列です。
     */
    private final double[] delta;

    /**
     * 密度 n_{i1} の配列です。
     */
    private final double[] n1;

    /**
     * 密度 n_{i2} の配列です。
     */
    private final double[] n2;

    /**
     * 化学ポテンシャル μ1 です。
     */
    @Setter
    private double mu1;

    /**
     * 化学ポテンシャル μ2 です。
     */
    @Setter
    private double mu2;

    /**
     * 平均場状態を生成します。
     *
     * @param siteCount サイト数です（1 以上）
     * @throws IllegalArgumentException siteCount が 1 未満の場合に発生します
     */
    public MeanFieldState(int siteCount) {
        if (siteCount <= 0) {
            throw new IllegalArgumentException("siteCount は 1 以上が必要です: " + siteCount);
        }
        this.siteCount = siteCount;
        this.delta = new double[siteCount];
        this.n1 = new double[siteCount];
        this.n2 = new double[siteCount];
    }

    /**
     * 成分1の総粒子数 N1 (= Σ n_{i1}) を返します。
     * 
     * @return 成分1の総粒子数です
     */
    public double totalN1() {
        return sum(n1);
    }

    /**
     * 成分2の総粒子数 N2 (= Σ n_{i2}) を返します。
     * 
     * @return 成分2の総粒子数です
     */
    public double totalN2() {
        return sum(n2);
    }

    /**
     * 全粒子数 N (= N1 + N2) を返します。
     * 
     * @return 全粒子数です
     */
    public double totalN() {
        return totalN1() + totalN2();
    }

    /**
     * 配列の総和を返します。
     *
     * @param array 配列です
     * @return 総和です
     */
    private double sum(double[] array) {
        double s = 0.0;
        for (double v : array) {
            s += v;
        }
        return s;
    }
}
