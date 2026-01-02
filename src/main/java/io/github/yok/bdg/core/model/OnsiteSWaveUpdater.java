package io.github.yok.bdg.core.model;

import io.github.yok.bdg.core.linearalgebra.EigenDecompositionBackend.EigenDecompositionResult;
import io.github.yok.bdg.core.state.MeanFieldIncrements;
import io.github.yok.bdg.core.state.MeanFieldState;
import lombok.Getter;
import org.ejml.data.DMatrixRMaj;

/**
 * オンサイト s-wave の更新式（n1, n2, Δ）を計算するクラスです。
 *
 * <p>
 * 固有分解結果（固有値・固有ベクトル）から、mixing 前の更新候補値を返します。
 * </p>
 */
@Getter
public final class OnsiteSWaveUpdater implements MeanFieldUpdater {

    /**
     * 相互作用 U です。
     */
    private final double interactionU;

    /**
     * 温度 T です。
     */
    private final double temperature;

    /**
     * ゼロ近傍の固有値を除外するための閾値です。
     */
    private final double positiveEnergyEpsilon;

    /**
     * 更新式を生成します。
     *
     * @param interactionU 相互作用 U です
     * @param temperature 温度 T です（0 以上）
     */
    public OnsiteSWaveUpdater(double interactionU, double temperature) {
        this(interactionU, temperature, 1e-12);
    }

    /**
     * 更新式を生成します。
     *
     * @param interactionU 相互作用 U です
     * @param temperature 温度 T です（0 以上）
     * @param positiveEnergyEpsilon ゼロ近傍の固有値を除外するための閾値です（0 以上）
     * @throws IllegalArgumentException temperature または positiveEnergyEpsilon が負の場合に発生します
     */
    public OnsiteSWaveUpdater(double interactionU, double temperature,
            double positiveEnergyEpsilon) {
        if (temperature < 0) {
            throw new IllegalArgumentException("temperature は 0 以上が必要です: " + temperature);
        }
        if (positiveEnergyEpsilon < 0) {
            throw new IllegalArgumentException(
                    "positiveEnergyEpsilon は 0 以上が必要です: " + positiveEnergyEpsilon);
        }
        this.interactionU = interactionU;
        this.temperature = temperature;
        this.positiveEnergyEpsilon = positiveEnergyEpsilon;
    }

    /**
     * 固有分解結果から更新候補値（mixing 前）を計算します。
     *
     * <p>
     * 本実装では、以下の形で計算します。
     * </p>
     * <ul>
     * <li>n1(i) = Σ |u_{n,i}|^2 f(E_n)</li>
     * <li>n2(i) = Σ |v_{n,i}|^2 f(-E_n)</li>
     * <li>Δ(i) = U Σ u_{n,i} v_{n,i} f(E_n)</li>
     * </ul>
     *
     * @param state 現在の平均場状態です
     * @param eigen 固有分解結果です
     * @return 更新候補値です
     * @throws IllegalArgumentException state または eigen が null の場合に発生します
     */
    @Override
    public MeanFieldIncrements computeIncrements(MeanFieldState state,
            EigenDecompositionResult eigen) {
        if (state == null) {
            throw new IllegalArgumentException("state は null 不可です");
        }
        if (eigen == null) {
            throw new IllegalArgumentException("eigen は null 不可です");
        }

        int siteCount = state.getSiteCount();

        // 返却する mixing 前の計算結果（サイトごとの配列）です。
        double[] n1Calculated = new double[siteCount];
        double[] n2Calculated = new double[siteCount];
        double[] deltaCalculated = new double[siteCount];

        double[] energies = eigen.getEigenvalues();
        DMatrixRMaj vectors = eigen.getEigenvectors();

        // 固有ベクトル行列の表現は「列 = 固有ベクトル」です。
        // BdG の表記では、上半分が u、下半分が v（それぞれ長さ N）に対応します。
        for (int mode = 0; mode < energies.length; mode++) {
            double energy = energies[mode];

            // ゼロ近傍は数値誤差として除外します。
            if (Math.abs(energy) <= positiveEnergyEpsilon) {
                continue;
            }

            double fE = fermiDistribution(energy, temperature);
            double fMinusE = fermiDistribution(-energy, temperature);

            for (int site = 0; site < siteCount; site++) {
                double u = vectors.get(site, mode); // u_{n,i}
                double v = vectors.get(site + siteCount, mode); // v_{n,i}

                double uu = u * u;
                double vv = v * v;

                // 成分1粒子密度:
                n1Calculated[site] += uu * fE;

                // 成分2粒子密度:
                n2Calculated[site] += vv * fMinusE;

                // ギャップ（秩序パラメータ）:
                deltaCalculated[site] += interactionU * (u * v) * fE;
            }
        }
        return new MeanFieldIncrements(n1Calculated, n2Calculated, deltaCalculated);
    }

    /**
     * フェルミ分布 f(E) を返します。
     *
     * @param energy エネルギー E です
     * @param temperature 温度 T です
     * @return フェルミ分布値です
     */
    private static double fermiDistribution(double energy, double temperature) {
        if (temperature <= 0.0) {
            // T=0 のステップ関数です。
            if (energy > 0.0) {
                return 0.0;
            }
            if (energy < 0.0) {
                return 1.0;
            }
            return 0.5;
        }

        double x = energy / temperature;

        // exp のオーバーフロー/アンダーフロー回避です。
        if (x > 50.0) {
            return 0.0;
        }
        if (x < -50.0) {
            return 1.0;
        }

        return 1.0 / (Math.exp(x) + 1.0);
    }
}
