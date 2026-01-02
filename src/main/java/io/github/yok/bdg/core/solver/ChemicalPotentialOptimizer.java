package io.github.yok.bdg.core.solver;

import io.github.yok.bdg.core.state.MeanFieldState;
import lombok.Value;

/**
 * 固定した Zeeman 場 h のもとで、総粒子数条件を満たすように 平均化学ポテンシャル μ を探索するインタフェースです。
 *
 * <p>
 * 各評価では {@code mu1 = mu + h}, {@code mu2 = mu - h} を用います。
 * </p>
 */
public interface ChemicalPotentialOptimizer {

    /**
     * 固定した h のもとで、総粒子数 N が目標値に一致するように μ を探索します。
     *
     * @param initialState 初期状態（SCFの初期値）です
     * @param totalParticleNumber 目標の総粒子数 N (=N1+N2) です
     * @param zeemanField Zeeman 場 h（有効磁場、h=(mu1-mu2)/2）です
     * @return 最適化結果です
     */
    OptimizationResult optimize(MeanFieldState initialState, int totalParticleNumber,
            double zeemanField);

    /**
     * 最適化結果を表すクラスです。
     */
    @Value
    class OptimizationResult {

        /**
         * 最適化後の μ1 です（通常は mu + h）。
         */
        double mu1;

        /**
         * 最適化後の μ2 です（通常は mu - h）。
         */
        double mu2;

        /**
         * 最終的な平均場状態です。
         */
        MeanFieldState state;

        /**
         * 実行した μ 探索の反復回数です。
         */
        int iterations;
    }
}
