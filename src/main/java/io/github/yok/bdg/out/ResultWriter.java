package io.github.yok.bdg.out;

import io.github.yok.bdg.core.solver.ChemicalPotentialOptimizer;
import io.github.yok.bdg.core.state.MeanFieldState;

/**
 * 計算結果を出力する処理のインタフェースです。
 *
 * <p>
 * Zeeman 場 h をスキャンして計算することを前提とし、<br>
 * 出力の命名規約に必要な {@code N} と {@code h}、および結果として得られた {@code P_out} を受け取ります。
 * </p>
 */
public interface ResultWriter {

    /**
     * 最適化結果と平均場状態を出力します。
     *
     * @param state 平均場状態です
     * @param optimization 化学ポテンシャル探索結果です
     * @param totalParticleNumber 目標の総粒子数 N（N1+N2）です
     * @param zeemanField Zeeman 場 h（有効磁場、h=(mu1-mu2)/2）です
     * @param outputPolarization 結果として得られた偏極率 P_out（(N1-N2)/(N1+N2)）です
     */
    void write(MeanFieldState state, ChemicalPotentialOptimizer.OptimizationResult optimization,
            int totalParticleNumber, double zeemanField, double outputPolarization);
}
