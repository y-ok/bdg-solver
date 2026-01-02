package io.github.yok.bdg.core.model;

import io.github.yok.bdg.core.linearalgebra.EigenDecompositionBackend.EigenDecompositionResult;
import io.github.yok.bdg.core.state.MeanFieldIncrements;
import io.github.yok.bdg.core.state.MeanFieldState;

/**
 * 固有分解結果から、平均場（n1, n2, Δ）の更新候補値を計算するインタフェースです。
 *
 * <p>
 * 秩序パラメータの定義や更新式の違い（例: s-wave / d-wave / FFLO）を差し替えるためのインタフェースです。
 * </p>
 */
public interface MeanFieldUpdater {

    /**
     * 更新候補値（mixing 前）を計算して返します。
     *
     * @param state 現在の平均場状態です
     * @param eigen 固有分解結果です
     * @return 更新候補値です
     */
    MeanFieldIncrements computeIncrements(MeanFieldState state, EigenDecompositionResult eigen);
}
