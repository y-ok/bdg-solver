package io.github.yok.bdg.core.model;

import io.github.yok.bdg.core.state.MeanFieldState;
import org.ejml.data.DMatrixRMaj;

/**
 * 平均場状態から BdG 行列を構築するモデルを表すインタフェースです。
 *
 * <p>
 * 系（相互作用・秩序パラメータ・格子・外部ポテンシャル）を差し替えるための境界です。
 * </p>
 */
public interface BdgModel {

    /**
     * サイト数 N を返します。
     *
     * @return サイト数です
     */
    int siteCount();

    /**
     * BdG 行列を実対称として扱えるかを返します。
     *
     * @return 実対称として扱える場合は true です
     */
    boolean isRealSymmetric();

    /**
     * 2N×2N の BdG 行列を構築して返します。
     *
     * @param state 平均場状態です
     * @return BdG 行列です
     */
    DMatrixRMaj buildBdgMatrix(MeanFieldState state);
}
