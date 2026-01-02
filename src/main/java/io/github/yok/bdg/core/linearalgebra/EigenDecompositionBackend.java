package io.github.yok.bdg.core.linearalgebra;

import lombok.Value;
import org.ejml.data.DMatrixRMaj;

/**
 * 固有分解を提供するバックエンドを表すインタフェースです。
 *
 * <p>
 * 使用するライブラリや複素数対応を差し替えやすくするためのインタフェースです。
 * </p>
 */
public interface EigenDecompositionBackend {

    /**
     * 実対称行列を固有分解し、固有値昇順の結果を返します。
     *
     * @param symmetricMatrix 実対称行列です
     * @return 固有値昇順の固有分解結果です
     * @throws IllegalStateException 固有分解に失敗した場合に発生します
     */
    EigenDecompositionResult decomposeSymmetricAndSort(DMatrixRMaj symmetricMatrix);

    /**
     * 固有分解の結果（固有値・固有ベクトル）を保持するクラスです。
     *
     * <p>
     * 固有ベクトル行列は「列が固有ベクトル」である前提です。
     * </p>
     */
    @Value
    class EigenDecompositionResult {

        /**
         * 固有値配列です。
         */
        double[] eigenvalues;

        /**
         * 固有ベクトル行列です（列が固有ベクトルです）。
         */
        DMatrixRMaj eigenvectors;
    }
}
