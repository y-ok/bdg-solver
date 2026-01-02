package io.github.yok.bdg.core.linearalgebra;

import java.util.Arrays;
import org.ejml.data.DMatrixRMaj;
import org.ejml.dense.row.factory.DecompositionFactory_DDRM;
import org.ejml.interfaces.decomposition.EigenDecomposition_F64;

/**
 * EJML を用いて、実対称行列の固有分解を行うクラスです。
 *
 * <p>
 * 固有値を昇順に並べ替え、固有ベクトルも同じ順序に揃えて返します。
 * </p>
 */
public final class EjmlSymmetricEigenDecompositionBackend implements EigenDecompositionBackend {

    /**
     * 実対称行列を固有分解し、固有値昇順の結果を返します。
     *
     * @param symmetricMatrix 実対称行列です
     * @return 固有値昇順の固有分解結果です
     * @throws IllegalArgumentException symmetricMatrix が null の場合に発生します
     * @throws IllegalStateException 固有分解に失敗した場合に発生します
     */
    @Override
    public EigenDecompositionResult decomposeSymmetricAndSort(DMatrixRMaj symmetricMatrix) {
        if (symmetricMatrix == null) {
            throw new IllegalArgumentException("symmetricMatrix は null 不可です");
        }

        // 行列次元（2N×2N など）を取得します。
        int dim = symmetricMatrix.numRows;

        // 実対称行列の固有分解器を生成します。
        EigenDecomposition_F64<DMatrixRMaj> decomposition =
                DecompositionFactory_DDRM.eig(dim, true, true);

        // 固有分解を実行します。
        if (!decomposition.decompose(symmetricMatrix)) {
            throw new IllegalStateException("固有分解に失敗しました（EJML）");
        }

        // EJML から固有値・固有ベクトルを取り出して配列/行列にまとめます。
        // ここでは「列 = 固有ベクトル」という表現に正規化します。
        double[] eigenvalues = new double[dim];
        DMatrixRMaj eigenvectors = new DMatrixRMaj(dim, dim);

        for (int col = 0; col < dim; col++) {
            // 実対称を想定しているため、固有値は実数部のみを使います。
            eigenvalues[col] = decomposition.getEigenvalue(col).getReal();

            // 固有ベクトルは dim×1 の列ベクトルで返る想定です。
            DMatrixRMaj vec = decomposition.getEigenVector(col);
            if (vec == null) {
                throw new IllegalStateException("固有ベクトルが取得できません: col=" + col);
            }

            // 返却用の固有ベクトル行列に、col 列としてコピーします。
            for (int row = 0; row < dim; row++) {
                eigenvectors.set(row, col, vec.get(row, 0));
            }
        }

        // 固有値を昇順に並べ替えるため、argsort（並び替えインデックス）を作ります。
        int[] order = argsortAscending(eigenvalues);

        // 固有値を昇順にし、固有ベクトルも同じ順序で並べ替えます。
        double[] sortedValues = new double[dim];
        DMatrixRMaj sortedVectors = new DMatrixRMaj(dim, dim);

        for (int newCol = 0; newCol < dim; newCol++) {
            int oldCol = order[newCol];

            sortedValues[newCol] = eigenvalues[oldCol];

            // oldCol 列を newCol 列へコピーします（固有ベクトルも並び替え）。
            for (int row = 0; row < dim; row++) {
                sortedVectors.set(row, newCol, eigenvectors.get(row, oldCol));
            }
        }

        return new EigenDecompositionResult(sortedValues, sortedVectors);
    }

    /**
     * 配列を昇順ソートしたときのインデックス順（argsort）を返します。
     *
     * @param values 対象配列です
     * @return 昇順のインデックス配列です
     */
    private static int[] argsortAscending(double[] values) {
        // ソート対象は「値」ではなく「インデックス」です。
        Integer[] indices = new Integer[values.length];
        for (int i = 0; i < values.length; i++) {
            indices[i] = i;
        }

        // values の大小でインデックス配列を並べ替えます。
        Arrays.sort(indices, (i, j) -> Double.compare(values[i], values[j]));

        // int[] に詰め替えます。
        int[] order = new int[values.length];
        for (int i = 0; i < values.length; i++) {
            order[i] = indices[i];
        }
        return order;
    }
}
