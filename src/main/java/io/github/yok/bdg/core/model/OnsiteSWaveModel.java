package io.github.yok.bdg.core.model;

import io.github.yok.bdg.core.lattice.Lattice;
import io.github.yok.bdg.core.state.MeanFieldState;
import lombok.Getter;
import org.ejml.data.DMatrixRMaj;

/**
 * オンサイト引力（U&lt;0）とオンサイト s-wave（Δ_i をサイトごとに持つ）を仮定した BdG モデルです。
 *
 * <p>
 * 有効ポテンシャルは W_{iσ} = V_i + U n_{i\barσ} - μ_σ とし、 BdG 行列は [ H1 Δ ; Δ -H2 ] を構築します（Δ は実数を想定します）。
 * </p>
 */
@Getter
public final class OnsiteSWaveModel implements BdgModel {

    /**
     * 格子です。
     */
    private final Lattice lattice;

    /**
     * 外部ポテンシャルです。
     */
    private final Potential potential;

    /**
     * ホッピング t です。
     */
    private final double hoppingT;

    /**
     * オンサイト相互作用 U です。
     */
    private final double interactionU;

    /**
     * モデルを生成します。
     *
     * @param lattice 格子です（null 不可）
     * @param potential 外部ポテンシャルです（null 不可）
     * @param hoppingT ホッピング t です
     * @param interactionU オンサイト相互作用 U です
     * @throws IllegalArgumentException lattice または potential が null の場合に発生します
     */
    public OnsiteSWaveModel(Lattice lattice, Potential potential, double hoppingT,
            double interactionU) {
        if (lattice == null) {
            throw new IllegalArgumentException("lattice は null 不可です");
        }
        if (potential == null) {
            throw new IllegalArgumentException("potential は null 不可です");
        }
        this.lattice = lattice;
        this.potential = potential;
        this.hoppingT = hoppingT;
        this.interactionU = interactionU;
    }

    /**
     * サイト数 N を返します。
     *
     * @return サイト数です
     */
    @Override
    public int siteCount() {
        return lattice.siteCount();
    }

    /**
     * BdG 行列を実対称として扱えるかを返します。
     *
     * @return 実対称として扱えるため true です
     */
    @Override
    public boolean isRealSymmetric() {
        return true;
    }

    /**
     * 2N×2N の BdG 行列を構築して返します。
     *
     * @param state 平均場状態です
     * @return BdG 行列です
     */
    @Override
    public DMatrixRMaj buildBdgMatrix(MeanFieldState state) {

        int siteCount = lattice.siteCount();
        DMatrixRMaj bdgMatrix = new DMatrixRMaj(2 * siteCount, 2 * siteCount);

        double[] density1 = state.getN1();
        double[] density2 = state.getN2();
        double[] delta = state.getDelta();

        double mu1 = state.getMu1();
        double mu2 = state.getMu2();

        for (int site = 0; site < siteCount; site++) {
            double externalV = potential.valueAt(site);

            double effectiveW1 = externalV + interactionU * density2[site] - mu1;
            double effectiveW2 = externalV + interactionU * density1[site] - mu2;

            int uIndex = site;
            int vIndex = site + siteCount;

            bdgMatrix.set(uIndex, uIndex, effectiveW1);
            bdgMatrix.set(vIndex, vIndex, -effectiveW2);

            double deltaSite = delta[site];
            bdgMatrix.set(uIndex, vIndex, deltaSite);
            bdgMatrix.set(vIndex, uIndex, deltaSite);

            for (int neighbor : lattice.neighborsOf(site)) {
                // H1 の最近接ホッピング（-t）
                bdgMatrix.add(uIndex, neighbor, -hoppingT);

                // -H2 の最近接ホッピング（+t）
                bdgMatrix.add(vIndex, neighbor + siteCount, +hoppingT);
            }
        }

        return bdgMatrix;
    }
}
