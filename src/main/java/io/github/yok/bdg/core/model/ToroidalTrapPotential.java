package io.github.yok.bdg.core.model;

import io.github.yok.bdg.core.lattice.SquareLattice2D;
import lombok.Getter;

/**
 * トロイダルトラップの外部ポテンシャルを計算するクラスです。
 *
 * <p>
 * V(r) = 1/2 * ω_ho * (r/r0)^2 + ω_tr * exp(-r/ξ)
 * </p>
 */
@Getter
public final class ToroidalTrapPotential implements Potential {

    /**
     * 2D 格子です。
     */
    private final SquareLattice2D lattice2D;

    /**
     * 半径スケール r0 です。
     */
    private final double r0;

    /**
     * 調和トラップ係数 ω_ho です。
     */
    private final double omegaHo;

    /**
     * トロイダル項係数 ω_tr です。
     */
    private final double omegaTr;

    /**
     * 減衰長 ξ です。
     */
    private final double xi;

    /**
     * 中心 x 座標です。
     */
    private final double centerX;

    /**
     * 中心 y 座標です。
     */
    private final double centerY;

    /**
     * トロイダルトラップを生成します。
     *
     * @param lattice2D 2D 格子です（null 不可）
     * @param r0 半径スケール r0 です（正）
     * @param omegaHo 調和トラップ係数 ω_ho です
     * @param omegaTr トロイダル項係数 ω_tr です
     * @param xi 減衰長 ξ です（正）
     * @throws IllegalArgumentException lattice2D が null、または r0/xi が正でない場合に発生します
     */
    public ToroidalTrapPotential(SquareLattice2D lattice2D, double r0, double omegaHo,
            double omegaTr, double xi) {
        if (lattice2D == null) {
            throw new IllegalArgumentException("lattice2D は null 不可です");
        }
        if (r0 <= 0) {
            throw new IllegalArgumentException("r0 は正である必要があります: " + r0);
        }
        if (xi <= 0) {
            throw new IllegalArgumentException("xi は正である必要があります: " + xi);
        }
        this.lattice2D = lattice2D;
        this.r0 = r0;
        this.omegaHo = omegaHo;
        this.omegaTr = omegaTr;
        this.xi = xi;

        // 0 始まり格子の中心（偶数サイズなら .5 になる）
        this.centerX = (lattice2D.lx() - 1) / 2.0;
        this.centerY = (lattice2D.ly() - 1) / 2.0;
    }

    /**
     * 指定サイトの外部ポテンシャルを返します。
     *
     * @param siteIndex サイトインデックスです
     * @return 外部ポテンシャルです
     */
    @Override
    public double valueAt(int siteIndex) {
        double dx = lattice2D.xOf(siteIndex) - centerX;
        double dy = lattice2D.yOf(siteIndex) - centerY;
        double r = Math.hypot(dx, dy);

        double rr0 = r / r0;
        double harmonic = 0.5 * omegaHo * rr0 * rr0;
        double toroidal = omegaTr * Math.exp(-r / xi);
        return harmonic + toroidal;
    }
}
