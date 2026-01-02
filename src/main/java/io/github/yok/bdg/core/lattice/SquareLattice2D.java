package io.github.yok.bdg.core.lattice;

import io.github.yok.bdg.app.BdgProperties;

/**
 * 2 次元の正方格子（4 近傍）を表すクラスです。
 *
 * <p>
 * インデックス変換は {@code index = y * lx + x} です。
 * </p>
 */
public final class SquareLattice2D implements Lattice {

    /**
     * x 方向サイズです。
     */
    private final int lx;

    /**
     * y 方向サイズです。
     */
    private final int ly;

    /**
     * 境界条件です。
     */
    private final BdgProperties.Lattice.Boundary boundary;

    /**
     * 2 次元正方格子を生成します。
     *
     * @param lx x 方向サイズです（1 以上）
     * @param ly y 方向サイズです（1 以上）
     * @param boundary 境界条件です（null 不可）
     * @throws IllegalArgumentException lx/ly が 1 未満、または boundary が null の場合に発生します
     */
    public SquareLattice2D(int lx, int ly, BdgProperties.Lattice.Boundary boundary) {
        if (lx <= 0 || ly <= 0) {
            throw new IllegalArgumentException("lx/ly は 1 以上が必要です: " + lx + "x" + ly);
        }
        if (boundary == null) {
            throw new IllegalArgumentException("boundary は null 不可です");
        }
        this.lx = lx;
        this.ly = ly;
        this.boundary = boundary;
    }

    /**
     * x 方向サイズ（Lx）を返します。
     *
     * @return x 方向サイズです
     */
    public int lx() {
        return lx;
    }

    /**
     * y 方向サイズ（Ly）を返します。
     *
     * @return y 方向サイズです
     */
    public int ly() {
        return ly;
    }

    /**
     * 座標 (x, y) をサイトインデックスに変換します。
     *
     * @param x x 座標です（0 以上 lx 未満）
     * @param y y 座標です（0 以上 ly 未満）
     * @return サイトインデックスです
     */
    public int indexOf(int x, int y) {
        return y * lx + x;
    }

    /**
     * サイトインデックスから x 座標を返します。
     *
     * @param siteIndex サイトインデックスです（0 以上 N 未満）
     * @return x 座標です
     */
    public int xOf(int siteIndex) {
        return siteIndex % lx;
    }

    /**
     * サイトインデックスから y 座標を返します。
     *
     * @param siteIndex サイトインデックスです（0 以上 N 未満）
     * @return y 座標です
     */
    public int yOf(int siteIndex) {
        return siteIndex / lx;
    }

    /**
     * サイト数 N を返します。
     *
     * @return サイト数です
     */
    @Override
    public int siteCount() {
        return lx * ly;
    }

    /**
     * 4 近傍の隣接サイトインデックスを返します。
     *
     * @param siteIndex サイトインデックスです（0 以上 N 未満）
     * @return 隣接サイトインデックス配列です（null にはしません）
     */
    @Override
    public int[] neighborsOf(int siteIndex) {
        int x = xOf(siteIndex);
        int y = yOf(siteIndex);

        int[] buf = new int[4];
        int count = 0;

        // 左右隣接サイト
        if (boundary == BdgProperties.Lattice.Boundary.PERIODIC) {
            buf[count++] = indexOf((x - 1 + lx) % lx, y);
            buf[count++] = indexOf((x + 1) % lx, y);
        } else {
            if (x > 0)
                buf[count++] = indexOf(x - 1, y);
            if (x < lx - 1)
                buf[count++] = indexOf(x + 1, y);
        }

        // 上下隣接サイト
        if (boundary == BdgProperties.Lattice.Boundary.PERIODIC) {
            buf[count++] = indexOf(x, (y - 1 + ly) % ly);
            buf[count++] = indexOf(x, (y + 1) % ly);
        } else {
            if (y > 0)
                buf[count++] = indexOf(x, y - 1);
            if (y < ly - 1)
                buf[count++] = indexOf(x, y + 1);
        }

        int[] out = new int[count];
        System.arraycopy(buf, 0, out, 0, count);
        return out;
    }
}
