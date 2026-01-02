package io.github.yok.bdg.core.lattice;

/**
 * BdG 計算で使用する格子（サイト集合）を表すインタフェースです。
 *
 * <p>
 * ソルバ側は格子の具体形状を意識せず、サイト数と隣接情報のみを利用します。
 * </p>
 */
public interface Lattice {

    /**
     * サイト数 N を返します。
     *
     * @return サイト数です
     */
    int siteCount();

    /**
     * 指定サイトに隣接するサイトのインデックス配列を返します。
     *
     * @param siteIndex サイトインデックスです（0 以上 N 未満）
     * @return 隣接サイトのインデックス配列です（null にはしません）
     */
    int[] neighborsOf(int siteIndex);
}
