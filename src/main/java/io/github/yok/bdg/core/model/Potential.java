package io.github.yok.bdg.core.model;

/**
 * 外部ポテンシャル V_i を提供するインタフェースです。
 *
 * <p>
 * 格子形状やトラップ形状の差し替えを容易にするための境界です。
 * </p>
 */
public interface Potential {

    /**
     * 指定サイトの外部ポテンシャル V_i を返します。
     *
     * @param siteIndex サイトインデックスです（0 以上 N 未満）
     * @return 外部ポテンシャルです
     */
    double valueAt(int siteIndex);
}
