package io.github.yok.bdg.core.state;

/**
 * 平均場状態の初期値を生成するインターフェースです。
 */
public interface MeanFieldStateInitializer {

    /**
     * 初期化済みの平均場状態を生成します。
     *
     * @return 初期化済みの平均場状態です
     */
    MeanFieldState create();
}
