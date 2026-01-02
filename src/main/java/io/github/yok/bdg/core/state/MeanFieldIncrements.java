package io.github.yok.bdg.core.state;

import lombok.Value;

/**
 * BdG の固有分解から得られる「更新候補値（mixing 前）」を保持するクラスです。
 *
 * <p>
 * SCF ループでは、この値を mixing して MeanFieldState に反映します。
 * </p>
 */
@Value
public class MeanFieldIncrements {

    /**
     * 密度 n1 の計算結果です。
     */
    double[] n1Calculated;

    /**
     * 密度 n2 の計算結果です。
     */
    double[] n2Calculated;

    /**
     * 秩序パラメータ Δ の計算結果です。
     */
    double[] deltaCalculated;
}
