package io.github.yok.bdg.app;

import io.github.yok.bdg.core.solver.ChemicalPotentialOptimizer;
import io.github.yok.bdg.core.state.MeanFieldState;
import io.github.yok.bdg.core.state.MeanFieldStateInitializer;
import io.github.yok.bdg.out.ResultWriter;
import java.util.List;
import lombok.RequiredArgsConstructor;
import org.springframework.boot.CommandLineRunner;
import org.springframework.stereotype.Component;

/**
 * CLI で bdg-solver を実行するクラスです。
 *
 * <p>
 * Zeeman 場 h をスキャンし、各 h について総粒子数 N が目標値になるように 平均化学ポテンシャル μ（結果としての μ1, μ2）を探索します。
 * </p>
 */
@Component
@RequiredArgsConstructor
public class BdgCliRunner implements CommandLineRunner {

    /**
     * bdg-solver の設定値（bdg.*）です。
     */
    private final BdgProperties properties;

    /**
     * 平均場状態の初期値生成ロジックです。
     */
    private final MeanFieldStateInitializer meanFieldStateInitializer;

    /**
     * 化学ポテンシャル探索ロジックです。
     */
    private final ChemicalPotentialOptimizer chemicalPotentialOptimizer;

    /**
     * 結果出力ロジックです。
     */
    private final ResultWriter resultWriter;

    /**
     * CLI 実行を開始します。
     *
     * @param args 起動引数です
     */
    @Override
    public void run(String... args) {
        System.out.println("=== bdg-solver start: optimize chemical potentials ===");
        System.out.print(properties.toMultilineString());

        int totalN = properties.getModel().getTotalParticleNumber();
        if (totalN <= 0) {
            throw new IllegalStateException("model.totalParticleNumber は 1 以上を指定してください: " + totalN);
        }

        // h の一覧（Zeeman 場）
        List<Double> zeemanFields = properties.getModel().getZeemanFieldScan().getValues();
        if (zeemanFields == null || zeemanFields.isEmpty()) {
            throw new IllegalStateException("zeemanFieldScan.values は必須です（h の一覧を指定してください）");
        }

        boolean warmStart = properties.getChemicalPotentialSearch().isWarmStart();

        MeanFieldState previousState = null;

        // Zeeman 場ごとに化学ポテンシャル探索を実行
        for (int i = 0; i < zeemanFields.size(); i++) {
            Double hObj = zeemanFields.get(i);
            if (hObj == null) {
                throw new IllegalStateException("zeemanFieldScan.values に null が含まれています");
            }
            double h = hObj.doubleValue();

            // warmStart=true の場合は直前の収束状態を再利用、そうでなければ初期値を新規生成
            MeanFieldState seed = (warmStart && previousState != null) ? previousState
                    : meanFieldStateInitializer.create();

            System.out.println("=== Zeeman場ごとの計算 ===");
            System.out.println("入力: h=" + fmt5(h) + ", 目標N=" + totalN + "（warmStart=" + warmStart
                    + ", step=" + (i + 1) + "/" + zeemanFields.size() + "）");

            ChemicalPotentialOptimizer.OptimizationResult opt =
                    chemicalPotentialOptimizer.optimize(seed, totalN, h);

            // optimizer が返した「最終状態」を採用（実装差し替え等で state が返らない場合の保険あり）
            MeanFieldState solved = opt.getState();
            if (solved == null) {
                seed.setMu1(opt.getMu1());
                seed.setMu2(opt.getMu2());
                solved = seed;
            }

            // 結果（粒子数・偏極）
            double actualN1 = solved.totalN1();
            double actualN2 = solved.totalN2();
            double actualN = actualN1 + actualN2;
            double pOut = (actualN != 0.0) ? (actualN1 - actualN2) / actualN : 0.0;

            // 出力（N, h, P_out を渡す）
            resultWriter.write(solved, opt, totalN, h, pOut);

            System.out.println("結果: N=" + fmt5(actualN) + " (N1=" + fmt5(actualN1) + ", N2="
                    + fmt5(actualN2) + ")" + ", P_out=" + fmt5(pOut));
            System.out.println("結果: mu1=" + fmt5(opt.getMu1()) + ", mu2=" + fmt5(opt.getMu2())
                    + ", iterations=" + opt.getIterations());

            if (warmStart) {
                previousState = solved;
            }
        }
    }

    /**
     * 数値を小数点以下5桁までの文字列に整形します。
     *
     * @param v 数値です
     * @return 整形した文字列です
     */
    private static String fmt5(double v) {
        return String.format(java.util.Locale.ROOT, "%.5f", v);
    }
}
