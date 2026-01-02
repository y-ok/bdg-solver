package io.github.yok.bdg.out;

import io.github.yok.bdg.core.solver.ChemicalPotentialOptimizer;
import io.github.yok.bdg.core.state.MeanFieldState;
import java.io.IOException;
import java.io.Writer;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Locale;
import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVPrinter;

/**
 * 計算結果を CSV に出力するクラスです。
 *
 * <p>
 * 出力ファイル名は、以下の命名規約に従います（N は総粒子数、h は Zeeman 場、P_out は結果偏極率）。
 * </p>
 *
 * <ul>
 * <li>{@code bdg_orderParameter_N=800_h=0.20.csv}</li>
 * <li>{@code bdg_magnetization_N=800_h=0.20.csv}</li>
 * <li>{@code bdg_density_N=800_h=0.20.csv}</li>
 * <li>{@code bdg_meta_N=800_h=0.20.csv}（μ1/μ2 等の補助情報、P_out など）</li>
 * </ul>
 *
 * <p>
 * 2次元空間（x,y）で可視化する前提のため、各CSVは x,y,i を必ず含みます。
 * </p>
 */
public final class CsvResultWriter implements ResultWriter {

    /**
     * ファイル名の先頭固定文字列です。
     *
     * <p>
     * 命名規約により {@code "bdg"} 固定とします。
     * </p>
     */
    private static final String FILE_HEAD = "bdg";

    /**
     * 出力先ディレクトリです。
     */
    private final Path outputDir;

    /**
     * x方向格子点数（Lx）です。
     */
    private final int lx;

    /**
     * y方向格子点数（Ly）です。
     */
    private final int ly;

    /**
     * CSV 出力を生成します。
     *
     * @param outputDir 出力先ディレクトリです
     * @param lx x方向格子点数（1以上）
     * @param ly y方向格子点数（1以上）
     * @throws IllegalArgumentException 引数が不正な場合に発生します
     */
    public CsvResultWriter(String outputDir, int lx, int ly) {
        if (outputDir == null || outputDir.isEmpty()) {
            throw new IllegalArgumentException("output.dir は必須です");
        }
        if (lx <= 0 || ly <= 0) {
            throw new IllegalArgumentException("lx/ly は 1 以上を指定してください: " + lx + "x" + ly);
        }
        this.outputDir = Paths.get(outputDir);
        this.lx = lx;
        this.ly = ly;
    }

    /**
     * 最適化結果と平均場状態を出力します。
     *
     * @param state 平均場状態です
     * @param optimization 化学ポテンシャル探索結果です
     * @param totalParticleNumber 目標の総粒子数 N（N1+N2）です
     * @param zeemanField Zeeman 場 h（有効磁場、h=(mu1-mu2)/2）です
     * @param outputPolarization 結果として得られた偏極率 P_out（(N1-N2)/(N1+N2)）です
     * @throws IllegalArgumentException 引数が不正な場合に発生します
     * @throws IllegalStateException 出力に失敗した場合に発生します
     */
    @Override
    public void write(MeanFieldState state,
            ChemicalPotentialOptimizer.OptimizationResult optimization, int totalParticleNumber,
            double zeemanField, double outputPolarization) {

        if (state == null) {
            throw new IllegalArgumentException("state は null 不可です");
        }
        if (optimization == null) {
            throw new IllegalArgumentException("optimization は null 不可です");
        }
        if (totalParticleNumber <= 0) {
            throw new IllegalArgumentException("総粒子数 N は 1 以上を指定してください: " + totalParticleNumber);
        }
        if (!Double.isFinite(zeemanField)) {
            throw new IllegalArgumentException("Zeeman 場 h は有限値を指定してください: " + zeemanField);
        }
        if (!Double.isFinite(outputPolarization)) {
            throw new IllegalArgumentException("P_out は有限値を指定してください: " + outputPolarization);
        }
        if (state.getSiteCount() != lx * ly) {
            throw new IllegalArgumentException("siteCount と lx*ly が一致しません: siteCount="
                    + state.getSiteCount() + ", lx*ly=" + (lx * ly));
        }

        try {
            Files.createDirectories(outputDir);

            // 1) Δ（秩序パラメータ）
            writeOrderParameterCsv(state, totalParticleNumber, zeemanField);

            // 2) 磁化（m=n1-n2）
            writeMagnetizationCsv(state, totalParticleNumber, zeemanField);

            // 3) 密度（n=n1+n2）
            writeDensityCsv(state, totalParticleNumber, zeemanField);

            // 4) メタ（μ1/μ2、反復、P_out、N1/N2合計など）
            writeMetaCsv(state, optimization, totalParticleNumber, zeemanField, outputPolarization);

        } catch (IOException e) {
            throw new IllegalStateException("CSV 出力に失敗しました: " + outputDir, e);
        }
    }

    /**
     * 秩序パラメータ（Δ）を出力します。
     *
     * @param state 平均場状態です
     * @param totalParticleNumber 総粒子数 N です
     * @param zeemanField Zeeman 場 h です
     * @throws IOException 出力に失敗した場合に発生します
     */
    private void writeOrderParameterCsv(MeanFieldState state, int totalParticleNumber,
            double zeemanField) throws IOException {

        Path file = outputDir
                .resolve(buildFileName("orderParameter", totalParticleNumber, zeemanField));

        double[] delta = state.getDelta();

        try (Writer w = Files.newBufferedWriter(file, StandardCharsets.UTF_8);
                CSVPrinter pr = CSVFormat.Builder.create(CSVFormat.DEFAULT)
                        .setHeader("x", "y", "i", "orderParameter").build().print(w)) {

            for (int y = 0; y < ly; y++) {
                for (int x = 0; x < lx; x++) {
                    int i = y * lx + x;
                    pr.printRecord(x, y, i, delta[i]);
                }
            }
        }
    }

    /**
     * 磁化（m=n1-n2）を出力します。
     *
     * @param state 平均場状態です
     * @param totalParticleNumber 総粒子数 N です
     * @param zeemanField Zeeman 場 h です
     * @throws IOException 出力に失敗した場合に発生します
     */
    private void writeMagnetizationCsv(MeanFieldState state, int totalParticleNumber,
            double zeemanField) throws IOException {

        Path file =
                outputDir.resolve(buildFileName("magnetization", totalParticleNumber, zeemanField));

        double[] n1 = state.getN1();
        double[] n2 = state.getN2();

        try (Writer w = Files.newBufferedWriter(file, StandardCharsets.UTF_8);
                CSVPrinter pr = CSVFormat.Builder.create(CSVFormat.DEFAULT)
                        .setHeader("x", "y", "i", "magnetization").build().print(w)) {

            for (int y = 0; y < ly; y++) {
                for (int x = 0; x < lx; x++) {
                    int i = y * lx + x;
                    double m = n1[i] - n2[i];
                    pr.printRecord(x, y, i, m);
                }
            }
        }
    }

    /**
     * 密度（n=n1+n2）を出力します。
     *
     * @param state 平均場状態です
     * @param totalParticleNumber 総粒子数 N です
     * @param zeemanField Zeeman 場 h です
     * @throws IOException 出力に失敗した場合に発生します
     */
    private void writeDensityCsv(MeanFieldState state, int totalParticleNumber, double zeemanField)
            throws IOException {

        Path file = outputDir.resolve(buildFileName("density", totalParticleNumber, zeemanField));

        double[] n1 = state.getN1();
        double[] n2 = state.getN2();

        try (Writer w = Files.newBufferedWriter(file, StandardCharsets.UTF_8);
                CSVPrinter pr = CSVFormat.Builder.create(CSVFormat.DEFAULT)
                        .setHeader("x", "y", "i", "density").build().print(w)) {

            for (int y = 0; y < ly; y++) {
                for (int x = 0; x < lx; x++) {
                    int i = y * lx + x;
                    double n = n1[i] + n2[i];
                    pr.printRecord(x, y, i, n);
                }
            }
        }
    }

    /**
     * メタ情報を出力します。
     *
     * <p>
     * N1/N2 の総和および P_out は、出力された n1/n2 から計算した値も併記します。
     * </p>
     *
     * @param state 平均場状態です
     * @param optimization 最適化結果です
     * @param totalParticleNumber 総粒子数 N（入力値）です
     * @param zeemanField Zeeman 場 h（入力値）です
     * @param outputPolarization 結果として得られた偏極率 P_out（入力として渡された値）です
     * @throws IOException 出力に失敗した場合に発生します
     */
    private void writeMetaCsv(MeanFieldState state,
            ChemicalPotentialOptimizer.OptimizationResult optimization, int totalParticleNumber,
            double zeemanField, double outputPolarization) throws IOException {

        Path file = outputDir.resolve(buildFileName("meta", totalParticleNumber, zeemanField));

        double n1Total = sum(state.getN1());
        double n2Total = sum(state.getN2());
        double pComputed = computePolarization(n1Total, n2Total);

        try (Writer w = Files.newBufferedWriter(file, StandardCharsets.UTF_8);
                CSVPrinter pr = CSVFormat.Builder.create(CSVFormat.DEFAULT)
                        .setHeader("key", "value").build().print(w)) {

            pr.printRecord("input.N", totalParticleNumber);
            pr.printRecord("input.h", zeemanField);

            pr.printRecord("sum.N1", n1Total);
            pr.printRecord("sum.N2", n2Total);
            pr.printRecord("sum.N", n1Total + n2Total);

            pr.printRecord("output.P_out", outputPolarization);
            pr.printRecord("computed.P_out", pComputed);

            pr.printRecord("mu1", optimization.getMu1());
            pr.printRecord("mu2", optimization.getMu2());
            pr.printRecord("iterations", optimization.getIterations());

            pr.printRecord("lx", lx);
            pr.printRecord("ly", ly);
            pr.printRecord("siteCount", state.getSiteCount());
        }
    }

    /**
     * 命名規約に従ってファイル名を作成します。
     *
     * <p>
     * 例: {@code bdg_density_N=800_h=0.20.csv}
     * </p>
     *
     * @param kind 量の識別子（orderParameter/magnetization/density/meta）
     * @param totalParticleNumber 総粒子数 N です
     * @param zeemanField Zeeman 場 h です
     * @return ファイル名です
     */
    private static String buildFileName(String kind, int totalParticleNumber, double zeemanField) {
        return FILE_HEAD + "_" + kind + "_N=" + totalParticleNumber + "_h=" + formatH(zeemanField)
                + ".csv";
    }

    /**
     * Zeeman 場 h を小数点以下2桁に整形します（ファイル名用）。
     *
     * @param h Zeeman 場です
     * @return 整形文字列（例: 0.20）
     */
    private static String formatH(double h) {
        return String.format(Locale.ROOT, "%.2f", h);
    }

    /**
     * 配列の総和を計算します。
     *
     * @param v 配列です
     * @return 総和です
     */
    private static double sum(double[] v) {
        if (v == null) {
            return 0.0;
        }
        double s = 0.0;
        for (double x : v) {
            s += x;
        }
        return s;
    }

    /**
     * 偏極率 P_out を計算します。
     *
     * @param n1Total N1 の総和です
     * @param n2Total N2 の総和です
     * @return P_out=(N1-N2)/(N1+N2) です
     */
    private static double computePolarization(double n1Total, double n2Total) {
        double denom = n1Total + n2Total;
        if (denom == 0.0) {
            return 0.0;
        }
        return (n1Total - n2Total) / denom;
    }
}
