package io.github.yok.bdg.app;

import java.util.List;
import javax.validation.Valid;
import javax.validation.constraints.NotEmpty;
import lombok.Data;
import lombok.Getter;
import lombok.Setter;
import lombok.ToString;
import org.springframework.boot.context.properties.ConfigurationProperties;
import org.springframework.validation.annotation.Validated;

/**
 * bdg-solver の設定値（bdg.*）を保持するクラスです。
 *
 * <p>
 * application.yml などから読み込まれ、CLI 実行時の初期化に使用します。
 * </p>
 */
@Data
@ToString(onlyExplicitlyIncluded = true)
@Validated
@ConfigurationProperties(prefix = "bdg")
public class BdgProperties {

    /**
     * 格子設定です。
     */
    private Lattice lattice = new Lattice();

    /**
     * モデル設定です。
     */
    private Model model = new Model();

    /**
     * トラップ設定です。
     */
    private Trap trap = new Trap();

    /**
     * 自己無撞着（SCF）設定です。
     */
    private Scf scf = new Scf();

    /**
     * 化学ポテンシャル探索（μ/h 二分法）の設定です。
     */
    private ChemicalPotentialSearch chemicalPotentialSearch = new ChemicalPotentialSearch();

    /**
     * 出力設定です。
     */
    private Output output = new Output();

    /**
     * 設定値を YAML 風の複数行文字列に整形して返します。
     *
     * @return 設定値の整形文字列です
     */
    @ToString.Include(name = "bdg")
    public String toMultilineString() {
        String nl = System.lineSeparator();

        Lattice l = getLattice();
        Model m = getModel();
        Trap t = getTrap();
        Scf s = getScf();
        ChemicalPotentialSearch cps = getChemicalPotentialSearch();
        Output o = getOutput();

        // 先頭改行を入れて、ログの可読性を上げます。
        StringBuilder sb = new StringBuilder(256).append(nl);

        appendSection(sb, nl, "lattice",
                // lx: x方向の格子点数
                "lx", l.getLx(),
                // ly: y方向の格子点数
                "ly", l.getLy(),
                // boundary: 境界条件（OPEN/PERIODIC）
                "boundary", l.getBoundary());

        appendSection(sb, nl, "model",
                // t: ホッピング（最近接の移動振幅）
                "t", m.getT(),
                // u: 相互作用強度（オンサイト）
                "u", m.getU(),
                // temperature: 温度（エネルギー単位）
                "temperature", m.getTemperature(),
                // totalParticleNumber: 総粒子数 N
                "totalParticleNumber", m.getTotalParticleNumber(),
                // muInitial: 探索用の化学ポテンシャル初期値
                "muInitial", m.getMuInitial(),
                // muInitialEstimation.enabled: mu-initial 自動推定を有効化するかどうか
                "muInitialEstimation.enabled", m.getMuInitialEstimation().isEnabled(),
                // muInitialEstimation.tolerance: |N(mu)-Ntarget| <= tolerance で推定を打ち切る
                "muInitialEstimation.tolerance", m.getMuInitialEstimation().getTolerance(),
                // muInitialEstimation.maxIterations: 二分法の最大反復回数
                "muInitialEstimation.maxIterations", m.getMuInitialEstimation().getMaxIterations(),
                // muInitialEstimation.bracketMargin: 挟み込み区間の余裕（探索レンジ拡張用マージン）
                "muInitialEstimation.bracketMargin", m.getMuInitialEstimation().getBracketMargin(),
                // muInitialEstimation.hartreeEnabled: Δなし Hartree-LDA を推定に含めるかどうか
                "muInitialEstimation.hartreeEnabled", m.getMuInitialEstimation().isHartreeEnabled(),
                // muInitialEstimation.hartreeMaxIterations: Hartree 固定点反復の最大反復回数
                "muInitialEstimation.hartreeMaxIterations",
                m.getMuInitialEstimation().getHartreeMaxIterations(),
                // muInitialEstimation.hartreeTolerance: Hartree 固定点反復の収束判定閾値
                "muInitialEstimation.hartreeTolerance",
                m.getMuInitialEstimation().getHartreeTolerance(),
                // muInitialEstimation.hartreeMixing: Hartree 固定点反復の緩和係数
                "muInitialEstimation.hartreeMixing", m.getMuInitialEstimation().getHartreeMixing(),
                // zeemanFieldScan.values: 計算する Zeeman 場 h の一覧
                "zeemanFieldScan.values", m.getZeemanFieldScan().getValues());

        appendSection(sb, nl, "trap",
                // r0: 基準半径
                "r0", t.getR0(),
                // omegaHo: 調和トラップ係数
                "omegaHo", t.getOmegaHo(),
                // omegaTr: トロイダル項係数
                "omegaTr", t.getOmegaTr(),
                // xi: 減衰長
                "xi", t.getXi());

        appendSection(sb, nl, "scf",
                // maxIter: 最大反復回数
                "maxIter", s.getMaxIter(),
                // mixing.densityAlpha: 密度（n1, n2）に適用する mixing 係数
                "mixing.densityAlpha", s.getMixing().getDensityAlpha(),
                // mixing.orderParameterAlpha: 秩序パラメータ（Δ）に適用する mixing 係数
                "mixing.orderParameterAlpha", s.getMixing().getOrderParameterAlpha(),
                // tolerance.densityAbs: 密度（n1, n2）の残差（絶対値）の許容上限
                "tolerance.densityAbs", s.getTolerance().getDensityAbs(),
                // tolerance.orderParameterAbs: 秩序パラメータ（Δ）の残差（絶対値）の許容上限
                "tolerance.orderParameterAbs", s.getTolerance().getOrderParameterAbs(),
                // tolerance.orderParameterRmsAbs: 秩序パラメータ（Δ）の残差RMS（絶対値）の許容上限
                "tolerance.orderParameterRmsAbs", s.getTolerance().getOrderParameterRmsAbs());

        appendSection(sb, nl, "chemicalPotentialSearch",
                // muMaxIterations: μ 探索の最大反復回数
                "muMaxIterations", cps.getMuMaxIterations(),
                // particleTolerance: 粒子数の許容誤差
                "particleTolerance", cps.getParticleTolerance(),
                // warmStart: 直前の収束状態を初期値として再利用するかどうか
                "warmStart", cps.isWarmStart());

        appendSection(sb, nl, "output",
                // dir: 出力先ディレクトリ
                "dir", o.getDir());

        return sb.toString();
    }

    /**
     * セクション名と (key, value) ペア列を、YAML 風の複数行テキストとして追記します。
     *
     * <pre>
     *   section:
     *     key: value
     * </pre>
     *
     * @param sb 追記先バッファです
     * @param nl 改行文字列です
     * @param section セクション名です
     * @param kvPairs key1, value1, key2, value2, ... の順で渡すペア列です
     */
    private static void appendSection(StringBuilder sb, String nl, String section,
            Object... kvPairs) {
        sb.append("  ").append(section).append(":").append(nl);
        for (int i = 0; i < kvPairs.length; i += 2) {
            String key = String.valueOf(kvPairs[i]);
            Object val = (i + 1 < kvPairs.length) ? kvPairs[i + 1] : null;
            sb.append("    ").append(key).append(": ").append(val).append(nl);
        }
    }

    @Data
    public static class Lattice {

        /**
         * x方向サイズ（Lx）です。
         */
        private int lx = 12;

        /**
         * y方向サイズ（Ly）です。
         */
        private int ly = 12;

        /**
         * 境界条件です。
         */
        private Boundary boundary = Boundary.OPEN;

        public enum Boundary {
            OPEN, PERIODIC
        }
    }

    @Data
    public static class Model {

        /**
         * ホッピング（最近接の移動振幅）です。
         */
        private double t = 1.0;

        /**
         * 相互作用強度（オンサイト）です。
         */
        private double u = -4.0;

        /**
         * 温度（エネルギー単位）です。
         */
        private double temperature = 0.001;

        /**
         * 総粒子数 N (=N1+N2) です。
         */
        private int totalParticleNumber = 40;

        /**
         * 探索用の化学ポテンシャル初期値です。
         */
        private double muInitial = 0.0;

        /**
         * 化学ポテンシャルの自動推定設定です。
         */
        private MuInitialEstimation muInitialEstimation = new MuInitialEstimation();

        @Getter
        @Setter
        public static class MuInitialEstimation {

            /**
             * mu-initial の自動推定を有効化するかどうかです。
             *
             * <p>
             * true の場合、設定値 {@code bdg.model.mu-initial} は推定結果で上書きされます。
             * </p>
             */
            private boolean enabled = false;

            /**
             * 推定の打ち切り許容誤差です。
             *
             * <p>
             * {@code |N(mu)-Ntarget| <= tolerance} を満たした時点で推定を終了します。
             * </p>
             */
            private double tolerance = 0.5;

            /**
             * 二分法の最大反復回数です。
             */
            private int maxIterations = 80;

            /**
             * 挟み込み区間の余裕（探索レンジ拡張用マージン）です。
             */
            private double bracketMargin = 8.0;

            /**
             * Hartree（Un）を LDA 推定に含めるかどうかです。
             *
             * <p>
             * true の場合、Δなしの Hartree-LDA により {@code N(mu)} を構成します。 BdG-SCF の N(μ)
             * に寄るため、初期μのズレが減る効果が期待できます。
             * </p>
             */
            private boolean hartreeEnabled = true;

            /**
             * Hartree の固定点反復の最大反復回数です。
             */
            private int hartreeMaxIterations = 30;

            /**
             * Hartree の固定点反復の収束判定閾値です。
             */
            private double hartreeTolerance = 1e-4;

            /**
             * Hartree 固定点反復の緩和係数です（0より大きく1以下）。
             */
            private double hartreeMixing = 0.5;
        }

        /**
         * Zeeman 場（有効磁場）h を指定して計算を繰り返す設定です。
         *
         * <p>
         * ここでの h は {@code h=(mu1-mu2)/2} を表します。
         * </p>
         */
        @Valid
        private ZeemanFieldScan zeemanFieldScan = new ZeemanFieldScan();

        @Data
        public static class ZeemanFieldScan {

            /**
             * 計算する Zeeman 場 h の一覧です。
             */
            @NotEmpty
            private List<Double> values = List.of();
        }
    }

    @Data
    public static class Trap {

        /**
         * 基準半径です。
         */
        private double r0 = 6.0;

        /**
         * 調和トラップ係数です。
         */
        private double omegaHo = 12.0;

        /**
         * トロイダル項係数です。
         */
        private double omegaTr = 8.0;

        /**
         * 減衰長です。
         */
        private double xi = 5.0;
    }

    @Data
    public static class Scf {

        /**
         * 最大反復回数です。
         */
        private int maxIter = 50;

        /**
         * mixing 設定です。
         */
        @Valid
        private Mixing mixing = new Mixing();

        /**
         * 収束判定（残差）の設定です。
         */
        @Valid
        private Tolerance tolerance = new Tolerance();

        @Data
        public static class Mixing {

            /**
             * 密度（n1, n2）に適用する mixing 係数です。
             */
            private double densityAlpha = 0.1;

            /**
             * 秩序パラメータ（Δ）に適用する mixing 係数です。
             */
            private double orderParameterAlpha = 0.1;
        }

        @Data
        public static class Tolerance {

            /**
             * 密度（n1, n2）の残差（絶対値）の許容上限です。
             */
            private double densityAbs = 1e-6;

            /**
             * 秩序パラメータ（Δ）の残差（絶対値）の許容上限です。
             */
            private double orderParameterAbs = 1e-6;

            /**
             * 秩序パラメータ（Δ）の収束判定に用いる許容誤差（RMS 残差の絶対値）です。
             * <p>
             * Δ が空間依存（FFLO など）して符号反転しても、差分のRMSは相殺しないため安定です。
             * </p>
             */
            private double orderParameterRmsAbs;
        }
    }

    /**
     * 化学ポテンシャル探索（μ/h 二分法）の設定です。
     */
    @Data
    public static class ChemicalPotentialSearch {

        /**
         * μ 探索の最大反復回数です。
         */
        private int muMaxIterations = 30;

        /**
         * 粒子数の許容誤差です。
         */
        private double particleTolerance = 1e-3;

        /**
         * 直前の収束状態を初期値として再利用するかどうかです。
         */
        private boolean warmStart = true;
    }

    @Data
    public static class Output {

        /**
         * 出力先ディレクトリです。
         */
        private String dir = "./out";
    }
}
