package io.github.yok.bdg.app;

import io.github.yok.bdg.core.lattice.SquareLattice2D;
import io.github.yok.bdg.core.linearalgebra.EigenDecompositionBackend;
import io.github.yok.bdg.core.linearalgebra.EjmlSymmetricEigenDecompositionBackend;
import io.github.yok.bdg.core.model.BdgModel;
import io.github.yok.bdg.core.model.MeanFieldUpdater;
import io.github.yok.bdg.core.model.OnsiteSWaveModel;
import io.github.yok.bdg.core.model.OnsiteSWaveUpdater;
import io.github.yok.bdg.core.model.ToroidalTrapPotential;
import io.github.yok.bdg.core.solver.ChemicalPotentialBisectionOptimizer;
import io.github.yok.bdg.core.solver.ChemicalPotentialOptimizer;
import io.github.yok.bdg.core.solver.SelfConsistencySolver;
import io.github.yok.bdg.core.state.MeanFieldStateInitializer;
import io.github.yok.bdg.core.state.ToroidalLdaRingInitializer;
import io.github.yok.bdg.out.CsvResultWriter;
import io.github.yok.bdg.out.ResultWriter;
import lombok.RequiredArgsConstructor;
import org.springframework.context.annotation.Bean;
import org.springframework.context.annotation.Configuration;

/**
 * トロイダルトラップ + オンサイト s-wave + SCF 解法の Bean 定義を行う設定クラスです。
 *
 * <p>
 * ToroidalTrapPotential と OnsiteSWaveModel を用いて、平均場 SCF ソルバ一式を組み立てます。
 * </p>
 */
@Configuration
@RequiredArgsConstructor
public class ToroidalOnsiteSWaveScfConfiguration {

    /**
     * bdg-solver の設定値（bdg.*）です。
     */
    private final BdgProperties p;

    /**
     * 格子を生成します。
     *
     * @return 格子です
     */
    @Bean
    public SquareLattice2D lattice() {
        return new SquareLattice2D(p.getLattice().getLx(), p.getLattice().getLy(),
                p.getLattice().getBoundary());
    }

    /**
     * トロイダルトラップポテンシャルを生成します。
     *
     * @param lattice 格子です
     * @return トラップポテンシャルです
     */
    @Bean
    public ToroidalTrapPotential toroidalTrapPotential(SquareLattice2D lattice) {
        BdgProperties.Trap t = p.getTrap();
        return new ToroidalTrapPotential(lattice, t.getR0(), t.getOmegaHo(), t.getOmegaTr(),
                t.getXi());
    }

    /**
     * オンサイト s-wave BdG モデルを生成します。
     *
     * @param lattice 格子です
     * @param trapPotential トラップポテンシャルです
     * @return BdG モデルです
     */
    @Bean
    public BdgModel onsiteSWaveModel(SquareLattice2D lattice, ToroidalTrapPotential trapPotential) {
        return new OnsiteSWaveModel(lattice, trapPotential, p.getModel().getT(),
                p.getModel().getU());
    }

    /**
     * 平均場状態の初期値生成ロジックを生成します。
     *
     * @param model BdG モデルです
     * @param trapPotential トラップポテンシャルです
     * @return 初期値生成ロジックです
     */
    @Bean
    public MeanFieldStateInitializer meanFieldStateInitializer(BdgModel model,
            ToroidalTrapPotential trapPotential) {
        return new ToroidalLdaRingInitializer(p, model, trapPotential);
    }

    /**
     * 固有分解バックエンドを生成します。
     *
     * @return 固有分解バックエンドです
     */
    @Bean
    public EigenDecompositionBackend eigenDecompositionBackend() {
        return new EjmlSymmetricEigenDecompositionBackend();
    }

    /**
     * 平均場更新式（密度と秩序変数の更新）を生成します。
     *
     * @return 更新式です
     */
    @Bean
    public MeanFieldUpdater meanFieldUpdater() {
        return new OnsiteSWaveUpdater(p.getModel().getU(), p.getModel().getTemperature());
    }

    /**
     * 平均場 SCF ソルバを生成します。
     *
     * @param model BdG モデルです
     * @param eigen 固有分解バックエンドです
     * @param updater 更新式です
     * @return SCF ソルバです
     */
    @Bean
    public SelfConsistencySolver selfConsistencySolver(BdgModel model,
            EigenDecompositionBackend eigen, MeanFieldUpdater updater) {
        BdgProperties.Scf scf = p.getScf();
        return new SelfConsistencySolver(model, eigen, updater, scf);
    }

    /**
     * 粒子数条件に合わせて μ1, μ2 を探索する最適化ロジックを生成します。
     *
     * @param scfSolver 平均場 SCF ソルバです
     * @return 化学ポテンシャル探索ロジックです
     */
    @Bean
    public ChemicalPotentialOptimizer chemicalPotentialOptimizer(SelfConsistencySolver scfSolver) {
        BdgProperties.ChemicalPotentialSearch s = p.getChemicalPotentialSearch();
        return new ChemicalPotentialBisectionOptimizer(scfSolver, s.getMuMaxIterations(),
                s.getParticleTolerance());
    }

    /**
     * 結果出力ロジックを生成します。
     *
     * @return 結果出力ロジックです
     */
    @Bean
    public ResultWriter resultWriter() {
        return new CsvResultWriter(p.getOutput().getDir(), p.getLattice().getLx(),
                p.getLattice().getLy());
    }
}
