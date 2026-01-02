package io.github.yok.bdg;

import org.springframework.boot.SpringApplication;
import org.springframework.boot.autoconfigure.SpringBootApplication;
import org.springframework.boot.context.properties.ConfigurationPropertiesScan;

/**
 * bdg-solver のエントリポイントです。
 *
 * <p>
 * 設定クラス（@ConfigurationProperties）をスキャンし、CLI 実行を開始します。
 * </p>
 */
@SpringBootApplication
@ConfigurationPropertiesScan(basePackages = "io.github.yok.bdg")
public class BdgSolverApplication {

    /**
     * Spring Boot アプリケーションを起動します。
     *
     * @param args 起動引数です
     */
    public static void main(String[] args) {
        SpringApplication.run(BdgSolverApplication.class, args);
    }
}
