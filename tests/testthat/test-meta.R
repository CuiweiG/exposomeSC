# tests/testthat/test-meta.R

test_that("run_meta_exwas combines results correctly", {
    r1 <- S4Vectors::DataFrame(
        gene = c("A", "B", "C"),
        celltype = "Mono",
        log2FoldChange = c(0.5, -0.3, 0.1),
        lfcSE = c(0.1, 0.15, 0.2),
        exposure = "PM2.5")
    r2 <- S4Vectors::DataFrame(
        gene = c("A", "B"),
        celltype = "Mono",
        log2FoldChange = c(0.4, -0.2),
        lfcSE = c(0.12, 0.18),
        exposure = "PM2.5")

    meta <- run_meta_exwas(list(r1, r2),
        cohort_names = c("HELIX", "ENVR"))

    expect_s4_class(meta, "DataFrame")
    ## Only genes in >= 2 cohorts
    expect_true(nrow(meta) >= 2)
    expect_true("meta_effect" %in% colnames(meta))
    expect_true("I2" %in% colnames(meta))
    expect_true("direction" %in% colnames(meta))
    ## Gene C only in 1 cohort -> excluded
    expect_false("C" %in% meta$gene)
})

test_that("meta-analysis with random effects works", {
    set.seed(42)
    r1 <- S4Vectors::DataFrame(
        gene = paste0("G", 1:5), celltype = "Mono",
        log2FoldChange = rnorm(5, 0.5, 0.1),
        lfcSE = rep(0.1, 5), exposure = "Pb")
    r2 <- S4Vectors::DataFrame(
        gene = paste0("G", 1:5), celltype = "Mono",
        log2FoldChange = rnorm(5, 0.3, 0.2),
        lfcSE = rep(0.15, 5), exposure = "Pb")
    r3 <- S4Vectors::DataFrame(
        gene = paste0("G", 1:5), celltype = "Mono",
        log2FoldChange = rnorm(5, 0.4, 0.15),
        lfcSE = rep(0.12, 5), exposure = "Pb")

    if (!requireNamespace("metafor", quietly = TRUE)) {
        expect_error(
            run_meta_exwas(list(r1, r2, r3), method = "random"),
            "metafor"
        )
    } else {
        meta <- run_meta_exwas(list(r1, r2, r3), method = "random")
        expect_equal(nrow(meta), 5L)
        expect_true(all(meta$n_cohorts == 3L))
        expect_true(all(meta$I2 >= 0 & meta$I2 <= 1))
        expect_true(all(meta$method == "REML_Hartung-Knapp"))
        expect_true(all(meta$meta_df == 2L))
    }
})

test_that("direction string encodes signs correctly", {
    r1 <- S4Vectors::DataFrame(
        gene = "A", celltype = "T",
        log2FoldChange = 0.5, lfcSE = 0.1,
        exposure = "X")
    r2 <- S4Vectors::DataFrame(
        gene = "A", celltype = "T",
        log2FoldChange = -0.3, lfcSE = 0.1,
        exposure = "X")

    meta <- run_meta_exwas(list(r1, r2))
    ## Should show discordant direction in ASCII
    expect_true(nchar(meta$direction[1]) == 2)
    expect_setequal(strsplit(meta$direction[1], "")[[1]], c("+", "-"))
})

test_that("fixed-effect estimate and uncertainty are exact", {
    r1 <- S4Vectors::DataFrame(
        gene = "A", celltype = "T", log2FC = 0.5, se = 0.1,
        exposure = "X"
    )
    r2 <- S4Vectors::DataFrame(
        gene = "A", celltype = "T", log2FC = 0.2, se = 0.2,
        exposure = "X"
    )
    result <- run_meta_exwas(
        list(r1, r2),
        cohort_names = c("C1", "C2")
    )
    weights <- c(100, 25)
    expected_effect <- sum(weights * c(0.5, 0.2)) / sum(weights)
    expect_equal(result$meta_effect, expected_effect, tolerance = 1e-12)
    expect_equal(result$meta_se, sqrt(1 / sum(weights)), tolerance = 1e-12)
    expect_equal(result$cohorts, "C1;C2")
    expect_equal(result$method, "inverse_variance_fixed")
})

test_that("current and legacy backend schemas can be synthesised together", {
    current <- S4Vectors::DataFrame(
        gene = "A", celltype = "T", exposure = "X",
        log2FC = 0.5, se = 0.1
    )
    legacy <- S4Vectors::DataFrame(
        gene = "A", celltype = "T", exposure = "X",
        log2FoldChange = 0.2, lfcSE = 0.2
    )
    result <- run_meta_exwas(
        list(current, legacy),
        cohort_names = c("current", "legacy")
    )

    expect_equal(result$meta_effect, 0.44, tolerance = 1e-12)
    expect_equal(result$n_cohorts, 2L)
    expect_equal(result$cohorts, "current;legacy")
})

test_that("exposure is part of the exact synthesis key", {
    cohort <- function(offset) {
        S4Vectors::DataFrame(
            gene = c("A", "A"),
            celltype = c("T", "T"),
            exposure = c("PM2.5", "Pb"),
            log2FC = c(0.4, -0.2) + offset,
            se = c(0.1, 0.1)
        )
    }
    result <- run_meta_exwas(list(cohort(0), cohort(0.1)))
    expect_equal(nrow(result), 2L)
    expect_setequal(result$exposure, c("PM2.5", "Pb"))
    expect_true(result$meta_effect[result$exposure == "PM2.5"] > 0)
    expect_true(result$meta_effect[result$exposure == "Pb"] < 0)
})

test_that("duplicate keys and malformed cohort labels are rejected", {
    duplicate <- S4Vectors::DataFrame(
        gene = c("A", "A"),
        celltype = c("T", "T"),
        exposure = c("X", "X"),
        log2FC = c(0.2, 0.3),
        se = c(0.1, 0.1)
    )
    valid <- duplicate[1, ]
    expect_error(
        run_meta_exwas(list(duplicate, valid)),
        "Duplicate gene-celltype-exposure key"
    )
    expect_error(
        run_meta_exwas(
            list(valid, valid),
            cohort_names = c("same", "same")
        ),
        "unique"
    )
})

test_that("edgeR QL results cannot enter inverse-variance synthesis", {
    edgeR_proxy <- S4Vectors::DataFrame(
        gene = "A",
        celltype = "T",
        exposure = "X",
        log2FC = 0.5,
        se = 0.1,
        statistic = 5,
        pvalue = 0.01,
        method = "edgeR_robust_QL"
    )
    edgeR_current <- edgeR_proxy
    edgeR_current$se <- NA_real_

    expect_error(
        run_meta_exwas(list(edgeR_proxy, edgeR_proxy)),
        "edgeR quasi-likelihood.*does not supply coefficient standard errors"
    )
    expect_error(
        run_meta_exwas(list(edgeR_current, edgeR_current)),
        "edgeR quasi-likelihood.*does not supply coefficient standard errors"
    )
})
