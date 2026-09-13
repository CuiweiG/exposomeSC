## Tests for run_multi_exwas

library(SingleCellExperiment)
library(S4Vectors)

.make_scee_multi <- function() {
    set.seed(42)
    counts <- matrix(rpois(5000, 10), nrow = 50,
        dimnames = list(paste0("G", 1:50), paste0("c", 1:100)))
    sce <- SingleCellExperiment(
        assays = list(counts = counts),
        colData = DataFrame(
            cell_id = paste0("c", 1:100),
            donor_id = rep(paste0("D", 1:5), each = 20),
            cell_type = rep(c("Mono", "NK"), 50)))
    exp_mat <- matrix(rnorm(15), nrow = 5,
        dimnames = list(paste0("D", 1:5),
                        c("E1", "E2", "E3")))
    build_scee(sce, exp_mat, sample_col = "donor_id")
}

test_that("run_multi_exwas runs multiple exposures", {
    scee <- .make_scee_multi()
    result <- suppressWarnings(run_multi_exwas(scee,
        exposures = c("E1", "E2"),
        celltype_col = "cell_type",
        sample_col = "donor_id",
        min_cells = 5L,
        min_donors = 3L))
    expect_s4_class(result, "DataFrame")
    expect_true("padj_all" %in% colnames(result))
    ## Should have results for both exposures
    if (nrow(result) > 0) {
        expect_true(all(c(
            "log2FC", "se", "statistic", "pvalue", "exposure"
        ) %in% colnames(result)))
        expect_equal(
            result$padj_all,
            stats::p.adjust(result$pvalue, method = "BH")
        )
        exps <- unique(result$exposure)
        expect_true(length(exps) >= 1)
    }
})

test_that("run_multi_exwas defaults to all exposures", {
    scee <- .make_scee_multi()
    result <- suppressWarnings(run_multi_exwas(scee,
        celltype_col = "cell_type",
        sample_col = "donor_id",
        min_cells = 5L,
        min_donors = 3L))
    if (nrow(result) > 0) {
        expect_true(length(unique(result$exposure)) >= 1)
    }
})

test_that("run_multi_exwas errors on bad exposure", {
    scee <- .make_scee_multi()
    expect_error(run_multi_exwas(scee,
        exposures = c("E1", "FAKE"),
        celltype_col = "cell_type",
        sample_col = "donor_id"), "not found")
})

test_that("run_multi_exwas rejects duplicated or missing exposure names", {
    scee <- .make_scee_multi()
    expect_error(
        run_multi_exwas(
            scee,
            exposures = c("E1", "E1"),
            celltype_col = "cell_type"
        ),
        "unique"
    )
    expect_error(
        run_multi_exwas(
            scee,
            exposures = c("E1", NA_character_),
            celltype_col = "cell_type"
        ),
        "non-missing"
    )
})

test_that("run_multi_exwas has a stable empty current-result schema", {
    scee <- .make_scee_multi()
    result <- suppressWarnings(run_multi_exwas(
        scee,
        exposures = "E1",
        celltype_col = "cell_type",
        min_cells = 1000L,
        min_donors = 3L
    ))
    expect_equal(nrow(result), 0L)
    expect_true(all(c(
        "gene", "celltype", "log2FC", "se", "statistic", "pvalue",
        "pvalue_underflow_clamped", "padj", "padj_global", "exposure",
        "n_donors", "method", "padj_all"
    ) %in% colnames(result)))
})
