library(SingleCellExperiment)
library(S4Vectors)

.make_scee <- function(n_donors = 10, n_cells = 200) {
    set.seed(42)
    counts <- matrix(rpois(100 * n_cells, 8), nrow = 100,
        dimnames = list(paste0("G", 1:100),
                        paste0("c", seq_len(n_cells))))
    sce <- SingleCellExperiment(
        assays = list(counts = counts),
        colData = DataFrame(
            cell_id = colnames(counts),
            donor_id = rep(paste0("D", seq_len(n_donors)),
                           each = n_cells / n_donors),
            cell_type = rep(c("Mono", "NK"), n_cells / 2)))
    exp_mat <- matrix(rnorm(n_donors * 3), nrow = n_donors,
        dimnames = list(paste0("D", seq_len(n_donors)),
                        c("E1", "E2", "E3")))
    build_scee(sce, exp_mat, sample_col = "donor_id")
}

test_that("run_sc_mixture returns weights that sum to 1", {
    scee <- .make_scee()
    mix <- run_sc_mixture(scee,
        exposures = c("E1", "E2", "E3"),
        celltype = "Mono",
        celltype_col = "cell_type",
        sample_col = "donor_id",
        min_cells = 3L)
    expect_type(mix, "list")
    expect_equal(length(mix$weights), 3L)
    expect_true(abs(sum(mix$weights) - 1) < 0.01)
    expect_equal(mix$celltype, "Mono")
})

test_that("run_sc_mixture errors on bad exposure", {
    scee <- .make_scee()
    expect_error(run_sc_mixture(scee,
        exposures = c("E1", "FAKE"),
        celltype = "Mono",
        celltype_col = "cell_type",
        sample_col = "donor_id"), "not found")
})

test_that("run_sc_mixture errors on too few donors", {
    scee <- .make_scee(n_donors = 2, n_cells = 20)
    expect_error(run_sc_mixture(scee,
        exposures = c("E1", "E2"),
        celltype = "Mono",
        celltype_col = "cell_type",
        sample_col = "donor_id",
        min_cells = 1L), "Need >= 5")
})

test_that("run_sc_mixture errors on bad celltype", {
    scee <- .make_scee()
    expect_error(run_sc_mixture(scee,
        exposures = c("E1", "E2"),
        celltype = "NONEXIST",
        celltype_col = "cell_type",
        sample_col = "donor_id"), "not found")
})

test_that("run_sc_mixture returns mixture_coef", {
    scee <- .make_scee()
    mix <- run_sc_mixture(scee,
        exposures = c("E1", "E2", "E3"),
        celltype = "NK",
        celltype_col = "cell_type",
        sample_col = "donor_id",
        min_cells = 3L)
    expect_true("mixture_coef" %in% names(mix))
    expect_true(is.numeric(mix$mixture_coef))
})

test_that("run_sc_mixture with target_genes", {
    scee <- .make_scee()
    mix <- run_sc_mixture(scee,
        exposures = c("E1", "E2"),
        celltype = "Mono",
        celltype_col = "cell_type",
        sample_col = "donor_id",
        target_genes = c("G1", "G2", "G3"),
        min_cells = 3L)
    expect_equal(mix$n_genes, 3L)
    expect_s4_class(mix$gene_results, "DataFrame")
})

test_that("run_sc_mixture with covariates", {
    scee <- .make_scee()
    mix <- run_sc_mixture(scee,
        exposures = c("E1", "E2"),
        celltype = "Mono",
        celltype_col = "cell_type",
        sample_col = "donor_id",
        covariates = "E3",
        min_cells = 3L)
    expect_type(mix, "list")
    expect_equal(length(mix$weights), 2L)
})

test_that("run_sc_mixture method field is quantile_linear", {
    scee <- .make_scee()
    mix <- run_sc_mixture(scee,
        exposures = c("E1", "E2"),
        celltype = "Mono",
        celltype_col = "cell_type",
        sample_col = "donor_id",
        min_cells = 3L)
    expect_equal(mix$method, "quantile_linear")
})
