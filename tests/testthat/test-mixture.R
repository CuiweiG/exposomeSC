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

test_that("run_sc_mixture weights average absolute per-gene coefficients", {
    scee <- .make_scee()
    genes <- paste0("G", 1:3)
    mix <- run_sc_mixture(scee, exposures = c("E1", "E2"), celltype = "Mono",
        celltype_col = "cell_type", sample_col = "donor_id",
        target_genes = genes, min_cells = 3L)
    cd <- SummarizedExperiment::colData(scee)
    pb <- exposomeSC:::.pseudobulk_aggregate(
        SummarizedExperiment::assay(scee, "counts"),
        as.character(cd$donor_id), as.character(cd$cell_type), "Mono",
        min_cells = 3L)
    log_cpm <- exposomeSC:::.log_cpm(pb$pb_mat)
    e <- exposureData(scee)[pb$valid_donors, c("E1", "E2")]
    score <- function(v) as.integer(cut(v, unique(stats::quantile(v,
        seq(0, 1, length.out = 5))), include.lowest = TRUE))
    q_df <- data.frame(E1 = score(e[, "E1"]), E2 = score(e[, "E2"]))
    coefs <- t(vapply(genes, function(g)
        stats::coef(lm(log_cpm[g, pb$valid_donors] ~ E1 + E2, data = q_df))[
            c("E1", "E2")], numeric(2)))
    expected <- colMeans(abs(coefs))
    expect_equal(unname(mix$weights), unname(expected / sum(expected)))
})

test_that("run_sc_mixture handles tied exposures and rejects unknown covariates", {
    scee <- .make_scee()
    exposure_matrix <- exposureData(scee)
    exposure_matrix[1:6, "E3"] <- 0.01
    exposureData(scee) <- exposure_matrix
    mix <- run_sc_mixture(scee, exposures = c("E1", "E3"), celltype = "Mono",
        celltype_col = "cell_type", sample_col = "donor_id", min_cells = 3L)
    expect_equal(sum(mix$weights), 1)
    expect_error(run_sc_mixture(scee, exposures = c("E1", "E2"),
        celltype = "Mono", celltype_col = "cell_type",
        sample_col = "donor_id", covariates = "nope", min_cells = 3L),
        "Covariates not found")
})
