# tests/testthat/test-gam.R

test_that("run_dose_response_gam works", {
    skip_if_not_installed("mgcv")

    set.seed(42)
    counts <- matrix(rpois(10000, 10), nrow = 100,
        dimnames = list(paste0("G", 1:100), paste0("c", 1:100)))
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = counts),
        colData = S4Vectors::DataFrame(
            cell_id = paste0("c", 1:100),
            donor_id = rep(paste0("D", 1:10), each = 10),
            cell_type = rep(c("Mono", "NK"), 50)))
    exp_mat <- matrix(rnorm(20), nrow = 10,
        dimnames = list(paste0("D", 1:10), c("E1", "E2")))

    scee <- build_scee(sce, exp_mat, sample_col = "donor_id")

    dr <- run_dose_response_gam(scee, exposure = "E1",
        celltype = "Mono", celltype_col = "cell_type",
        min_cells = 3L,
        target_genes = paste0("G", 1:10))

    expect_s4_class(dr, "DataFrame")
    expect_true(nrow(dr) > 0)
    expect_true("edf" %in% colnames(dr))
    expect_true("R2_loocv_linear" %in% colnames(dr))
    expect_true("R2_loocv_gam" %in% colnames(dr))
    expect_true("deviance_explained" %in% colnames(dr))
})

test_that("GAM LOOCV R² is computed", {
    skip_if_not_installed("mgcv")

    set.seed(42)
    counts <- matrix(rpois(10000, 10), nrow = 100,
        dimnames = list(paste0("G", 1:100), paste0("c", 1:100)))
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = counts),
        colData = S4Vectors::DataFrame(
            cell_id = paste0("c", 1:100),
            donor_id = rep(paste0("D", 1:10), each = 10),
            cell_type = rep(c("Mono", "NK"), 50)))
    exp_mat <- matrix(rnorm(20), nrow = 10,
        dimnames = list(paste0("D", 1:10), c("E1", "E2")))

    scee <- build_scee(sce, exp_mat, sample_col = "donor_id")

    dr <- run_dose_response_gam(scee, exposure = "E1",
        celltype = "Mono", celltype_col = "cell_type",
        min_cells = 3L,
        target_genes = paste0("G", 1:5), loocv = TRUE)

    expect_true(all(!is.na(dr$R2_loocv_linear)))
    ## LOOCV R² should be <= apparent R²
    expect_true(all(dr$R2_loocv_linear <= dr$R2_linear + 0.01))
})

test_that("run_dose_response_gam leaves donors without a dose out", {
    skip_if_not_installed("mgcv")
    set.seed(1)
    donor <- rep(paste0("D", 1:8), each = 10)
    counts <- matrix(rpois(50 * 80, 8), nrow = 50,
        dimnames = list(paste0("G", 1:50), paste0("c", 1:80)))
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = counts),
        colData = S4Vectors::DataFrame(donor_id = donor, cell_type = "Mono"))
    dose <- matrix(c(1, 2, NA, 4:8), ncol = 1,
        dimnames = list(paste0("D", 1:8), "E1"))
    scee <- build_scee(sce, dose, sample_col = "donor_id")
    genes <- paste0("G", 1:3)
    expect_message(
        with_missing <- run_dose_response_gam(scee, exposure = "E1",
            celltype = "Mono", min_cells = 3L, target_genes = genes),
        "Leaving 1 donor\\(s\\) with a missing 'E1' value")
    expect_true(all(with_missing$n_donors == 7L))
    expect_false(anyNA(with_missing$R2_loocv_linear))
    kept <- donor != "D3"
    without <- run_dose_response_gam(
        build_scee(sce[, kept], dose[-3, , drop = FALSE],
            sample_col = "donor_id"),
        exposure = "E1", celltype = "Mono", min_cells = 3L,
        target_genes = genes)
    expect_equal(with_missing, without)
})

test_that("run_dose_response_gam checks k against the exposure values", {
    skip_if_not_installed("mgcv")
    set.seed(1)
    donors <- paste0("D", 1:10)
    counts <- matrix(rpois(50 * 100, 8), nrow = 50,
        dimnames = list(paste0("G", 1:50), paste0("c", 1:100)))
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = counts),
        colData = S4Vectors::DataFrame(donor_id = rep(donors, each = 10),
            cell_type = "Mono"))
    gam_for <- function(values, ...) {
        dose <- matrix(values, ncol = 1, dimnames = list(donors, "dose"))
        run_dose_response_gam(build_scee(sce, dose, sample_col = "donor_id"),
            exposure = "dose", celltype = "Mono", min_cells = 3L,
            target_genes = c("G1", "G2"), ...)
    }
    expect_error(gam_for(rep(c(0, 1), each = 5)),
        "at least 4 distinct exposure values; found 2")
    ## with three values k would be 3, for which the test of nonlinearity
    ## cannot be fitted
    expect_error(gam_for(rep(0:2, length.out = 10)),
        "at least 4 distinct exposure values; found 3")
    graded <- gam_for(rep(1:4, length.out = 10))
    expect_identical(nrow(graded), 2L)
    expect_false(anyNA(graded$p_nonlinear))
    expect_false(anyNA(graded$R2_loocv_gam))
    ## one donor alone at the highest dose: its refit sees three values
    lone <- gam_for(c(0, 0, 0, 1, 1, 1, 2, 2, 2, 5))
    expect_identical(nrow(lone), 2L)
    expect_false(anyNA(lone$R2_loocv_gam))
    expect_error(gam_for(rnorm(10), k = 20L),
        "exceeds the number of distinct exposure values")
    expect_error(gam_for(rnorm(10), k = 2.5),
        "k must be NULL or one whole number >= 4")
    expect_error(gam_for(rnorm(10), k = 3L),
        "k must be NULL or one whole number >= 4")
})
