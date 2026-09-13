# tests/testthat/test-dose-response-models.R

.donor_fixture <- function(n_donors = 16, curvature = 0, slope_b = 0,
                           seed = 21) {
    set.seed(seed)
    donors <- sprintf("D%02d", seq_len(n_donors))
    cells_per_type <- 25
    donor <- rep(rep(donors, each = cells_per_type), times = 2)
    cell_type <- rep(c("A", "B"), each = n_donors * cells_per_type)
    exposure <- stats::setNames(seq(-1.5, 1.5, length.out = n_donors), donors)
    counts <- matrix(stats::rpois(600 * length(donor), 30), nrow = 600,
        dimnames = list(paste0("G", 1:600), paste0("c", seq_along(donor))))
    e <- exposure[donor]
    mean_g1 <- 30 * exp(0.6 * e + curvature * e^2 +
                        ifelse(cell_type == "B", slope_b * e, 0))
    mean_g2 <- 30 * exp(-0.6 * e)
    counts[1, ] <- stats::rpois(length(donor), mean_g1)
    counts[2, ] <- stats::rpois(length(donor), mean_g2)
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = counts),
        colData = S4Vectors::DataFrame(donor_id = donor,
                                       cell_type = cell_type))
    exp_mat <- cbind(E = exposure, age = seq(30, 60, length.out = n_donors))
    rownames(exp_mat) <- donors
    build_scee(sce, exp_mat, sample_col = "donor_id")
}

test_that("run_interaction_test uses covariates and detects slope differences", {
    scee <- .donor_fixture(slope_b = -1.2)
    expect_error(run_interaction_test(scee, exposure = "E",
        covariates = "nope"), "Covariate")
    res <- as.data.frame(run_interaction_test(scee, exposure = "E",
        target_genes = c("G1", "G2")))
    expect_lt(res$p_interaction[res$gene == "G1"], 0.001)
    expect_equal(res$df_interaction[res$gene == "G1"], 1L)
    adjusted <- as.data.frame(run_interaction_test(scee, exposure = "E",
        covariates = "age", target_genes = "G1"))
    expect_false(isTRUE(all.equal(adjusted$p_interaction,
                                  res$p_interaction[res$gene == "G1"])))
})

test_that("run_dose_response tests the pre-specified polynomial", {
    scee <- .donor_fixture(curvature = 0.8)
    res <- as.data.frame(run_dose_response(scee, exposure = "E",
        celltype = "A", target_genes = c("G1", "G2")))
    expect_true(all(res$nonlinear_model == "cubic"))
    expect_lt(res$p_nonlinear[res$gene == "G1"], 0.01)
})

test_that("run_dose_response_gam separates curvature from a linear trend", {
    skip_if_not_installed("mgcv")
    scee <- .donor_fixture(n_donors = 20, curvature = 0.8)
    res <- as.data.frame(run_dose_response_gam(scee, exposure = "E",
        celltype = "A", target_genes = c("G1", "G2"), loocv = FALSE))
    expect_lt(res$p_nonlinear[res$gene == "G1"], 0.01)
    expect_gt(res$p_nonlinear[res$gene == "G2"], 0.05)
})
