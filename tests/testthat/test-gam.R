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

    for (f in list.files("../../R", full.names = TRUE))
        source(f, local = TRUE)

    scee <- build_scee(sce, exp_mat, sample_col = "donor_id")

    dr <- run_dose_response_gam(scee, exposure = "E1",
        celltype = "Mono", celltype_col = "cell_type",
        sample_col = "donor_id", min_cells = 3L,
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

    for (f in list.files("../../R", full.names = TRUE))
        source(f, local = TRUE)

    scee <- build_scee(sce, exp_mat, sample_col = "donor_id")

    dr <- run_dose_response_gam(scee, exposure = "E1",
        celltype = "Mono", celltype_col = "cell_type",
        sample_col = "donor_id", min_cells = 3L,
        target_genes = paste0("G", 1:5), loocv = TRUE)

    expect_true(all(!is.na(dr$R2_loocv_linear)))
    ## LOOCV R² should be <= apparent R²
    expect_true(all(dr$R2_loocv_linear <= dr$R2_linear + 0.01))
})
