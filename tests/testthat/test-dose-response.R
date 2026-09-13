library(SingleCellExperiment)
library(S4Vectors)

.make_scee_dr <- function() {
    set.seed(42)
    counts <- matrix(rpois(5000, 8), nrow = 100,
        dimnames = list(paste0("G", 1:100),
                        paste0("c", 1:50)))
    sce <- SingleCellExperiment(
        assays = list(counts = counts),
        colData = DataFrame(
            cell_id = paste0("c", 1:50),
            donor_id = rep(paste0("D", 1:10), each = 5),
            cell_type = rep(c("Mono", "NK"), 25)))
    exp_mat <- matrix(rnorm(20), nrow = 10,
        dimnames = list(paste0("D", 1:10), c("E1", "E2")))
    build_scee(sce, exp_mat, sample_col = "donor_id")
}

test_that("run_dose_response returns valid DataFrame", {
    scee <- .make_scee_dr()
    dr <- run_dose_response(scee,
        exposure = "E1", celltype = "Mono",
        celltype_col = "cell_type",
        sample_col = "donor_id",
        min_cells = 3L, max_degree = 2L)
    expect_s4_class(dr, "DataFrame")
    expect_true("gene" %in% colnames(dr))
    expect_true("best_model" %in% colnames(dr))
    expect_true("AIC_linear" %in% colnames(dr))
    expect_true("p_nonlinear" %in% colnames(dr))
    expect_true(all(dr$best_model %in%
        c("linear", "quadratic", "cubic")))
})

test_that("run_dose_response errors on bad exposure", {
    scee <- .make_scee_dr()
    expect_error(run_dose_response(scee,
        exposure = "FAKE", celltype = "Mono",
        celltype_col = "cell_type",
        sample_col = "donor_id"), "not found")
})

test_that("run_dose_response errors on too few donors", {
    set.seed(1)
    counts <- matrix(rpois(400, 8), nrow = 100,
        dimnames = list(paste0("G", 1:100),
                        paste0("c", 1:4)))
    sce <- SingleCellExperiment(
        assays = list(counts = counts),
        colData = DataFrame(
            cell_id = paste0("c", 1:4),
            donor_id = rep(c("D1", "D2"), each = 2),
            cell_type = rep("A", 4)))
    exp_mat <- matrix(rnorm(2), nrow = 2,
        dimnames = list(c("D1", "D2"), "E1"))
    scee <- build_scee(sce, exp_mat, sample_col = "donor_id")
    expect_error(run_dose_response(scee,
        exposure = "E1", celltype = "A",
        celltype_col = "cell_type",
        sample_col = "donor_id",
        min_cells = 1L, max_degree = 2L), "donors")
})
