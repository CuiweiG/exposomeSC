library(SingleCellExperiment)
library(S4Vectors)

.make_scee_ix <- function() {
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

test_that("run_interaction_test returns DataFrame", {
    scee <- .make_scee_ix()
    ix <- run_interaction_test(scee,
        exposure = "E1",
        celltype_col = "cell_type",
        sample_col = "donor_id",
        min_cells = 3L, min_donors = 3L)
    expect_s4_class(ix, "DataFrame")
    expect_true("gene" %in% colnames(ix))
    expect_true("p_interaction" %in% colnames(ix))
    expect_true("padj_interaction" %in% colnames(ix))
    expect_true(all(ix$p_interaction >= 0 &
                    ix$p_interaction <= 1))
})

test_that("interaction test errors on single celltype", {
    set.seed(1)
    counts <- matrix(rpois(500, 8), nrow = 100,
        dimnames = list(paste0("G", 1:100),
                        paste0("c", 1:5)))
    sce <- SingleCellExperiment(
        assays = list(counts = counts),
        colData = DataFrame(
            cell_id = paste0("c", 1:5),
            donor_id = rep(paste0("D", 1:5), each = 1),
            cell_type = rep("A", 5)))
    exp_mat <- matrix(rnorm(5), nrow = 5,
        dimnames = list(paste0("D", 1:5), "E1"))
    scee <- build_scee(sce, exp_mat, sample_col = "donor_id")
    expect_error(run_interaction_test(scee,
        exposure = "E1",
        celltype_col = "cell_type",
        sample_col = "donor_id",
        min_cells = 1L, min_donors = 1L),
        "2 cell types")
})
