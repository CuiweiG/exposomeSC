## Tests for S4 validity constraints

library(SingleCellExperiment)
library(S4Vectors)

test_that("validity catches mismatched exposureInfo nrow", {
    set.seed(1)
    counts <- matrix(rpois(100, 5), nrow = 10,
        dimnames = list(paste0("G", 1:10), paste0("c", 1:10)))
    sce <- SingleCellExperiment(
        assays = list(counts = counts),
        colData = DataFrame(
            cell_id = paste0("c", 1:10),
            donor = rep(c("D1", "D2"), each = 5),
            ct = rep("A", 10)))
    exp_mat <- matrix(rnorm(4), nrow = 2,
        dimnames = list(c("D1", "D2"), c("E1", "E2")))
    scee <- build_scee(sce, exp_mat, sample_col = "donor")

    ## Try to set exposureInfo with wrong number of rows
    bad_info <- DataFrame(exposure = "X", family = "Y",
                          unit = "Z")
    expect_error(exposureInfo(scee) <- bad_info)
})

test_that("validity catches sampleMap sample_id not in exposureData", {
    set.seed(1)
    counts <- matrix(rpois(100, 5), nrow = 10,
        dimnames = list(paste0("G", 1:10), paste0("c", 1:10)))
    sce <- SingleCellExperiment(
        assays = list(counts = counts),
        colData = DataFrame(
            cell_id = paste0("c", 1:10),
            donor = rep(c("D1", "D2"), each = 5),
            ct = rep("A", 10)))
    exp_mat <- matrix(rnorm(4), nrow = 2,
        dimnames = list(c("D1", "D2"), c("E1", "E2")))
    scee <- build_scee(sce, exp_mat, sample_col = "donor")

    ## Replace exposureData with different sample IDs
    new_mat <- matrix(rnorm(4), nrow = 2,
        dimnames = list(c("X1", "X2"), c("E1", "E2")))
    expect_error(exposureData(scee) <- new_mat)
})
