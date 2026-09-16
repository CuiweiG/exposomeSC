## Tests for [ subsetting

library(SingleCellExperiment)
library(S4Vectors)

.make_scee_sub <- function() {
    set.seed(42)
    counts <- matrix(rpois(500, 10), nrow = 50,
        dimnames = list(paste0("G", 1:50), paste0("c", 1:10)))
    sce <- SingleCellExperiment(
        assays = list(counts = counts),
        colData = DataFrame(
            cell_id = paste0("c", 1:10),
            donor_id = rep(c("D1", "D2"), each = 5),
            cell_type = rep(c("A", "B"), 5)))
    exp_mat <- matrix(rnorm(4), nrow = 2,
        dimnames = list(c("D1", "D2"), c("E1", "E2")))
    build_scee(sce, exp_mat, sample_col = "donor_id")
}

test_that("subsetting by genes preserves exposure data", {
    scee <- .make_scee_sub()
    sub <- scee[1:10, ]
    expect_s4_class(sub, "SingleCellExposomeExperiment")
    expect_equal(nrow(sub), 10L)
    expect_equal(ncol(sub), 10L)
    expect_equal(nrow(exposureData(sub)), 2L)
    expect_equal(nrow(cellSampleMap(sub)), 10L)
})

test_that("subsetting by cells updates sampleMap", {
    scee <- .make_scee_sub()
    ## Keep only first 5 cells (D1 only)
    sub <- scee[, 1:5]
    expect_s4_class(sub, "SingleCellExposomeExperiment")
    expect_equal(ncol(sub), 5L)
    expect_equal(nrow(cellSampleMap(sub)), 5L)
    ## Only D1 should remain in exposureData
    expect_equal(nrow(exposureData(sub)), 1L)
    expect_equal(rownames(exposureData(sub)), "D1")
})

test_that("subsetting by both genes and cells works", {
    scee <- .make_scee_sub()
    sub <- scee[1:20, 1:5]
    expect_equal(nrow(sub), 20L)
    expect_equal(ncol(sub), 5L)
    expect_equal(nrow(exposureData(sub)), 1L)
})
