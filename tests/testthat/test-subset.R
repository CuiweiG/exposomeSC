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

test_that("subsetting by cells drops donors without cells", {
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

test_that("rbind needs the same exposure slots in every object", {
    set.seed(3)
    donor <- rep(c("D1", "D2", "D3"), each = 4)
    make_sce <- function(genes) {
        counts <- matrix(stats::rpois(length(genes) * 12, 5),
            nrow = length(genes),
            dimnames = list(genes, paste0("c", 1:12)))
        SingleCellExperiment::SingleCellExperiment(
            assays = list(counts = counts),
            colData = S4Vectors::DataFrame(donor_id = donor,
                cell_type = "Mono"))
    }
    exposures <- matrix(c(0.1, 0.5, 0.9), ncol = 1,
        dimnames = list(c("D1", "D2", "D3"), "PM2.5"))
    a <- build_scee(make_sce(c("G1", "G2")), exposures, "donor_id")
    b <- build_scee(make_sce(c("G3", "G4")), exposures, "donor_id")
    combined <- rbind(a, b)
    expect_s4_class(combined, "SingleCellExposomeExperiment")
    expect_identical(rownames(combined), c("G1", "G2", "G3", "G4"))
    expect_identical(exposureData(combined), exposureData(a))
    shifted <- exposures
    shifted["D2", "PM2.5"] <- 2
    c2 <- build_scee(make_sce(c("G3", "G4")), shifted, "donor_id")
    expect_error(rbind(a, c2), "same exposure values for each donor")
    ## the same values in another row order or storage mode combine
    reordered <- build_scee(make_sce(c("G3", "G4")),
        exposures[3:1, , drop = FALSE], "donor_id")
    expect_identical(exposureData(rbind(a, reordered)), exposureData(a))
    whole <- matrix(1:3, ncol = 1,
        dimnames = list(c("D1", "D2", "D3"), "PM2.5"))
    as_integer <- build_scee(make_sce(c("G1", "G2")), whole, "donor_id")
    as_double <- build_scee(make_sce(c("G3", "G4")), whole * 1, "donor_id")
    expect_s4_class(rbind(as_integer, as_double),
        "SingleCellExposomeExperiment")
    ## the sample column and exposureInfo must agree, and every argument
    ## must be an exposome object
    sce2 <- make_sce(c("G3", "G4"))
    sce2$donor2 <- sce2$donor_id
    expect_error(rbind(a, build_scee(sce2, exposures, "donor2")),
        "same sample column")
    relabelled <- b
    info <- exposureInfo(relabelled)
    info$family <- "gaussian"
    exposureInfo(relabelled) <- info
    expect_error(rbind(a, relabelled), "same sample column and exposureInfo")
    expect_error(
        suppressWarnings(rbind(a, methods::as(b, "SingleCellExperiment"))),
        "objects only")
})
