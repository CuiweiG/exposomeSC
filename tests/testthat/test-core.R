library(SingleCellExperiment)
library(S4Vectors)

.make_test_scee <- function(n_donors = 5, n_cells_per = 20,
                             n_genes = 50) {
    set.seed(42)
    n_total <- n_donors * n_cells_per
    counts <- matrix(rpois(n_genes * n_total, lambda = 10),
                      nrow = n_genes)
    rownames(counts) <- paste0("Gene", seq_len(n_genes))
    colnames(counts) <- paste0("cell_", seq_len(n_total))
    sce <- SingleCellExperiment(
        assays = list(counts = counts),
        colData = DataFrame(
            cell_id = colnames(counts),
            donor_id = rep(paste0("D", seq_len(n_donors)),
                           each = n_cells_per),
            cell_type = rep(c("Mono", "NK"),
                            length.out = n_total)))

    exp_mat <- matrix(rnorm(n_donors * 3), nrow = n_donors,
        dimnames = list(paste0("D", seq_len(n_donors)),
                        c("PM2.5", "Pb", "BPA")))
    build_scee(sce, exp_mat, sample_col = "donor_id")
}

## ---- build_scee ---------------------------------------------------

test_that("build_scee constructs valid SCEE", {
    scee <- .make_test_scee()
    expect_s4_class(scee, "SingleCellExposomeExperiment")
    expect_s4_class(scee, "SingleCellExperiment")
    expect_equal(ncol(exposureData(scee)), 3L)
    expect_equal(nrow(exposureData(scee)), 5L)
    expect_equal(nrow(cellSampleMap(scee)), 100L)
    expect_true(validObject(scee))
})

test_that("accessors work", {
    scee <- .make_test_scee()
    expect_equal(exposureVariables(scee), c("PM2.5", "Pb", "BPA"))
    expect_s4_class(exposureInfo(scee), "DataFrame")
    expect_s4_class(cellSampleMap(scee), "DataFrame")
})

test_that("show method works", {
    scee <- .make_test_scee()
    expect_message(show(scee), "exposureData")
    expect_message(show(scee), "sampleMap")
})

test_that("build_scee errors on missing sample_col", {
    set.seed(1)
    sce <- SingleCellExperiment(
        assays = list(counts = matrix(1L, 5, 10)),
        colData = DataFrame(x = rep("A", 10)))
    exp_mat <- matrix(1, nrow = 1,
        dimnames = list("A", "E1"))
    expect_error(build_scee(sce, exp_mat,
        sample_col = "nonexistent"), "not found")
})

test_that("build_scee errors on mismatched samples", {
    set.seed(1)
    sce <- SingleCellExperiment(
        assays = list(counts = matrix(1L, 5, 10)),
        colData = DataFrame(donor = rep("A", 10)))
    exp_mat <- matrix(1, nrow = 1,
        dimnames = list("B", "E1"))
    expect_error(build_scee(sce, exp_mat,
        sample_col = "donor"), "not found in")
})

test_that("build_scee errors on missing rownames", {
    set.seed(1)
    sce <- SingleCellExperiment(
        assays = list(counts = matrix(1L, 5, 10)),
        colData = DataFrame(donor = rep("A", 10)))
    exp_mat <- matrix(1, nrow = 1, ncol = 1)
    colnames(exp_mat) <- "E1"
    expect_error(build_scee(sce, exp_mat,
        sample_col = "donor"), "row names")
})

test_that("build_scee accepts custom exposure_info", {
    set.seed(1)
    counts <- matrix(rpois(100, 5), nrow = 10,
        dimnames = list(paste0("G", 1:10), paste0("c", 1:10)))
    sce <- SingleCellExperiment(
        assays = list(counts = counts),
        colData = DataFrame(
            cell_id = paste0("c", 1:10),
            donor = rep(c("D1", "D2"), each = 5),
            ct = rep(c("A", "B"), 5)))
    exp_mat <- matrix(rnorm(4), nrow = 2,
        dimnames = list(c("D1", "D2"), c("E1", "E2")))
    ei <- DataFrame(exposure = c("E1", "E2"),
                    family = c("air", "metal"),
                    unit = c("ug/m3", "ug/dL"))
    scee <- build_scee(sce, exp_mat, sample_col = "donor",
                        exposure_info = ei)
    expect_equal(nrow(exposureInfo(scee)), 2L)
    expect_equal(exposureInfo(scee)$family, c("air", "metal"))
})

## ---- Setter methods -----------------------------------------------

test_that("exposureData<- works", {
    scee <- .make_test_scee()
    new_mat <- exposureData(scee) * 2
    exposureData(scee) <- new_mat
    expect_equal(exposureData(scee), new_mat)
})

test_that("exposureInfo<- works", {
    scee <- .make_test_scee()
    new_info <- DataFrame(
        exposure = c("PM2.5", "Pb", "BPA"),
        family = c("air", "metal", "chemical"),
        unit = c("ug/m3", "ug/dL", "ng/mL"))
    exposureInfo(scee) <- new_info
    expect_equal(exposureInfo(scee)$family,
                 c("air", "metal", "chemical"))
})

test_that("exposureData<- rejects non-matrix", {
    scee <- .make_test_scee()
    expect_error(exposureData(scee) <- data.frame(x = 1))
})

## ---- Validity -----------------------------------------------------

test_that("validity catches mismatched exposureInfo rows", {
    scee <- .make_test_scee()
    bad_info <- DataFrame(exposure = "X", family = "Y",
                          unit = "Z")
    expect_error(exposureInfo(scee) <- bad_info)
})

## ---- run_sc_exwas -------------------------------------------------

test_that("run_sc_exwas pseudobulk works", {
    scee <- .make_test_scee()
    result <- run_sc_exwas(scee,
        exposure = "PM2.5",
        celltype_col = "cell_type",
        sample_col = "donor_id",
        min_cells = 5L,
        min_donors = 3L)
    expect_s4_class(result, "DataFrame")
    expect_true("gene" %in% colnames(result))
    expect_true("celltype" %in% colnames(result))
    expect_true("log2FC" %in% colnames(result))
    expect_true("se" %in% colnames(result))
    expect_true("statistic" %in% colnames(result))
    expect_true("pvalue" %in% colnames(result))
    expect_true("padj_global" %in% colnames(result))
    expect_true("n_donors" %in% colnames(result))
    expect_true(all(is.na(result$se)))
    expect_identical(
        S4Vectors::metadata(result)$parameters$statistic_type,
        "signed_sqrt_qlf"
    )
})

test_that("run_sc_exwas errors on bad exposure name", {
    scee <- .make_test_scee()
    expect_error(run_sc_exwas(scee,
        exposure = "NONEXISTENT",
        celltype_col = "cell_type",
        sample_col = "donor_id"), "not found")
})

test_that("run_sc_exwas respects celltypes parameter", {
    scee <- .make_test_scee()
    result <- run_sc_exwas(scee,
        exposure = "PM2.5",
        celltype_col = "cell_type",
        celltypes = "Mono",
        sample_col = "donor_id",
        min_cells = 5L,
        min_donors = 3L)
    expect_true(all(result$celltype == "Mono"))
})

test_that("run_sc_exwas with covariates", {
    scee <- .make_test_scee(n_donors = 8)
    result <- run_sc_exwas(scee,
        exposure = "PM2.5",
        celltype_col = "cell_type",
        sample_col = "donor_id",
        covariates = "Pb",
        min_cells = 5L,
        min_donors = 3L)
    expect_s4_class(result, "DataFrame")
    expect_true(nrow(result) > 0)
})

test_that("run_sc_exwas rejects an underdetermined donor design", {
    scee <- .make_test_scee(n_donors = 5, n_cells_per = 10)
    observed_warnings <- character()
    result <- withCallingHandlers(
        run_sc_exwas(scee,
            exposure = "PM2.5",
            celltype_col = "cell_type",
            sample_col = "donor_id",
            covariates = c("Pb", "BPA"),
            min_cells = 3L,
            min_donors = 3L),
        warning = function(condition) {
            observed_warnings <<- c(
                observed_warnings,
                conditionMessage(condition)
            )
            invokeRestart("muffleWarning")
        }
    )
    expect_equal(nrow(result), 0L)
    expect_length(observed_warnings, 2L)
    expect_true(all(grepl(
        "insufficient residual degrees of freedom",
        observed_warnings,
        fixed = TRUE
    )))
})

test_that("run_sc_exwas returns empty DataFrame if no valid celltypes", {
    scee <- .make_test_scee(n_donors = 2, n_cells_per = 5)
    result <- suppressWarnings(run_sc_exwas(scee,
        exposure = "PM2.5",
        celltype_col = "cell_type",
        sample_col = "donor_id",
        min_cells = 5L,
        min_donors = 5L))
    expect_s4_class(result, "DataFrame")
    expect_equal(nrow(result), 0L)
})

test_that("run_sc_exwas filter_genes=FALSE allows all genes", {
    scee <- .make_test_scee()
    r1 <- run_sc_exwas(scee,
        exposure = "PM2.5",
        celltype_col = "cell_type",
        sample_col = "donor_id",
        min_cells = 5L, min_donors = 3L,
        filter_genes = TRUE)
    r2 <- run_sc_exwas(scee,
        exposure = "PM2.5",
        celltype_col = "cell_type",
        sample_col = "donor_id",
        min_cells = 5L, min_donors = 3L,
        filter_genes = FALSE)
    expect_true(nrow(r2) >= nrow(r1))
})
