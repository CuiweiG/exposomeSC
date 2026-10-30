## The cell-to-donor mapping lives only in colData

.make_coldata_scee <- function(sample_col = "donor_id", n_donors = 3L,
                               cells_per_donor = 4L) {
    set.seed(11)
    donors <- paste0("D", seq_len(n_donors))
    n_cells <- n_donors * cells_per_donor
    cells <- paste0("c", seq_len(n_cells))
    counts <- matrix(stats::rpois(8 * n_cells, 20), nrow = 8,
        dimnames = list(paste0("G", seq_len(8)), cells))
    col_data <- S4Vectors::DataFrame(
        cell_id = cells,
        donor = rep(donors, each = cells_per_donor),
        cell_type = rep(c("A", "B"), length.out = n_cells))
    names(col_data)[2] <- sample_col
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = counts), colData = col_data)
    exposures <- cbind(E1 = seq_len(n_donors),
                       age = 40 + 10 * seq_len(n_donors))
    rownames(exposures) <- donors
    build_scee(sce, exposures, sample_col = sample_col)
}

test_that("build_scee stores the donor column, not a copy of the mapping", {
    scee <- .make_coldata_scee("subject")
    expect_identical(sampleColumn(scee), "subject")
    expect_identical(as.character(cellSampleMap(scee)$sample_id),
                     as.character(scee$subject))
})

test_that("cbind of cell subsets restores the original object", {
    scee <- .make_coldata_scee()
    first <- scee[, 1:4]
    second <- scee[, 5:12]
    expect_identical(rownames(exposureData(first)), "D1")

    combined <- cbind(first, second)
    expect_s4_class(combined, "SingleCellExposomeExperiment")
    expect_true(methods::validObject(combined))
    expect_identical(colnames(combined), colnames(scee))
    expect_identical(exposureData(combined), exposureData(scee))
    expect_identical(exposureInfo(combined), exposureInfo(scee))
    expect_identical(sampleColumn(combined), "donor_id")
    expect_identical(cellSampleMap(combined), cellSampleMap(scee))
    expect_identical(SummarizedExperiment::assay(combined, "counts"),
                     SummarizedExperiment::assay(scee, "counts"))
})

test_that("cbind accepts a shared donor only with identical exposures", {
    scee <- .make_coldata_scee()
    ## cells 1-4 and 5-6 share no donor; 5-6 and 7-8 are both D2
    combined <- cbind(scee[, 1:6], scee[, 7:8])
    expect_identical(rownames(exposureData(combined)), c("D1", "D2"))
    expect_true(methods::validObject(combined))

    other <- scee[, 7:8]
    changed <- exposureData(other)
    changed["D2", "E1"] <- 99
    exposureData(other) <- changed
    expect_error(cbind(scee[, 1:6], other),
                 "different exposure values")
})

test_that("cbind rejects objects that do not fit together", {
    scee <- .make_coldata_scee()
    ## mixing classes also makes S4 dispatch warn about the ambiguity
    suppressWarnings(expect_error(
        cbind(scee[, 1:4], as(scee[, 5:12], "SingleCellExperiment")),
        "combines SingleCellExposomeExperiment objects only"))

    renamed <- .make_coldata_scee("subject")
    colnames(renamed) <- paste0("k", seq_len(ncol(renamed)))
    expect_error(cbind(scee, renamed), "same sample column")

    fewer <- .make_coldata_scee()
    colnames(fewer) <- paste0("k", seq_len(ncol(fewer)))
    exposureInfo(fewer)$unit <- c("ug", "years")
    expect_error(cbind(scee, fewer), "same exposureInfo")
})

test_that("renaming cells keeps the object valid", {
    scee <- .make_coldata_scee()
    colnames(scee) <- paste0("cell_", seq_len(ncol(scee)))
    expect_true(methods::validObject(scee))
    map <- cellSampleMap(scee)
    expect_identical(as.character(map$cell_id), colnames(scee))
    expect_identical(as.character(map$sample_id),
                     rep(c("D1", "D2", "D3"), each = 4))
})

test_that("colData edits are checked against exposureData", {
    scee <- .make_coldata_scee()

    ## moving cells between known donors is allowed, and the map follows
    moved <- scee
    moved$donor_id <- rep(c("D3", "D2", "D1"), each = 4)
    expect_identical(as.character(cellSampleMap(moved)$sample_id),
                     rep(c("D3", "D2", "D1"), each = 4))

    ## other columns can be edited freely
    edited <- scee
    edited$batch <- rep(1:2, 6)
    expect_true(methods::validObject(edited))

    expect_error(scee$donor_id <- rep(c("D1", "D2", "D9"), each = 4),
                 "not in exposureData")
    expect_error(scee$donor_id <- rep(c("D1", "D2"), each = 6),
                 "without cells")
    expect_error(scee$donor_id <- NULL, "not a column of colData")
    expect_error(scee[["donor_id"]] <- rep(NA_character_, 12),
                 "non-missing")
    expect_error(
        SummarizedExperiment::colData(scee) <-
            S4Vectors::DataFrame(cell_type = scee$cell_type),
        "not a column of colData")
    ## the failed edits left the object unchanged
    expect_identical(cellSampleMap(scee), cellSampleMap(.make_coldata_scee()))
})

test_that("the donor column can be switched to an equivalent column", {
    scee <- .make_coldata_scee()
    scee$subject <- scee$donor_id
    sampleColumn(scee) <- "subject"
    expect_identical(sampleColumn(scee), "subject")
    expect_true(methods::validObject(scee))

    scee$other_ids <- paste0("X", scee$donor_id)
    expect_error(sampleColumn(scee) <- "other_ids", "not in exposureData")
    expect_error(sampleColumn(scee) <- "absent", "not a column of colData")
    expect_error(sampleColumn(scee) <- c("a", "b"), "single column name")
})

test_that("subsetting keeps the donor column and the exposure rows aligned", {
    scee <- .make_coldata_scee("subject")
    sub <- scee[, c(12, 1, 5)]
    expect_identical(sampleColumn(sub), "subject")
    expect_identical(rownames(exposureData(sub)), c("D3", "D1", "D2"))
    expect_identical(as.character(cellSampleMap(sub)$sample_id),
                     c("D3", "D1", "D2"))
    genes <- scee[1:3, ]
    expect_identical(exposureData(genes), exposureData(scee))
})

test_that("analysis functions read the donor from the stored column", {
    by_donor_id <- .make_coldata_scee("donor_id", n_donors = 8L,
                                      cells_per_donor = 10L)
    by_subject <- .make_coldata_scee("subject", n_donors = 8L,
                                     cells_per_donor = 10L)
    fit_a <- run_sc_exwas(by_donor_id, exposure = "E1", min_cells = 3L,
        min_donors = 5L, filter_genes = FALSE)
    fit_b <- run_sc_exwas(by_subject, exposure = "E1", min_cells = 3L,
        min_donors = 5L, filter_genes = FALSE)
    expect_gt(nrow(fit_a), 0L)
    expect_equal(as.data.frame(fit_a), as.data.frame(fit_b))

    expect_error(run_sc_exwas(by_subject, exposure = "E1",
        sample_col = "subject"), "Unused argument")
    expect_error(run_dose_response(by_subject, exposure = "E1",
        celltype = "A", sample_col = "subject"), "unused argument")
})
