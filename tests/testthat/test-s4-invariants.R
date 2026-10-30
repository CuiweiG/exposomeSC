.make_scee_invariant <- function() {
    set.seed(42)
    counts <- matrix(
        stats::rpois(500, 10),
        nrow = 50,
        dimnames = list(paste0("G", 1:50), paste0("c", 1:10))
    )
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = counts),
        colData = S4Vectors::DataFrame(
            cell_id = paste0("c", 1:10),
            donor_id = rep(c("D1", "D2"), each = 5),
            cell_type = rep(c("A", "B"), 5)
        )
    )
    exposures <- matrix(
        stats::rnorm(4),
        nrow = 2,
        dimnames = list(c("D1", "D2"), c("E1", "E2"))
    )
    build_scee(sce, exposures, sample_col = "donor_id")
}

test_that("the cell-to-sample map is derived from colData", {
    scee <- .make_scee_invariant()
    expect_identical(sampleColumn(scee), "donor_id")
    map <- cellSampleMap(scee)
    expect_identical(as.character(map$cell_id), colnames(scee))
    expect_identical(as.character(map$sample_id),
                     as.character(scee$donor_id))

    ## cells need unique, non-missing names
    expect_error(colnames(scee) <- rep("c1", ncol(scee)), "unique")
    expect_error(colnames(scee) <- NULL, "non-missing names")
})

test_that("cell subsetting keeps cell order and rejects duplicates", {
    scee <- .make_scee_invariant()
    selected <- c(8L, 2L, 7L, 1L)
    subset <- scee[, selected]
    expect_identical(as.character(cellSampleMap(subset)$cell_id), colnames(subset))
    expect_identical(
        rownames(exposureData(subset)),
        unique(as.character(cellSampleMap(subset)$sample_id))
    )
    expect_true(methods::validObject(subset))

    expect_error(
        scee[, c(1L, 1L, 3L)],
        "Duplicated cell selection"
    )

    empty <- scee[, FALSE]
    expect_equal(ncol(empty), 0L)
    expect_equal(nrow(cellSampleMap(empty)), 0L)
    expect_equal(nrow(exposureData(empty)), 0L)
    expect_true(methods::validObject(empty))
})

test_that("exposure metadata names must match exposure columns", {
    scee <- .make_scee_invariant()
    info <- exposureInfo(scee)
    info$exposure <- rev(info$exposure)
    expect_error(exposureInfo(scee) <- info, "must match")
})

test_that("exposure data enforce numeric finite and exact sample support", {
    scee <- .make_scee_invariant()

    character_data <- matrix(
        c("0", "1", "2", "3"),
        nrow = 2,
        dimnames = dimnames(exposureData(scee))
    )
    expect_error(
        exposureData(scee) <- character_data,
        "numeric matrix"
    )

    infinite_data <- exposureData(scee)
    infinite_data[1, 1] <- Inf
    expect_error(exposureData(scee) <- infinite_data, "not NaN or infinite")

    scee <- .make_scee_invariant()
    extra_data <- rbind(
        exposureData(scee),
        D3 = c(E1 = 0, E2 = 0)
    )
    expect_error(exposureData(scee) <- extra_data, "without cells")
})

test_that("build_scee rejects malformed matrices and extra samples", {
    scee <- .make_scee_invariant()
    sce <- as(scee, "SingleCellExperiment")
    exposure <- exposureData(scee)

    character_data <- matrix(
        as.character(exposure),
        nrow = nrow(exposure),
        dimnames = dimnames(exposure)
    )
    expect_error(
        build_scee(sce, character_data, sample_col = "donor_id"),
        "numeric matrix"
    )

    extra_data <- rbind(exposure, D3 = c(E1 = 0, E2 = 0))
    expect_error(
        build_scee(sce, extra_data, sample_col = "donor_id"),
        "not represented"
    )

    bad_info <- S4Vectors::DataFrame(exposure = rev(colnames(exposure)))
    expect_error(
        build_scee(
            sce,
            exposure,
            sample_col = "donor_id",
            exposure_info = bad_info
        ),
        "must match"
    )
})

test_that("CelltypeNetworkResult reports a dimension mismatch as one message", {
    message_text <- tryCatch(
        CelltypeNetworkResult(precision = diag(2),
            adjacency = matrix(0L, 3L, 3L),
            node_info = S4Vectors::DataFrame(feature = c("G1", "G2")),
            celltype = "Mono"),
        error = conditionMessage)
    expect_match(message_text,
        "adjacency_matrix dimensions must match precision_matrix",
        fixed = TRUE)
    expect_false(grepl("2: precision_matrix", message_text, fixed = TRUE))
})
