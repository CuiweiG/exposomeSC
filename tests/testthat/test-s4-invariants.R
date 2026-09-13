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

test_that("sampleMap requires a complete ordered one-to-one cell mapping", {
    scee <- .make_scee_invariant()

    missing_columns <- S4Vectors::DataFrame(
        wrong_cell = colnames(scee),
        wrong_sample = rep(c("D1", "D2"), each = 5)
    )
    methods::slot(scee, "sampleMap") <- missing_columns
    expect_error(methods::validObject(scee), "missing required column")

    scee <- .make_scee_invariant()
    duplicate_map <- sampleMap(scee)
    duplicate_map$cell_id[[2]] <- duplicate_map$cell_id[[1]]
    methods::slot(scee, "sampleMap") <- duplicate_map
    expect_error(methods::validObject(scee), "must be unique")

    scee <- .make_scee_invariant()
    reversed_map <- sampleMap(scee)[
        rev(seq_len(nrow(sampleMap(scee)))),
        ,
        drop = FALSE
    ]
    methods::slot(scee, "sampleMap") <- reversed_map
    expect_error(methods::validObject(scee), "match colnames")
})

test_that("cell subsetting preserves sampleMap order and rejects duplicates", {
    scee <- .make_scee_invariant()
    selected <- c(8L, 2L, 7L, 1L)
    subset <- scee[, selected]
    expect_identical(as.character(sampleMap(subset)$cell_id), colnames(subset))
    expect_identical(
        rownames(exposureData(subset)),
        unique(as.character(sampleMap(subset)$sample_id))
    )
    expect_true(methods::validObject(subset))

    expect_error(
        scee[, c(1L, 1L, 3L)],
        "Duplicated cell selection"
    )

    empty <- scee[, FALSE]
    expect_equal(ncol(empty), 0L)
    expect_equal(nrow(sampleMap(empty)), 0L)
    expect_equal(nrow(exposureData(empty)), 0L)
    expect_true(methods::validObject(empty))
})

test_that("exposure metadata names must match exposure columns", {
    scee <- .make_scee_invariant()
    info <- exposureInfo(scee)
    info$exposure <- rev(info$exposure)
    methods::slot(scee, "exposureInfo") <- info
    expect_error(methods::validObject(scee), "must match")
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
    methods::slot(scee, "exposureData") <- infinite_data
    expect_error(methods::validObject(scee), "not NaN or infinite")

    scee <- .make_scee_invariant()
    extra_data <- rbind(
        exposureData(scee),
        D3 = c(E1 = 0, E2 = 0)
    )
    methods::slot(scee, "exposureData") <- extra_data
    expect_error(methods::validObject(scee), "not represented in sampleMap")
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
