test_that("sparse pseudobulk aggregation equals direct donor sums", {
    set.seed(20260712)
    counts <- matrix(
        rpois(8L * 18L, lambda = 3),
        nrow = 8L,
        dimnames = list(paste0("g", seq_len(8L)), paste0("c", seq_len(18L)))
    )
    donor <- rep(c("d3", "d1", "d2"), times = c(7L, 5L, 6L))
    celltype <- rep(c("A", "B", "A"), times = c(7L, 5L, 6L))

    expected <- cbind(
        d2 = rowSums(counts[, donor == "d2" & celltype == "A", drop = FALSE]),
        d3 = rowSums(counts[, donor == "d3" & celltype == "A", drop = FALSE])
    )

    dense_result <- exposomeSC:::.pseudobulk_aggregate(
        counts,
        donor,
        celltype,
        "A",
        min_cells = 6L
    )
    sparse_result <- exposomeSC:::.pseudobulk_aggregate(
        methods::as(Matrix::Matrix(counts), "CsparseMatrix"),
        donor,
        celltype,
        "A",
        min_cells = 6L
    )

    expect_identical(dense_result$valid_donors, c("d2", "d3"))
    expect_equal(dense_result$pb_mat, expected, tolerance = 0)
    expect_equal(sparse_result$pb_mat, expected, tolerance = 0)
    expect_identical(dense_result$n_cells, c(6L, 7L))
})

test_that("pseudobulk aggregation handles absent and single-cell strata", {
    counts <- matrix(
        seq_len(12L),
        nrow = 3L,
        dimnames = list(paste0("g", 1:3), paste0("c", 1:4))
    )
    donor <- c("d1", "d2", "d2", "d3")
    celltype <- c("A", "A", "B", "B")

    expect_null(exposomeSC:::.pseudobulk_aggregate(
        counts,
        donor,
        celltype,
        "missing",
        min_cells = 1L
    ))
    result <- exposomeSC:::.pseudobulk_aggregate(
        counts,
        donor,
        celltype,
        "A",
        min_cells = 1L
    )
    expect_equal(result$pb_mat[, "d1"], counts[, "c1"], tolerance = 0)
    expect_equal(result$pb_mat[, "d2"], counts[, "c2"], tolerance = 0)
})

test_that("aggregation matches rowsum on sparse counts with missing labels", {
    set.seed(20261030)
    n_cells <- 600L
    counts <- Matrix::rsparsematrix(50L, n_cells, density = 0.1,
        rand.x = function(n) as.numeric(rpois(n, 3) + 1))
    dimnames(counts) <- list(paste0("g", 1:50), paste0("c", seq_len(n_cells)))
    donor <- sample(c("D10", "D2", "D1", "D3"), n_cells, replace = TRUE)
    celltype <- sample(c("A", "B", NA), n_cells, replace = TRUE,
        prob = c(0.6, 0.38, 0.02))

    result <- exposomeSC:::.pseudobulk_aggregate(
        counts, donor, celltype, "A", min_cells = 1L
    )
    in_a <- !is.na(celltype) & celltype == "A"
    expected <- t(rowsum(t(as.matrix(counts[, in_a])), donor[in_a]))
    expected <- expected[, names(table(donor[in_a])), drop = FALSE]
    expect_identical(result$valid_donors, colnames(expected))
    expect_identical(dimnames(result$pb_mat), dimnames(expected))
    expect_equal(result$pb_mat, expected, tolerance = 0)
    expect_identical(result$n_cells,
        as.integer(table(donor[in_a])[result$valid_donors]))
})

test_that("donor_pseudobulk returns one SummarizedExperiment per cell type", {
    set.seed(1)
    counts <- matrix(rpois(20L * 120L, 4), nrow = 20L,
        dimnames = list(paste0("G", 1:20), paste0("c", 1:120)))
    donor <- rep(paste0("D", 1:6), each = 20L)
    cell_type <- rep(c("Mono", "NK", "Mono", "B"), 30L)
    cell_type[donor == "D6" & cell_type == "B"] <- "Mono"
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = counts),
        colData = S4Vectors::DataFrame(donor_id = donor,
            cell_type = cell_type))
    exposure <- matrix(c(seq_len(6), c(NA, 2:6)), ncol = 2,
        dimnames = list(paste0("D", 1:6), c("E1", "age")))
    scee <- build_scee(sce, exposure, sample_col = "donor_id")

    pb <- donor_pseudobulk(scee, min_cells = 5L)
    expect_s4_class(pb, "SimpleList")
    expect_identical(names(pb), c("B", "Mono", "NK"))
    expect_identical(S4Vectors::metadata(pb)$skipped, character())
    mono <- pb[["Mono"]]
    expect_s4_class(mono, "SummarizedExperiment")
    expect_identical(colnames(mono), paste0("D", 1:6))
    expect_identical(
        colnames(SummarizedExperiment::colData(mono)),
        c("donor_id", "cell_type", "n_cells", "E1", "age")
    )
    expect_equal(mono$age, c(NA, 2, 3, 4, 5, 6))
    expected <- vapply(paste0("D", 1:6), function(d) {
        rowSums(counts[, donor == d & cell_type == "Mono", drop = FALSE])
    }, numeric(20L))
    expect_equal(SummarizedExperiment::assay(mono, "counts"), expected,
        tolerance = 0)
    expect_identical(mono$n_cells,
        vapply(paste0("D", 1:6), function(d) {
            sum(donor == d & cell_type == "Mono")
        }, integer(1), USE.NAMES = FALSE))
    expect_identical(colnames(pb[["B"]]), paste0("D", 1:5))

    ## a cell type in which no donor reaches min_cells is listed as skipped
    pb_strict <- donor_pseudobulk(scee, min_cells = 11L)
    expect_identical(names(pb_strict), "Mono")
    expect_identical(S4Vectors::metadata(pb_strict)$skipped, c("B", "NK"))

    expect_error(donor_pseudobulk(scee, celltype_col = "nope"),
        "not a column")
    expect_error(donor_pseudobulk(scee, min_cells = 0L), "positive")
})

test_that("donor_pseudobulk gives the counts that run_sc_exwas fits", {
    set.seed(3)
    counts <- matrix(rnbinom(30L * 160L, mu = 20, size = 5), nrow = 30L,
        dimnames = list(paste0("G", 1:30), paste0("c", 1:160)))
    donor <- rep(paste0("D", 1:8), each = 20L)
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = counts),
        colData = S4Vectors::DataFrame(donor_id = donor, cell_type = "T"))
    exposure <- matrix(rep(c(0, 1), each = 4L), ncol = 1,
        dimnames = list(paste0("D", 1:8), "E1"))
    scee <- build_scee(sce, exposure, sample_col = "donor_id")

    fit <- run_sc_exwas(scee, exposure = "E1", min_group_donors = 4L,
        filter_min_donors = 4L)
    pb <- donor_pseudobulk(scee)[["T"]]
    diagnostics <- S4Vectors::metadata(fit)$diagnostics$T
    expect_identical(diagnostics$donor_ids, colnames(pb))
    expect_equal(diagnostics$full_library_size,
        colSums(SummarizedExperiment::assay(pb, "counts")))
})
