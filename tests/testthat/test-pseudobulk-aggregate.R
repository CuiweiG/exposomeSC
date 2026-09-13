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
