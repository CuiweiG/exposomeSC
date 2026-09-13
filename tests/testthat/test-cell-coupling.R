# tests/testthat/test-cell-coupling.R

test_that("run_cell_coupling works on minimal data", {
    skip_if_not_installed("metafor")
    skip_if_not_installed("ppcor")

    ## Create minimal simulated data
    set.seed(42)
    n_d <- 15
    n_c <- 100
    n_g <- 2
    n_p <- 2
    donors <- paste0("D", seq_len(n_d))
    cell_donor <- rep(donors, each = n_c)
    exposure <- runif(n_d, 0, 3)

    ## Gene expression with exposure-dependent coupling
    gene_mat <- matrix(rpois(n_d * n_c * n_g, 10), nrow = n_g)
    prot_mat <- matrix(rpois(n_d * n_c * n_p, 5), nrow = n_p)

    cell_ids <- paste0("cell_", seq_len(n_d * n_c))
    colnames(gene_mat) <- cell_ids
    rownames(gene_mat) <- paste0("Gene", seq_len(n_g))
    colnames(prot_mat) <- cell_ids
    rownames(prot_mat) <- paste0("Prot", seq_len(n_p))

    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = Matrix::Matrix(gene_mat, sparse = TRUE)),
        colData = S4Vectors::DataFrame(
            donor = cell_donor, celltype = "T"))
    colnames(sce) <- cell_ids

    prot_se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = Matrix::Matrix(prot_mat, sparse = TRUE)))
    SingleCellExperiment::altExp(sce, "CITE") <- prot_se

    exp_mat <- matrix(exposure, ncol = 1,
                      dimnames = list(donors, "exposure"))
    scee <- build_scee(sce, exp_mat, sample_col = "donor")

    res <- run_cell_coupling(
        scee, genes = "Gene1", proteins = "Prot1",
        exposure = "exposure", celltype = "T",
        altexp_name = "CITE",
        min_cells = 10L, min_donors = 10L)

    expect_s3_class(res, "data.frame")
    expect_true(nrow(res) >= 1)
    expect_true("beta0" %in% names(res))
    expect_true("beta1" %in% names(res))
    expect_true("pval1" %in% names(res))
    expect_true("tau2" %in% names(res))
    expect_true("n_donors" %in% names(res))
    expect_true(res$n_donors >= 10)
})

test_that("run_cell_coupling validates inputs", {
    skip_if_not_installed("metafor")
    skip_if_not_installed("ppcor")

    set.seed(1)
    n_d <- 5; n_c <- 30
    donors <- paste0("D", 1:n_d)
    cell_ids <- paste0("c", seq_len(n_d * n_c))
    gene_mat <- Matrix::Matrix(rpois(2 * n_d * n_c, 10), nrow = 2, sparse = TRUE)
    colnames(gene_mat) <- cell_ids
    rownames(gene_mat) <- c("G1", "G2")
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = gene_mat),
        colData = S4Vectors::DataFrame(donor = rep(donors, each = n_c), celltype = "X"))
    colnames(sce) <- cell_ids
    prot_mat <- Matrix::Matrix(rpois(2 * n_d * n_c, 5), nrow = 2, sparse = TRUE)
    colnames(prot_mat) <- cell_ids
    rownames(prot_mat) <- c("P1", "P2")
    prot_se <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = prot_mat))
    SingleCellExperiment::altExp(sce, "CITE") <- prot_se
    exp_mat <- matrix(runif(n_d), ncol = 1, dimnames = list(donors, "exp"))
    scee <- build_scee(sce, exp_mat, sample_col = "donor")

    ## Too few donors -> warning + empty df
    expect_warning(
        res <- run_cell_coupling(scee, genes = "G1", proteins = "P1",
            exposure = "exp", celltype = "X", min_donors = 10L),
        "donors")
    expect_equal(nrow(res), 0)

    ## Wrong exposure name -> error
    expect_error(
        run_cell_coupling(scee, genes = "G1", proteins = "P1",
            exposure = "WRONG", celltype = "X"),
        "WRONG")

    ## Wrong celltype -> error
    expect_error(
        run_cell_coupling(scee, genes = "G1", proteins = "P1",
            exposure = "exp", celltype = "NONEXIST"),
        "NONEXIST")
})

test_that("run_cell_coupling with confounders", {
    skip_if_not_installed("metafor")
    skip_if_not_installed("ppcor")

    set.seed(99)
    n_d <- 15; n_c <- 100
    donors <- paste0("D", seq_len(n_d))
    cell_ids <- paste0("c", seq_len(n_d * n_c))
    cell_donor <- rep(donors, each = n_c)

    gene_mat <- Matrix::Matrix(rpois(2 * n_d * n_c, 10), nrow = 2, sparse = TRUE)
    colnames(gene_mat) <- cell_ids; rownames(gene_mat) <- c("G1", "G2")
    prot_mat <- Matrix::Matrix(rpois(2 * n_d * n_c, 5), nrow = 2, sparse = TRUE)
    colnames(prot_mat) <- cell_ids; rownames(prot_mat) <- c("P1", "P2")

    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = gene_mat),
        colData = S4Vectors::DataFrame(
            donor = cell_donor, celltype = "T",
            cc_score = runif(n_d * n_c)))
    colnames(sce) <- cell_ids
    prot_se <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = prot_mat))
    SingleCellExperiment::altExp(sce, "CITE") <- prot_se

    exp_mat <- matrix(runif(n_d, 0, 3), ncol = 1, dimnames = list(donors, "exp"))
    scee <- build_scee(sce, exp_mat, sample_col = "donor")

    res <- run_cell_coupling(scee, genes = "G1", proteins = "P1",
        exposure = "exp", celltype = "T",
        confounders = "cc_score",
        min_cells = 10L, min_donors = 10L)

    expect_s3_class(res, "data.frame")
    expect_true(nrow(res) >= 1)
    expect_equal(res$n_confounders, 3)  # 2 lib sizes + 1 cc_score
})
