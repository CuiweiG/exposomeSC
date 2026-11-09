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
            donor = cell_donor, cell_type = "T"))
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

    expect_s4_class(res, "DataFrame")
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
        colData = S4Vectors::DataFrame(donor = rep(donors, each = n_c),
            cell_type = "X"))
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
    expect_s4_class(res, "DataFrame")
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
            donor = cell_donor, cell_type = "T",
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

    expect_s4_class(res, "DataFrame")
    expect_true(nrow(res) >= 1)
    expect_equal(res$n_confounders, 3)  # 2 lib sizes + 1 cc_score
})

test_that("meta_moderators accepts a one-sided formula of column names", {
    skip_if_not_installed("metafor")
    skip_if_not_installed("ppcor")
    set.seed(1)
    donors <- sprintf("D%02d", 1:12)
    donor <- rep(donors, each = 60)
    gene <- matrix(stats::rpois(2 * length(donor), 10), nrow = 2,
        dimnames = list(c("Gene1", "Gene2"), paste0("c", seq_along(donor))))
    protein <- matrix(stats::rpois(2 * length(donor), 5), nrow = 2,
        dimnames = list(c("Prot1", "Prot2"), colnames(gene)))
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = gene),
        colData = S4Vectors::DataFrame(donor_id = donor, cell_type = "T"))
    SingleCellExperiment::altExp(sce, "CITE") <-
        SummarizedExperiment::SummarizedExperiment(
            assays = list(counts = protein))
    exp_mat <- cbind(exposure = seq(0, 2, length.out = 12),
        age = c(34, 51, 47, 29, 62, 55, 41, 38, 58, 45, 33, 49),
        bmi = c(22, 27, 31, 24, 26, 29, 23, 25, 30, 28, 21, 26))
    rownames(exp_mat) <- donors
    scee <- build_scee(sce, exp_mat, sample_col = "donor_id")
    couple <- function(moderators) {
        suppressMessages(run_cell_coupling(scee, genes = "Gene1",
            proteins = "Prot1", exposure = "exposure", celltype = "T",
            min_cells = 20L, min_donors = 10L, meta_moderators = moderators))
    }
    by_name <- couple(c("age", "bmi"))
    expect_equal(nrow(by_name), 1L)
    expect_equal(couple(~ age + bmi), by_name)
    for (moderators in list(~ age + I(age^2), ~ age:bmi, y ~ age, ~ age - 1,
                            ~ 1)) {
        expect_error(couple(moderators),
            "must be one-sided and only list columns of exposureData")
    }
})

test_that("run_cell_coupling leaves out donors without a moderator value", {
    skip_if_not_installed("metafor")
    skip_if_not_installed("ppcor")
    set.seed(1)
    donors <- sprintf("D%02d", 1:12)
    donor <- rep(donors, each = 60)
    gene <- matrix(stats::rpois(2 * length(donor), 10), nrow = 2,
        dimnames = list(c("Gene1", "Gene2"), paste0("c", seq_along(donor))))
    protein <- matrix(stats::rpois(2 * length(donor), 5), nrow = 2,
        dimnames = list(c("Prot1", "Prot2"), colnames(gene)))
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = gene),
        colData = S4Vectors::DataFrame(donor_id = donor, cell_type = "T"))
    SingleCellExperiment::altExp(sce, "CITE") <-
        SummarizedExperiment::SummarizedExperiment(
            assays = list(counts = protein))
    exp_mat <- cbind(exposure = seq(0, 2, length.out = 12),
        age = c(NA, 51, 47, 29, 62, 55, 41, 38, 58, 45, 33, 49))
    rownames(exp_mat) <- donors
    couple <- function(cells, min_donors = 10L) {
        scee <- build_scee(sce[, cells],
            exp_mat[unique(donor[cells]), , drop = FALSE],
            sample_col = "donor_id")
        run_cell_coupling(scee, genes = "Gene1", proteins = "Prot1",
            exposure = "exposure", celltype = "T", min_cells = 20L,
            min_donors = min_donors, meta_moderators = "age")
    }
    all_cells <- rep(TRUE, length(donor))
    expect_message(with_missing <- couple(all_cells),
        "Leaving 1 donor\\(s\\) with a missing exposure or moderator")
    expect_identical(with_missing$n_donors, 11L)
    without <- suppressMessages(couple(donor != "D01"))
    expect_equal(with_missing, without)
    ## eleven donors with an age are fewer than min_donors = 12
    expect_warning(too_few <- suppressMessages(couple(all_cells, 12L)),
        "Only 11 donors with >= 20 cells and complete exposure and moderator")
    expect_identical(nrow(too_few), 0L)
    ## donors without an exposure are counted the same way, also without
    ## moderators
    exp_mat[c("D01", "D02"), "exposure"] <- NA
    plain <- function(cells, min_donors) {
        scee <- build_scee(sce[, cells],
            exp_mat[unique(donor[cells]), , drop = FALSE],
            sample_col = "donor_id")
        run_cell_coupling(scee, genes = "Gene1", proteins = "Prot1",
            exposure = "exposure", celltype = "T", min_cells = 20L,
            min_donors = min_donors)
    }
    expect_warning(none <- suppressMessages(plain(all_cells, 11L)),
        "Only 10 donors with >= 20 cells and complete exposure")
    expect_identical(nrow(none), 0L)
    expect_message(ten <- plain(all_cells, 10L), "Leaving 2 donor\\(s\\)")
    expect_identical(ten$n_donors, 10L)
    expect_equal(ten,
        suppressMessages(plain(!donor %in% c("D01", "D02"), 10L)))
})

test_that("run_cell_coupling refuses an unordered confounder label", {
    skip_if_not_installed("metafor")
    skip_if_not_installed("ppcor")
    set.seed(1)
    donors <- sprintf("D%02d", 1:12)
    donor <- rep(donors, each = 60)
    gene <- matrix(stats::rpois(2 * length(donor), 10), nrow = 2,
        dimnames = list(c("Gene1", "Gene2"), paste0("c", seq_along(donor))))
    protein <- matrix(stats::rpois(2 * length(donor), 5), nrow = 2,
        dimnames = list(c("Prot1", "Prot2"), colnames(gene)))
    phase <- sample(c("G1", "S", "G2M"), length(donor), replace = TRUE)
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = gene),
        colData = S4Vectors::DataFrame(donor_id = donor, cell_type = "T",
            Phase = phase,
            Phase_ordered = factor(phase, levels = c("G1", "S", "G2M"),
                ordered = TRUE),
            Phase_factor = factor(phase),
            two_valued = ifelse(phase == "G1", "a", "b"),
            S_phase = as.integer(phase == "S"),
            G2M_phase = as.integer(phase == "G2M")))
    SingleCellExperiment::altExp(sce, "CITE") <-
        SummarizedExperiment::SummarizedExperiment(
            assays = list(counts = protein))
    exp_mat <- matrix(seq(0, 2, length.out = 12), ncol = 1,
        dimnames = list(donors, "exposure"))
    scee <- build_scee(sce, exp_mat, sample_col = "donor_id")
    couple <- function(confounders) {
        suppressMessages(run_cell_coupling(scee, genes = "Gene1",
            proteins = "Prot1", exposure = "exposure", celltype = "T",
            min_cells = 20L, min_donors = 10L, confounders = confounders))
    }
    expect_error(couple("Phase"),
        "more than two unordered categories: Phase")
    expect_error(couple("Phase_factor"),
        "more than two unordered categories: Phase_factor")
    ## a label with two values is used as before
    expect_identical(nrow(couple("two_valued")), 1L)
    expect_identical(nrow(couple("Phase_ordered")), 1L)
    expect_identical(nrow(couple(c("S_phase", "G2M_phase"))), 1L)
})
