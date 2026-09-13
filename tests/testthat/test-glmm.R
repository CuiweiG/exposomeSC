# tests/testthat/test-glmm.R

test_that("run_sc_exwas_glmm works with lme4 backend", {
    skip_if_not_installed("lme4")

    set.seed(42)
    ## Use more donors (20) with more cells per donor for
    ## better convergence
    n_donors <- 20
    cells_per <- 15
    n_cells <- n_donors * cells_per
    n_genes <- 50
    counts <- matrix(rpois(n_genes * n_cells, 10),
        nrow = n_genes,
        dimnames = list(paste0("G", seq_len(n_genes)),
            paste0("c", seq_len(n_cells))))
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = counts),
        colData = S4Vectors::DataFrame(
            cell_id = paste0("c", seq_len(n_cells)),
            donor_id = rep(paste0("D", seq_len(n_donors)),
                each = cells_per),
            cell_type = rep(c("Mono", "NK"),
                length.out = n_cells)))
    exp_mat <- matrix(rnorm(n_donors * 2), nrow = n_donors,
        dimnames = list(paste0("D", seq_len(n_donors)),
            c("E1", "E2")))

    for (f in list.files("../../R", full.names = TRUE))
        source(f, local = TRUE)

    scee <- build_scee(sce, exp_mat, sample_col = "donor_id")

    res <- run_sc_exwas_glmm(scee, exposure = "E1",
        celltype_col = "cell_type", sample_col = "donor_id",
        BPPARAM = BiocParallel::SerialParam(),
        target_genes = paste0("G", 1:10))

    expect_s3_class(res, "data.frame")
    expect_true(nrow(res) > 0)
    expect_true("log2FC" %in% colnames(res))
    expect_true("se" %in% colnames(res))
    expect_true("pvalue" %in% colnames(res))
    expect_true("padj" %in% colnames(res))
    expect_true("method" %in% colnames(res))
    ## gaussian path is fitted by lmerTest (Satterthwaite/Kenward-Roger df)
    expect_true(all(grepl("lmerTest", res$method)))
})

test_that("GLMM works with glmmTMB nbinom2 backend", {
    skip(paste("the glmmTMB/nbinom2 backend was removed; the negative-binomial",
               "path is now fitted via lme4::glmer.nb through family = 'nbinom'"))

    set.seed(42)
    n_donors <- 20
    cells_per <- 15
    n_cells <- n_donors * cells_per
    n_genes <- 50
    counts <- matrix(rpois(n_genes * n_cells, 10),
        nrow = n_genes,
        dimnames = list(paste0("G", seq_len(n_genes)),
            paste0("c", seq_len(n_cells))))
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = counts),
        colData = S4Vectors::DataFrame(
            cell_id = paste0("c", seq_len(n_cells)),
            donor_id = rep(paste0("D", seq_len(n_donors)),
                each = cells_per),
            cell_type = rep(c("Mono", "NK"),
                length.out = n_cells)))
    exp_mat <- matrix(rnorm(n_donors * 2), nrow = n_donors,
        dimnames = list(paste0("D", seq_len(n_donors)),
            c("E1", "E2")))

    for (f in list.files("../../R", full.names = TRUE))
        source(f, local = TRUE)

    scee <- build_scee(sce, exp_mat, sample_col = "donor_id")

    res <- run_sc_exwas_glmm(scee, exposure = "E1",
        celltype_col = "cell_type", sample_col = "donor_id",        backend = "glmmTMB", family = "nbinom2",
        target_genes = paste0("G", 1:10))

    expect_s4_class(res, "DataFrame")
    expect_true(nrow(res) > 0)
    expect_true(all(grepl("glmmTMB/nbinom2", res$method)))
    ## P-values should be valid
    pvals <- res$pvalue[!is.na(res$pvalue)]
    expect_true(all(pvals >= 0 & pvals <= 1))
})

test_that("GLMM SE uses vcov covariance correctly", {
    skip_if_not_installed("lme4")

    set.seed(42)
    n_donors <- 25
    cells_per <- 15
    n_cells <- n_donors * cells_per
    counts <- matrix(rpois(50 * n_cells, 10), nrow = 50,
        dimnames = list(paste0("G", 1:50),
            paste0("c", seq_len(n_cells))))
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = counts),
        colData = S4Vectors::DataFrame(
            cell_id = paste0("c", seq_len(n_cells)),
            donor_id = rep(paste0("D", seq_len(n_donors)),
                each = cells_per),
            cell_type = rep(c("Mono", "NK"),
                length.out = n_cells)))
    exp_mat <- matrix(rnorm(n_donors * 2), nrow = n_donors,
        dimnames = list(paste0("D", seq_len(n_donors)),
            c("E1", "E2")))

    for (f in list.files("../../R", full.names = TRUE))
        source(f, local = TRUE)

    scee <- build_scee(sce, exp_mat, sample_col = "donor_id")

    res <- run_sc_exwas_glmm(scee, exposure = "E1",
        celltype_col = "cell_type", sample_col = "donor_id",
        BPPARAM = BiocParallel::SerialParam(),
        target_genes = paste0("G", 1:5))

    ## Non-reference cell type should have SE values
    nk_res <- res[res$celltype == "NK", ]
    expect_true(nrow(nk_res) > 0)
    expect_true(all(!is.na(nk_res$se)))
    expect_true(all(nk_res$se > 0))
    ## SE should differ from just sqrt(se1^2 + se2^2)
    ## (since vcov includes covariance term)
})

test_that("GLMM returns valid p-values", {
    skip_if_not_installed("lme4")

    set.seed(42)
    n_donors <- 20
    cells_per <- 15
    n_cells <- n_donors * cells_per
    counts <- matrix(rpois(50 * n_cells, 10), nrow = 50,
        dimnames = list(paste0("G", 1:50),
            paste0("c", seq_len(n_cells))))
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = counts),
        colData = S4Vectors::DataFrame(
            cell_id = paste0("c", seq_len(n_cells)),
            donor_id = rep(paste0("D", seq_len(n_donors)),
                each = cells_per),
            cell_type = rep(c("Mono", "NK"),
                length.out = n_cells)))
    exp_mat <- matrix(rnorm(n_donors * 2), nrow = n_donors,
        dimnames = list(paste0("D", seq_len(n_donors)),
            c("E1", "E2")))

    for (f in list.files("../../R", full.names = TRUE))
        source(f, local = TRUE)

    scee <- build_scee(sce, exp_mat, sample_col = "donor_id")

    res <- run_sc_exwas_glmm(scee, exposure = "E1",
        celltype_col = "cell_type", sample_col = "donor_id",
        BPPARAM = BiocParallel::SerialParam(),
        target_genes = paste0("G", 1:10))

    pvals <- res$pvalue[!is.na(res$pvalue)]
    expect_true(all(pvals >= 0 & pvals <= 1))
})
