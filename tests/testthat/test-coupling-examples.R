# tests/testthat/test-coupling-examples.R

test_that("run_state_coupling returns one row per populated state bin", {
    skip_if_not_installed("metafor")
    set.seed(1)
    donors <- sprintf("D%02d", 1:12)
    donor <- rep(donors, each = 60)
    gene <- matrix(stats::rpois(2 * length(donor), 10), nrow = 2,
        dimnames = list(c("Gene1", "Gene2"), paste0("c", seq_along(donor))))
    protein <- matrix(stats::rpois(2 * length(donor), 5), nrow = 2,
        dimnames = list(c("Prot1", "Prot2"), colnames(gene)))
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = gene),
        colData = S4Vectors::DataFrame(donor = donor, celltype = "T",
            pseudotime = stats::runif(length(donor))))
    SingleCellExperiment::altExp(sce, "CITE") <-
        SummarizedExperiment::SummarizedExperiment(
            assays = list(counts = protein))
    scee <- build_scee(sce, matrix(seq(0, 2, length.out = 12), ncol = 1,
        dimnames = list(donors, "exposure")), sample_col = "donor")
    res <- run_state_coupling(scee, gene = "Gene1", protein = "Prot1",
        exposure = "exposure", celltype = "T", state_col = "pseudotime",
        n_bins = 3L, min_cells_per_bin = 10L, min_donors = 10L)
    expect_s3_class(res, "data.frame")
    expect_equal(nrow(res), 3L)
    slope <- attr(res, "slope_heterogeneity")
    expect_s3_class(slope, "data.frame")
    expect_true(is.finite(slope$pvalue))
    expect_error(run_state_coupling(scee, gene = "Gene1", protein = "Prot1",
        exposure = "nope", celltype = "T", state_col = "pseudotime"),
        "not in exposureData")
})

test_that("run_causal_mediation runs on a two-cell-type container", {
    set.seed(2)
    donors <- sprintf("D%02d", 1:12)
    n_mono <- 20:31
    donor <- c(rep(donors, times = n_mono), rep(donors, each = 25))
    cell_type <- c(rep("Mono", sum(n_mono)), rep("NK", 25 * 12))
    counts <- matrix(stats::rpois(30 * length(donor), 8), nrow = 30,
        dimnames = list(paste0("G", 1:30), paste0("c", seq_along(donor))))
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = counts),
        colData = S4Vectors::DataFrame(donor_id = donor,
            cell_type = cell_type))
    scee <- build_scee(sce, matrix(stats::rnorm(12), ncol = 1,
        dimnames = list(donors, "PM2.5")), sample_col = "donor_id")
    exposomeSC:::.reset_experimental_warnings()
    expect_warning(res <- run_causal_mediation(scee, exposure = "PM2.5",
        celltype = "Mono", genes = c("G1", "G2"), n_sims = 50L,
        sensitivity = TRUE), "experimental")
    expect_equal(nrow(res), 2L)
    expect_true(all(c("rho_at_zero", "method") %in% colnames(res)))
    expect_true(all(res$method %in% c("mediation", "difference_bootstrap")))
})

test_that("stability selection works when features outnumber donors", {
    set.seed(4)
    donors <- sprintf("D%02d", 1:20)
    donor <- rep(donors, each = 40)
    exposure <- stats::setNames(stats::rnorm(20), donors)
    counts <- matrix(stats::rpois(40 * length(donor), 20), nrow = 40,
        dimnames = list(paste0("G", 1:40), paste0("c", seq_along(donor))))
    effect <- c(0.8, 0.8, -0.8, -0.8)
    for (g in 1:4)
        counts[g, ] <- stats::rpois(length(donor),
            20 * exp(effect[g] * exposure[donor]))
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = counts),
        colData = S4Vectors::DataFrame(donor_id = donor, cell_type = "Mono"))
    scee <- build_scee(sce, matrix(exposure, ncol = 1,
        dimnames = list(donors, "PM2.5")), sample_col = "donor_id")
    metab <- matrix(stats::rnorm(60), nrow = 20,
        dimnames = list(donors, paste0("M", 1:3)))
    metab[, "M1"] <- metab[, "M1"] + 1.5 * exposure
    metab[, "M3"] <- 1
    net <- run_exposure_network(scee, metab, celltype = "Mono",
        exposure = "PM2.5", sample_col = "donor_id",
        selection_method = "stability", stability_q = 6L,
        network_method = "block_glasso")
    md <- net@metadata
    expect_equal(md$pfer_bound, 36 / (0.2 * 43))
    expect_true(md$transform %in% c("vst", "log_cpm"))
    expect_setequal(net@node_info$feature, c(paste0("G", 1:4), "M1"))
    expect_false("M3" %in% net@node_info$feature)
    expect_error(run_exposure_network(scee, metab, celltype = "Mono",
        exposure = "PM2.5", sample_col = "donor_id",
        selection_method = "stability", pi_threshold = 0.5), "pi_threshold")
    uni <- run_exposure_network(scee, metab, celltype = "Mono",
        exposure = "PM2.5", sample_col = "donor_id",
        selection_method = "univariate", fdr_threshold = 0.2,
        network_method = "block_glasso")
    expect_false("M3" %in% uni@node_info$feature)
})

test_that("run_exposure_network honours network_method", {
    set.seed(3)
    donors <- sprintf("D%02d", 1:20)
    donor <- rep(donors, each = 40)
    exposure <- stats::setNames(stats::rnorm(20), donors)
    counts <- matrix(stats::rpois(20 * length(donor), 20), nrow = 20,
        dimnames = list(paste0("G", 1:20), paste0("c", seq_along(donor))))
    effect <- c(0.8, 0.8, -0.8, -0.8)
    for (g in 1:4)
        counts[g, ] <- stats::rpois(length(donor),
            20 * exp(effect[g] * exposure[donor]))
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = counts),
        colData = S4Vectors::DataFrame(donor_id = donor, cell_type = "Mono"))
    scee <- build_scee(sce, matrix(exposure, ncol = 1,
        dimnames = list(donors, "PM2.5")), sample_col = "donor_id")
    metab <- matrix(stats::rnorm(60), nrow = 20,
        dimnames = list(donors, paste0("M", 1:3)))
    metab[, "M1"] <- metab[, "M1"] + 1.5 * exposure
    net <- run_exposure_network(scee, metab, celltype = "Mono",
        exposure = "PM2.5", sample_col = "donor_id",
        selection_method = "univariate", network_method = "block_glasso")
    expect_s4_class(net, "CelltypeNetworkResult")
    expect_setequal(net@node_info$feature, c(paste0("G", 1:4), "M1"))
    expect_identical(net@metadata$network_method, "block_glasso")
    expect_error(run_exposure_network(scee, metab, celltype = "Mono",
        exposure = "PM2.5", sample_col = "donor_id",
        selection_method = "univariate", network_method = "lasso"))
})
