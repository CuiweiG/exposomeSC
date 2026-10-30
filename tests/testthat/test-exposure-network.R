# tests/testthat/test-exposure-network.R

test_that("stability selection works when features outnumber donors", {
    skip_if_not_installed("glasso")
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
        exposure = "PM2.5",
        selection_method = "stability", stability_q = 6L,
        network_method = "block_glasso")
    md <- metadata(net)
    expect_equal(md$pfer_bound, 36 / (0.2 * 43))
    expect_true(md$transform %in% c("vst", "log_cpm"))
    expect_setequal(networkNodes(net)$feature, c(paste0("G", 1:4), "M1"))
    expect_false("M3" %in% networkNodes(net)$feature)
    expect_error(run_exposure_network(scee, metab, celltype = "Mono",
        exposure = "PM2.5",
        selection_method = "stability", pi_threshold = 0.5), "pi_threshold")
    uni <- run_exposure_network(scee, metab, celltype = "Mono",
        exposure = "PM2.5",
        selection_method = "univariate", fdr_threshold = 0.2,
        network_method = "block_glasso")
    expect_false("M3" %in% networkNodes(uni)$feature)
})

test_that("run_exposure_network honours network_method", {
    skip_if_not_installed("glasso")
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
        exposure = "PM2.5",
        selection_method = "univariate", network_method = "block_glasso")
    expect_s4_class(net, "CelltypeNetworkResult")
    expect_setequal(networkNodes(net)$feature, c(paste0("G", 1:4), "M1"))
    expect_identical(metadata(net)$network_method, "block_glasso")
    expect_error(run_exposure_network(scee, metab, celltype = "Mono",
        exposure = "PM2.5",
        selection_method = "univariate", network_method = "lasso"))
})
