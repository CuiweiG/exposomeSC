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
        network_method = "block_glasso", transform = "log_cpm")
    md <- metadata(net)
    ## the bound counts the 42 candidates that vary across the donors, not
    ## the constant metabolite M3, which can never be selected
    expect_identical(md$n_total_features, 43L)
    expect_identical(md$constant_features, "M3")
    expect_equal(md$pfer_bound, 36 / (0.2 * 42))
    expect_identical(md$transform, "log_cpm")
    expect_setequal(networkNodes(net)$feature, c(paste0("G", 1:4), "M1"))
    expect_false("M3" %in% networkNodes(net)$feature)
    ## so does the default stability_q: with a gene without counts as well,
    ## 41 candidates vary, and floor(sqrt(0.6 * 41)) is 4, where
    ## floor(sqrt(0.6 * 43)) would be 5
    counts_gap <- counts
    counts_gap["G40", ] <- 0L
    sce_gap <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = counts_gap),
        colData = S4Vectors::DataFrame(donor_id = donor, cell_type = "Mono"))
    scee_gap <- build_scee(sce_gap, matrix(exposure, ncol = 1,
        dimnames = list(donors, "PM2.5")), sample_col = "donor_id")
    by_default <- suppressMessages(run_exposure_network(scee_gap, metab,
        celltype = "Mono", exposure = "PM2.5",
        selection_method = "stability", pi_threshold = 0.8,
        network_method = "block_glasso", transform = "log_cpm"))
    expect_identical(metadata(by_default)$stability_q, 4L)
    expect_equal(metadata(by_default)$pfer_bound, 16 / (0.6 * 41))
    expect_setequal(metadata(by_default)$constant_features, c("G40", "M3"))
    expect_error(run_exposure_network(scee, metab, celltype = "Mono",
        exposure = "PM2.5", selection_method = "stability",
        stability_q = 43L, network_method = "block_glasso",
        transform = "log_cpm"),
        "the number of features that vary across the donors (42)",
        fixed = TRUE)
    ## features adjusted for the exposure cannot be screened against it
    expect_error(run_exposure_network(scee, metab, celltype = "Mono",
        exposure = "PM2.5", covariates = "PM2.5",
        selection_method = "univariate", network_method = "block_glasso",
        transform = "log_cpm"),
        "The exposure must not also be listed as a covariate.", fixed = TRUE)
    expect_error(run_exposure_network(scee, metab, celltype = "Mono",
        exposure = "PM2.5",
        selection_method = "stability", pi_threshold = 0.5,
        transform = "log_cpm"), "pi_threshold")
    uni <- run_exposure_network(scee, metab, celltype = "Mono",
        exposure = "PM2.5",
        selection_method = "univariate", fdr_threshold = 0.2,
        network_method = "block_glasso", transform = "log_cpm")
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
        selection_method = "univariate", network_method = "block_glasso",
        transform = "log_cpm")
    expect_s4_class(net, "CelltypeNetworkResult")
    expect_setequal(networkNodes(net)$feature, c(paste0("G", 1:4), "M1"))
    expect_identical(metadata(net)$network_method, "block_glasso")
    ## the object is labelled as a network restricted to the selected
    ## features; the estimator is in the metadata
    expect_identical(networkMethod(net), "exposure_restricted")
    expect_true("  Method: exposure_restricted" %in%
        utils::capture.output(show(net)))
    expect_error(run_exposure_network(scee, metab, celltype = "Mono",
        exposure = "PM2.5",
        selection_method = "univariate", network_method = "lasso",
        transform = "log_cpm"))
})

test_that("the switch from coglasso to block glasso is announced", {
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
    ## no metabolite is associated with the exposure, so only transcripts
    ## pass the univariate screen
    metab <- matrix(stats::rnorm(60), nrow = 20,
        dimnames = list(donors, paste0("M", 1:3)))
    expect_warning(
        net <- suppressMessages(run_exposure_network(scee, metab,
            celltype = "Mono", exposure = "PM2.5",
            selection_method = "univariate", network_method = "coglasso",
            transform = "log_cpm")),
        "one omic layer"
    )
    expect_identical(metadata(net)$network_method, "block_glasso")
    expect_true(all(networkNodes(net)$omic_layer == "transcript"))
})

test_that("run_exposure_network takes and records the glasso penalties", {
    skip_if_not_installed("glasso")
    set.seed(4)
    donors <- sprintf("D%02d", 1:20)
    donor <- rep(donors, each = 40)
    exposure <- stats::setNames(stats::rnorm(20), donors)
    counts <- matrix(stats::rpois(40 * length(donor), 20), nrow = 40,
        dimnames = list(paste0("G", 1:40), paste0("c", seq_along(donor))))
    for (g in 1:4)
        counts[g, ] <- stats::rpois(length(donor),
            20 * exp(c(0.8, 0.8, -0.8, -0.8)[g] * exposure[donor]))
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = counts),
        colData = S4Vectors::DataFrame(donor_id = donor, cell_type = "Mono"))
    scee <- build_scee(sce, matrix(exposure, ncol = 1,
        dimnames = list(donors, "PM2.5")), sample_col = "donor_id")
    metab <- matrix(stats::rnorm(60), nrow = 20,
        dimnames = list(donors, paste0("M", 1:3)))
    metab[, "M1"] <- metab[, "M1"] + 1.5 * exposure
    network <- function(...) {
        run_exposure_network(scee, metab, celltype = "Mono",
            exposure = "PM2.5", selection_method = "univariate",
            fdr_threshold = 0.2, network_method = "block_glasso", ...,
            transform = "log_cpm")
    }
    default <- network()
    expect_identical(metadata(default)$selection,
        list(lambda_w = 0.3, lambda_b = 0.5))
    sparse <- network(lambda_w = 0.9, lambda_b = 0.9)
    expect_identical(metadata(sparse)$selection,
        list(lambda_w = 0.9, lambda_b = 0.9))
    expect_lt(sum(networkAdjacency(sparse)), sum(networkAdjacency(default)))
    ## coglasso records the values it selected, c among them
    skip_if_not_installed("coglasso", minimum_version = "1.1.0")
    selected <- suppressMessages(run_exposure_network(scee, metab,
        celltype = "Mono", exposure = "PM2.5",
        selection_method = "univariate", fdr_threshold = 0.2,
        transform = "log_cpm"))
    expect_identical(metadata(selected)$network_method, "coglasso")
    expect_named(metadata(selected)$selection, c("nlambda_w", "nlambda_b",
        "selected_lambda_w", "selected_lambda_b", "selected_c",
        "stars_variability"))
    expect_identical(metadata(selected)$selection[1:2],
        list(nlambda_w = 10L, nlambda_b = 10L))
    expect_true(all(vapply(metadata(selected)$selection, function(value)
        is.numeric(value) && length(value) == 1L && is.finite(value),
        logical(1))))
})

test_that("the rank correlation without covariates is that of cor.test()", {
    values <- cbind(
        a = c(3, 1, 4, 1, 5, 9, 2, 6, 5, 3, 5, 8),
        b = c(2, 7, 1, 8, 2, 8, 1, 8, 2, 8, 4, 5),
        c = c(0, 0, 0, 2, 0, 0, 0, 0, 7, 0, 0, 1))
    response <- c(1.2, 0.4, 2.2, 0.9, 1.8, 3.1, 0.2, 2.5, 1.1, 0.7, 2.9, 1.6)
    result <- .partial_rank_correlation(values, response)
    expect_equal(result$rho,
        unname(stats::cor(values, response, method = "spearman")[, 1]),
        tolerance = 1e-12)
    by_test <- apply(values, 2L, function(v) {
        stats::cor.test(v, response, method = "spearman",
            exact = FALSE)$p.value
    })
    expect_equal(result$p_value, unname(by_test), tolerance = 1e-12)
})

test_that("the partial rank correlation is that of the ranks after lm()", {
    values <- cbind(
        a = c(3, 1, 4, 1, 5, 9, 2, 6, 5, 3, 5, 8),
        b = c(2, 7, 1, 8, 2, 8, 1, 8, 2, 8, 4, 5),
        c = c(0, 0, 0, 2, 0, 0, 0, 0, 7, 0, 0, 1))
    response <- c(1.2, 0.4, 2.2, 0.9, 1.8, 3.1, 0.2, 2.5, 1.1, 0.7, 2.9, 1.6)
    c1 <- c(5, 3, 8, 1, 9, 12, 2, 10, 6, 4, 11, 7)
    c2 <- c(1, 1, 2, 2, 1, 2, 1, 2, 1, 1, 2, 2)
    n <- length(response)
    one <- function(v) stats::residuals(stats::lm(rank(v) ~ rank(c1)))
    two <- function(v) {
        stats::residuals(stats::lm(rank(v) ~ rank(c1) + rank(c2)))
    }
    reference <- function(residuals_of, k) {
        rho <- apply(values, 2L, function(v) {
            stats::cor(residuals_of(v), residuals_of(response))
        })
        statistic <- rho * sqrt((n - 2 - k) / (1 - rho^2))
        list(rho = unname(rho),
            p_value = unname(2 * stats::pt(-abs(statistic), n - 2 - k)))
    }
    expect_equal(.partial_rank_correlation(values, response, cbind(c1)),
        reference(one, 1), tolerance = 1e-12)
    expect_equal(.partial_rank_correlation(values, response, cbind(c1, c2)),
        reference(two, 2), tolerance = 1e-12)
})

test_that("tied values are not ranked by the covariates", {
    ## a feature with a value in one donor only, and a covariate that
    ## ranks the donors almost as the response does
    covariate <- 1:12
    response <- c(1, 3, 2, 4, 6, 5, 7, 9, 8, 10, 12, 11)
    sparse <- c(0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0)
    ## the residuals of a linear adjustment rank the eleven zeros by the
    ## covariate: rho is -0.825
    linear <- stats::cor(stats::residuals(stats::lm(sparse ~ covariate)),
        response, method = "spearman")
    expect_lt(linear, -0.8)
    ## the partial rank correlation keeps them tied: rho is 0.0052
    partial <- .partial_rank_correlation(cbind(sparse), response,
        cbind(covariate))
    expect_lt(abs(partial$rho), 0.01)
    expect_gt(partial$p_value, 0.9)
})

test_that("covariates that repeat each other take one degree of freedom", {
    covariate <- c(5, 3, 8, 1, 9, 12, 2, 10, 6, 4, 11, 7)
    response <- c(1.2, 0.4, 2.2, 0.9, 1.8, 3.1, 0.2, 2.5, 1.1, 0.7, 2.9, 1.6)
    values <- cbind(a = c(3, 1, 4, 1, 5, 9, 2, 6, 5, 3, 5, 8),
        b = c(2, 7, 1, 8, 2, 8, 1, 8, 2, 8, 4, 5))
    once <- .partial_rank_correlation(values, response, cbind(covariate))
    ## the same covariate twice and a constant one
    twice <- .partial_rank_correlation(values, response,
        cbind(covariate, covariate, 7))
    expect_equal(twice, once)
})

test_that("what the adjustment reproduces has no rank correlation", {
    covariate <- c(5, 3, 8, 1, 9, 12, 2, 10, 6, 4, 11, 7)
    response <- c(1.2, 0.4, 2.2, 0.9, 1.8, 3.1, 0.2, 2.5, 1.1, 0.7, 2.9, 1.6)
    values <- cbind(copy = covariate, constant = 5,
        varying = c(3, 1, 4, 1, 5, 9, 2, 6, 5, 3, 5, 8))
    ## a copy of the covariate and a constant column leave rounding noise
    ## as residuals, which must not be correlated with anything
    result <- .partial_rank_correlation(values, response, cbind(covariate))
    expect_identical(is.na(result$rho), c(TRUE, TRUE, FALSE))
    expect_identical(is.na(result$p_value), c(TRUE, TRUE, FALSE))
    ## a constant column, without an adjustment
    unadjusted <- .partial_rank_correlation(values, response)
    expect_identical(is.na(unadjusted$rho), c(FALSE, TRUE, FALSE))
    expect_identical(is.na(unadjusted$p_value), c(FALSE, TRUE, FALSE))
    ## a response that the adjustment reproduces
    reproduced <- .partial_rank_correlation(values, 2 * covariate + 1,
        cbind(covariate))
    expect_true(all(is.na(reproduced$rho)))
    expect_true(all(is.na(reproduced$p_value)))
})

## Twenty donors with forty cells each and forty genes, of which G1 to G4
## follow the exposure, and three metabolites: M1 follows the exposure and
## M3 is constant. The exposure of D01 is not known, and G39 has counts in
## the cells of D01 only.
.missing_exposure_fixture <- function() {
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
    counts["G39", donor != "D01"] <- 0L
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = counts),
        colData = S4Vectors::DataFrame(donor_id = donor, cell_type = "Mono"))
    recorded <- exposure
    recorded["D01"] <- NA
    scee <- build_scee(sce, matrix(recorded, ncol = 1,
        dimnames = list(donors, "PM2.5")), sample_col = "donor_id")
    metab <- matrix(stats::rnorm(60), nrow = 20,
        dimnames = list(donors, paste0("M", 1:3)))
    metab[, "M1"] <- metab[, "M1"] + 1.5 * exposure
    metab[, "M3"] <- 1
    list(scee = scee, metab = metab, donors = donors)
}

test_that("the selection is made among the donors with an exposure", {
    skip_if_not_installed("glasso")
    fx <- .missing_exposure_fixture()
    stable <- suppressMessages(run_exposure_network(fx$scee, fx$metab,
        celltype = "Mono", exposure = "PM2.5",
        selection_method = "stability", pi_threshold = 0.8,
        network_method = "block_glasso", transform = "log_cpm"))
    md <- metadata(stable)
    ## G39, with counts in D01 only, is constant across the other donors
    expect_setequal(md$constant_features, c("G39", "M3"))
    ## 41 of the 43 candidates vary: floor(sqrt(0.6 * 41)) is 4, where
    ## floor(sqrt(0.6 * 42)) would be 5
    expect_identical(md$stability_q, 4L)
    expect_equal(md$pfer_bound, 16 / (0.6 * 41))
    expect_identical(md$n_donors, 20L)
    expect_identical(md$n_donors_exposure, 19L)
    expect_identical(md$donors, fx$donors)
    ## G39 is not tested, so its test gives no warning
    expect_no_warning(univariate <- suppressMessages(run_exposure_network(
        fx$scee, fx$metab, celltype = "Mono", exposure = "PM2.5",
        selection_method = "univariate", fdr_threshold = 0.2,
        network_method = "block_glasso", transform = "log_cpm")))
    expect_setequal(metadata(univariate)$constant_features, c("G39", "M3"))
    expect_identical(metadata(univariate)$n_donors_exposure, 19L)
    expect_identical(metadata(univariate)$donors, fx$donors)
})

test_that("no feature that varies with an exposure is an error", {
    donors <- sprintf("D%02d", 1:12)
    donor <- rep(donors, each = 12)
    ## the same counts in every cell, so every gene is constant
    counts <- matrix(10L, nrow = 4, ncol = length(donor),
        dimnames = list(paste0("G", 1:4), paste0("c", seq_along(donor))))
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = counts),
        colData = S4Vectors::DataFrame(donor_id = donor, cell_type = "Mono"))
    scee <- build_scee(sce, matrix(as.numeric(1:12), ncol = 1,
        dimnames = list(donors, "PM2.5")), sample_col = "donor_id")
    metab <- matrix(1, nrow = 12, ncol = 2,
        dimnames = list(donors, c("M1", "M2")))
    for (method in c("univariate", "stability")) {
        expect_error(
            suppressMessages(run_exposure_network(scee, metab,
                celltype = "Mono", exposure = "PM2.5",
                selection_method = method, network_method = "block_glasso",
                transform = "log_cpm")),
            paste("No candidate feature varies across the donors with a",
                "value of the exposure."), fixed = TRUE)
    }
})

## Twenty donors with twenty cells each. The exposure follows age, a
## covariate, and a deviation from it, which M1 to M3 follow; M4 follows
## neither and M5 follows age. The genes S1 to S5 have counts only in D09
## and D10, and G1 to G20 have none of these associations.
.covariate_fixture <- function() {
    set.seed(5)
    donors <- sprintf("D%02d", 1:20)
    donor <- rep(donors, each = 20)
    i <- 1:20
    age <- stats::setNames(as.numeric(i), donors)
    deviation <- stats::setNames(sin(i), donors)
    exposure <- age + 6 * deviation
    counts <- matrix(stats::rpois(25 * length(donor), 20), nrow = 25,
        dimnames = list(c(paste0("G", 1:20), paste0("S", 1:5)),
            paste0("c", seq_along(donor))))
    counts[21:25, !donor %in% c("D09", "D10")] <- 0L
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = counts),
        colData = S4Vectors::DataFrame(donor_id = donor, cell_type = "Mono"))
    scee <- build_scee(sce, cbind(PM2.5 = exposure, age = age),
        sample_col = "donor_id")
    metab <- cbind(M1 = deviation + 0.1 * cos(i),
        M2 = -deviation + 0.1 * sin(2 * i), M3 = deviation + 0.1 * cos(3 * i),
        M4 = cos(5 * i), M5 = age + 0.5 * cos(7 * i))
    rownames(metab) <- donors
    list(scee = scee, metab = metab, exposure = exposure, age = age)
}

test_that("sparse genes are not selected for a covariate of the exposure", {
    skip_if_not_installed("glasso")
    fx <- .covariate_fixture()
    expect_gt(stats::cor(fx$exposure, fx$age, method = "spearman"), 0.75)
    values <- .assemble_crossomic(fx$scee, fx$metab, "Mono",
        transform = "log_cpm", drop_constant = FALSE,
        keep_unadjusted = TRUE)$unadjusted
    sparse <- paste0("S", 1:5)
    ## ranked after a linear adjustment for age, the zeros of S1 to S5 are
    ## ordered by age, and the genes would be selected
    linear <- apply(values, 2L, function(v) {
        adjusted <- stats::residuals(stats::lm(v ~ fx$age))
        stats::cor.test(adjusted, fx$exposure, method = "spearman",
            exact = FALSE)$p.value
    })
    expect_true(all(sparse %in%
        colnames(values)[stats::p.adjust(linear, "BH") < 0.05]))
    adjusted <- suppressMessages(run_exposure_network(fx$scee, fx$metab,
        celltype = "Mono", exposure = "PM2.5", covariates = "age",
        selection_method = "univariate", network_method = "block_glasso",
        transform = "log_cpm"))
    expect_false(any(sparse %in% networkNodes(adjusted)$feature))
    expect_identical(networkNodes(adjusted)$feature, c("M1", "M2", "M3"))
})

test_that("without covariates the selection is that of cor.test()", {
    skip_if_not_installed("glasso")
    fx <- .covariate_fixture()
    values <- .assemble_crossomic(fx$scee, fx$metab, "Mono",
        transform = "log_cpm", drop_constant = FALSE,
        keep_unadjusted = TRUE)$unadjusted
    p_values <- apply(values, 2L, function(v) {
        stats::cor.test(v, fx$exposure, method = "spearman",
            exact = FALSE)$p.value
    })
    by_test <- colnames(values)[stats::p.adjust(p_values, "BH") < 0.5]
    selected <- suppressMessages(run_exposure_network(fx$scee, fx$metab,
        celltype = "Mono", exposure = "PM2.5", selection_method = "univariate",
        fdr_threshold = 0.5, network_method = "block_glasso",
        transform = "log_cpm"))
    expect_gte(length(by_test), 3L)
    expect_identical(networkNodes(selected)$feature, by_test)
})

test_that("a tied exposure is tested as by cor.test()", {
    values <- cbind(a = c(3, 1, 4, 1, 5, 9, 2, 6, 5, 3, 5, 8),
        c = c(0, 0, 0, 2, 0, 0, 0, 0, 7, 0, 0, 1))
    exposed <- c(0, 0, 1, 0, 1, 1, 0, 0, 1, 0, 1, 0)
    by_test <- apply(values, 2L, function(v) {
        stats::cor.test(v, exposed, method = "spearman",
            exact = FALSE)$p.value
    })
    expect_equal(.partial_rank_correlation(values, exposed)$p_value,
        unname(by_test), tolerance = 1e-12)
})

test_that("the selector of stability selection ranks given the covariates", {
    covariate <- 1:12
    deviation <- rep(c(0, 1, -1), 4)
    response <- covariate + deviation
    values <- cbind(follows_covariate = covariate,
        follows_deviation = deviation,
        neither = c(3, 1, 4, 1, 5, 9, 2, 6, 5, 3, 5, 8))
    ## unadjusted, the copy of the covariate is ranked first
    expect_identical(.top_marginal_features(values, response, 1L),
        c(1, 0, 0))
    ## given the covariate, the copy has nothing left, and the deviation,
    ## which is what the response adds to the covariate, is ranked first
    expect_identical(.top_marginal_features(values, response, 1L,
        cbind(covariate)), c(0, 1, 0))
})

test_that("stability selection is made given the covariates", {
    skip_if_not_installed("glasso")
    fx <- .covariate_fixture()
    values <- .assemble_crossomic(fx$scee, fx$metab, "Mono",
        transform = "log_cpm", drop_constant = FALSE,
        keep_unadjusted = TRUE)$unadjusted
    donors <- rownames(values)
    ## the half-samples of the function, drawn here with the same seed
    set.seed(11)
    frequency <- numeric(ncol(values))
    for (b in seq_len(100L)) {
        idx <- sample.int(20L, 10L)
        frequency <- frequency + .top_marginal_features(
            values[idx, , drop = FALSE], fx$exposure[donors][idx], 3L,
            cbind(fx$age[donors][idx]))
    }
    expected <- colnames(values)[frequency >= 60]
    expect_gte(length(expected), 3L)
    expect_false(any(paste0("S", 1:5) %in% expected))
    set.seed(11)
    selected <- suppressMessages(run_exposure_network(fx$scee, fx$metab,
        celltype = "Mono", exposure = "PM2.5", covariates = "age",
        selection_method = "stability", stability_q = 3L,
        network_method = "block_glasso", transform = "log_cpm"))
    expect_identical(networkNodes(selected)$feature, expected)
})

test_that("donors without an exposure are left out with their covariates", {
    skip_if_not_installed("glasso")
    fx <- .covariate_fixture()
    exposures <- exposureData(fx$scee)
    exposures["D20", "PM2.5"] <- NA
    scee <- build_scee(methods::as(fx$scee, "SingleCellExperiment"),
        exposures, sample_col = "donor_id")
    values <- .assemble_crossomic(scee, fx$metab, "Mono",
        transform = "log_cpm", drop_constant = FALSE,
        keep_unadjusted = TRUE)$unadjusted
    with_value <- setdiff(rownames(values), "D20")
    test <- .partial_rank_correlation(values[with_value, ],
        exposures[with_value, "PM2.5"],
        exposures[with_value, "age", drop = FALSE])
    by_hand <- colnames(values)[
        which(stats::p.adjust(test$p_value, "BH") < 0.5)]
    selected <- suppressMessages(run_exposure_network(scee, fx$metab,
        celltype = "Mono", exposure = "PM2.5", covariates = "age",
        selection_method = "univariate", fdr_threshold = 0.5,
        network_method = "block_glasso", transform = "log_cpm"))
    expect_gte(length(by_hand), 3L)
    expect_identical(networkNodes(selected)$feature, by_hand)
    expect_identical(metadata(selected)$n_donors_exposure, 19L)
})

test_that("too few donors with an exposure are an error, with covariates", {
    fx <- .covariate_fixture()
    with_exposure <- function(n) {
        exposures <- exposureData(fx$scee)
        exposures[-seq_len(n), "PM2.5"] <- NA
        build_scee(methods::as(fx$scee, "SingleCellExperiment"), exposures,
            sample_col = "donor_id")
    }
    select <- function(scee, method) {
        suppressMessages(run_exposure_network(scee, fx$metab,
            celltype = "Mono", exposure = "PM2.5", covariates = "age",
            selection_method = method, network_method = "block_glasso",
            transform = "log_cpm"))
    }
    expect_error(select(with_exposure(4L), "univariate"),
        "needs at least 5 donors with a value of the exposure; 4 have one",
        fixed = TRUE)
    expect_error(select(with_exposure(8L), "stability"),
        "needs at least 10 donors with an observed exposure", fixed = TRUE)
})

test_that("an exposure that the covariates give leaves nothing to select", {
    fx <- .covariate_fixture()
    select <- function(exposures, covariates) {
        scee <- build_scee(methods::as(fx$scee, "SingleCellExperiment"),
            exposures, sample_col = "donor_id")
        suppressMessages(run_exposure_network(scee, fx$metab,
            celltype = "Mono", exposure = "PM2.5", covariates = covariates,
            selection_method = "univariate", network_method = "block_glasso",
            transform = "log_cpm"))
    }
    exposures <- exposureData(fx$scee)
    ## the exposure on another scale, given as a covariate
    rescaled <- cbind(exposures, log_PM2.5 = log(exposures[, "PM2.5"] + 10))
    expect_error(select(rescaled, "log_PM2.5"), paste("does not vary across",
        "the donors with a value of it once the covariates are adjusted",
        "for, so nothing can be selected by it."), fixed = TRUE)
    ## an exposure with one value
    constant <- exposures
    constant[, "PM2.5"] <- 3
    expect_error(select(constant, NULL), paste("does not vary across the",
        "donors with a value of it, so nothing can be selected by it."),
        fixed = TRUE)
})
