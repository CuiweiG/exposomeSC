# tests/testthat/test-network-mediation.R

.mediation_fixture <- function(sign = 1) {
    set.seed(11)
    donors <- sprintf("D%02d", 1:30)
    donor <- rep(donors, each = 30)
    exposure <- stats::setNames(stats::rnorm(30), donors)
    mediator <- sign * 2 * exposure + stats::rnorm(30, sd = 0.3)
    counts <- matrix(stats::rpois(10 * length(donor), 20), nrow = 10,
        dimnames = list(paste0("G", 1:10), paste0("c", seq_along(donor))))
    lift <- exp(0.5 * mediator[donor])
    counts[1, ] <- stats::rpois(length(donor), 20 * lift)
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = counts),
        colData = S4Vectors::DataFrame(donor_id = donor, cell_type = "Mono"))
    scee <- build_scee(sce, matrix(exposure, ncol = 1,
        dimnames = list(donors, "E")), sample_col = "donor_id")
    features <- c("G1", "G2", "M1")
    prec <- diag(3)
    prec[1, 3] <- prec[3, 1] <- 0.4
    adj <- (prec != 0) * 1L
    diag(adj) <- 0L
    dimnames(prec) <- dimnames(adj) <- list(features, features)
    net <- new("CelltypeNetworkResult",
        precision_matrix = prec, adjacency_matrix = adj,
        stability_scores = matrix(nrow = 0, ncol = 0),
        node_info = S4Vectors::DataFrame(feature = features,
            omic_layer = c("transcript", "transcript", "metabolite"),
            block = c(1L, 1L, 2L)),
        celltype = "Mono", method = "precomputed",
        metadata = list(n_donors = 30L))
    list(scee = scee, metab = matrix(mediator, ncol = 1,
         dimnames = list(donors, "M1")), net = net)
}

test_that("run_network_mediation warns once that it is experimental", {
    fx <- .mediation_fixture()
    exposomeSC:::.reset_experimental_warnings()
    expect_warning(run_network_mediation(fx$scee, fx$metab, "Mono", "E",
        fx$net, sample_col = "donor_id", n_boot = 50L), "experimental")
    expect_no_warning(run_network_mediation(fx$scee, fx$metab, "Mono", "E",
        fx$net, sample_col = "donor_id", n_boot = 50L))
})

test_that("bootstrap p-value is two-sided for negative indirect effects", {
    fx <- .mediation_fixture(sign = -1)
    res <- suppressWarnings(run_network_mediation(fx$scee, fx$metab, "Mono",
        "E", fx$net, sample_col = "donor_id", n_boot = 200L))
    expect_lt(res$indirect_effect[1], 0)
    expect_lt(res$p_value[1], 0.05)
})

test_that("only topology-guided mediation is offered", {
    fx <- .mediation_fixture()
    expect_error(suppressWarnings(run_network_mediation(fx$scee, fx$metab,
        "Mono", "E", fx$net, sample_col = "donor_id", method = "penalized",
        n_boot = 10L)))
})
