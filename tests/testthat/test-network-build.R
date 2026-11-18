# tests/testthat/test-network-build.R

test_that("run_celltype_network works with block_glasso", {
    skip_if_not_installed("glasso")
    skip_if_not_installed("SingleCellExperiment")
    skip_if_not_installed("S4Vectors")

    library(SingleCellExperiment)
    library(S4Vectors)

    set.seed(42)
    n_donors <- 20
    n_cells_per <- 30
    n_genes <- 30
    n_metab <- 5

    # Build mock SCE - all cells are Mono for this test
    donor_ids <- paste0("D", seq_len(n_donors))
    cell_ids <- paste0("c", seq_len(n_donors * n_cells_per))
    donors <- rep(donor_ids, each = n_cells_per)
    ctypes <- rep("Mono", length.out =
                      n_donors * n_cells_per)

    counts <- matrix(rpois(n_genes * length(cell_ids), 10),
                     nrow = n_genes,
                     dimnames = list(paste0("G", seq_len(n_genes)),
                                    cell_ids))
    sce <- SingleCellExperiment(
        assays = list(counts = counts),
        colData = DataFrame(
            cell_id = cell_ids,
            donor_id = donors,
            cell_type = ctypes))

    exp_mat <- matrix(rnorm(n_donors * 2), nrow = n_donors,
                      dimnames = list(donor_ids,
                                      c("PM2.5", "Pb")))

    scee <- build_scee(sce, exp_mat, sample_col = "donor_id")

    # Metabolites
    metab <- matrix(rnorm(n_donors * n_metab), nrow = n_donors,
                    dimnames = list(donor_ids,
                                   paste0("M", seq_len(n_metab))))

    # Build network
    net <- run_celltype_network(
        scee, metab,
        celltype = "Mono",
        celltype_col = "cell_type",
        method = "block_glasso",
        stability = FALSE,
        top_var_genes = 10, transform = "log_cpm")

    expect_s4_class(net, "CelltypeNetworkResult")
    expect_equal(networkCelltype(net), "Mono")
    expect_equal(networkMethod(net), "block_glasso")
    expect_true(nrow(networkPrecision(net)) > 0)
    expect_equal(nrow(networkPrecision(net)),
                 ncol(networkPrecision(net)))
    expect_equal(nrow(networkNodes(net)),
                 nrow(networkPrecision(net)))

    # Check node_info
    expect_true("omic_layer" %in% colnames(networkNodes(net)))
    expect_true(any(networkNodes(net)$omic_layer == "transcript"))
    expect_true(any(networkNodes(net)$omic_layer == "metabolite"))

    # Show method
    expect_output(show(net), "CelltypeNetworkResult")
})

test_that("coglasso networks return edge-level selection frequencies", {
    skip_if_not_installed("coglasso", minimum_version = "1.1.0")
    skip_if_not_installed("SingleCellExperiment")

    set.seed(7)
    donor_ids <- paste0("D", seq_len(20))
    donors <- rep(donor_ids, each = 20)
    counts <- matrix(stats::rpois(20 * length(donors), 10), nrow = 20,
        dimnames = list(paste0("G", seq_len(20)),
                        paste0("c", seq_along(donors))))
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = counts),
        colData = S4Vectors::DataFrame(
            cell_id = colnames(counts), donor_id = donors,
            cell_type = "Mono"))
    scee <- build_scee(sce,
        matrix(stats::rnorm(20), ncol = 1,
               dimnames = list(donor_ids, "PM2.5")),
        sample_col = "donor_id")
    metab <- matrix(stats::rnorm(20 * 4), nrow = 20,
        dimnames = list(donor_ids, paste0("M", seq_len(4))))
    build <- function(...) run_celltype_network(scee, metab,
        celltype = "Mono", method = "coglasso",
        nlambda_w = 3L, nlambda_b = 3L, top_var_genes = 6L, ...,
        transform = "log_cpm")

    net <- build(rep_num = 4L, subsample_ratio = 0.7)
    stab <- networkStability(net)
    p <- nrow(networkPrecision(net))
    ## one frequency per pair of features, not a single summary number
    expect_equal(dim(stab), c(p, p))
    expect_identical(dimnames(stab), dimnames(networkPrecision(net)))
    expect_true(isSymmetric(unname(stab)))
    expect_true(all(stab >= 0 & stab <= 1))
    ## each frequency is a count out of rep_num subsamples
    expect_true(all(abs(stab * 4 - round(stab * 4)) < 1e-8))
    expect_equal(metadata(net)$selection$rep_num, 4L)
    expect_equal(metadata(net)$selection$subsample_ratio, 0.7)
    expect_length(metadata(net)$selection$stars_variability, 1L)

    no_stability <- build(rep_num = 4L, stability = FALSE)
    expect_equal(nrow(networkStability(no_stability)), 0L)
    expect_error(build(subsample_ratio = 1), "subsample_ratio")
    expect_error(build(rep_num = 1), "rep_num")
})

test_that("a precomputed network is accepted as a named list", {
    skip_if_not_installed("SingleCellExperiment")
    skip_if_not_installed("S4Vectors")

    set.seed(11)
    donor_ids <- paste0("D", seq_len(6L))
    donors <- rep(donor_ids, each = 5L)
    cell_ids <- paste0("c", seq_along(donors))
    counts <- matrix(stats::rpois(4L * length(cell_ids), 10), nrow = 4L,
        dimnames = list(paste0("G", seq_len(4L)), cell_ids))
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = counts),
        colData = S4Vectors::DataFrame(
            cell_id = cell_ids, donor_id = donors, cell_type = "Mono"))
    scee <- build_scee(sce,
        matrix(stats::rnorm(6L), ncol = 1L,
               dimnames = list(donor_ids, "PM2.5")),
        sample_col = "donor_id")

    adj <- matrix(0L, 3L, 3L)
    adj[1L, 2L] <- adj[2L, 1L] <- 1L
    adj[2L, 3L] <- adj[3L, 2L] <- 1L
    prec <- diag(3L)
    prec[1L, 2L] <- prec[2L, 1L] <- -0.4
    prec[2L, 3L] <- prec[3L, 2L] <- 0.3
    precomputed <- list(
        adjacency = adj,
        precision = prec,
        stability = matrix(0.5, 3L, 3L),
        feature_names = c("G1", "G2", "M1"),
        feature_blocks = c("transcript", "transcript", "metabolite"))

    expect_message(
        run_celltype_network(scee, celltype = "Mono",
                             precomputed_network = precomputed),
        "Wrapping precomputed network")
    net <- suppressMessages(
        run_celltype_network(scee, celltype = "Mono",
                             precomputed_network = precomputed))

    expect_s4_class(net, "CelltypeNetworkResult")
    expect_equal(networkCelltype(net), "Mono")
    expect_equal(networkMethod(net), "precomputed")
    expect_equal(metadata(net)$source, "list")
    expect_equal(metadata(net)$n_features, 3L)
    ## two undirected edges, counted once each
    expect_equal(metadata(net)$n_edges, 2L)
    expect_identical(networkNodes(net)$feature, c("G1", "G2", "M1"))
    expect_identical(networkNodes(net)$omic_layer,
                     c("transcript", "transcript", "metabolite"))
    expect_identical(networkNodes(net)$block, c(1L, 1L, 2L))
    ## the summary, the plot and the comparison read the omic layer
    expect_output(show(net), "2 transcripts, 1 metabolites")
    expect_output(show(net), "1 cross-omic")
    other <- suppressMessages(run_celltype_network(scee, celltype = "NK",
        precomputed_network = precomputed))
    comparison <- run_comparative_network(net, other)
    expect_s4_class(comparison, "NetworkComparison")
    expect_equal(nrow(sharedEdges(comparison)), 2L)
    if (requireNamespace("ggplot2", quietly = TRUE) &&
            requireNamespace("igraph", quietly = TRUE) &&
            requireNamespace("ggraph", quietly = TRUE)) {
        expect_s3_class(plot_celltype_network(net), "ggplot")
    }
    expect_identical(dimnames(networkAdjacency(net)),
                     list(c("G1", "G2", "M1"), c("G1", "G2", "M1")))
    expect_identical(dimnames(networkPrecision(net)),
                     dimnames(networkAdjacency(net)))
    ## metabolites and the estimation settings are not consulted
    expect_identical(
        suppressMessages(run_celltype_network(scee, celltype = "Mono",
            precomputed_network = precomputed, method = "coglasso")),
        net)

    ## optional elements may be omitted
    minimal <- suppressMessages(run_celltype_network(scee, celltype = "Mono",
        precomputed_network = list(adjacency = adj)))
    expect_equal(nrow(networkStability(minimal)), 0L)
    expect_true(all(is.na(networkPrecision(minimal))))
    expect_identical(networkNodes(minimal)$feature, paste0("F", seq_len(3L)))
    expect_identical(networkNodes(minimal)$omic_layer, rep("unknown", 3L))
    expect_identical(networkNodes(minimal)$block, rep(1L, 3L))

    expect_error(
        run_celltype_network(scee, celltype = "Mono",
                             precomputed_network = list(precision = prec)),
        "\\$adjacency")
    expect_error(
        run_celltype_network(scee, celltype = "Mono",
                             precomputed_network = seq_len(3L)),
        "must be a list")
    expect_error(
        run_celltype_network(scee, celltype = "Mono",
            precomputed_network = list(adjacency = adj,
                feature_names = c("A", "B"))),
        "unique, non-missing and non-empty name")
    expect_error(
        run_celltype_network(scee, celltype = "Mono",
            precomputed_network = list(adjacency = adj,
                feature_names = c("A", "A", "B"))),
        "unique, non-missing and non-empty name")
    expect_error(
        run_celltype_network(scee, celltype = "Mono",
            precomputed_network = list(adjacency = adj,
                feature_names = c("A", "", "B"))),
        "unique, non-missing and non-empty name")
    ## an empty adjacency gets no default names and an empty network
    expect_message(
        empty <- run_celltype_network(scee, celltype = "Mono",
            precomputed_network = list(adjacency = matrix(0L, 0, 0))),
        "p=0 features, 0 edges")
    expect_identical(nrow(networkNodes(empty)), 0L)
    expect_error(
        run_celltype_network(scee, celltype = "Mono",
            precomputed_network = list(adjacency = adj,
                feature_blocks = c("transcript", NA, "metabolite"))),
        "no missing or empty labels")
    expect_error(
        run_celltype_network(scee, celltype = "Mono",
            precomputed_network = list(adjacency = adj,
                block_sizes = c(transcript = 2L, metabolite = NA))),
        "whole numbers named by omic layer")
    expect_error(
        run_celltype_network(scee, celltype = "Mono",
            precomputed_network = list(adjacency = adj,
                block_sizes = c(transcript = 2L, 1L))),
        "whole numbers named by omic layer")
    expect_error(
        run_celltype_network(scee, celltype = "Mono",
            precomputed_network = list(adjacency = adj,
                precision = diag(2L))),
        "dimensions of 'adjacency'")
    expect_error(
        run_celltype_network(scee, celltype = "Mono",
            precomputed_network = list(adjacency = adj,
                block_sizes = c(transcript = 1.5, metabolite = 1.5))),
        "whole numbers named by omic layer")
    from_factor <- suppressMessages(run_celltype_network(scee,
        celltype = "Mono", precomputed_network = list(adjacency = adj,
            feature_names = factor(c("G1", "G2", "M1")))))
    expect_identical(networkNodes(from_factor)$feature, c("G1", "G2", "M1"))
    expect_error(
        run_celltype_network(scee, celltype = "Mono",
            precomputed_network = list(adjacency = adj,
                block_sizes = c(transcript = 2L, metabolite = 2L))),
        "label every feature")
    expect_error(
        run_celltype_network(scee, celltype = "Mono",
            precomputed_network = list(adjacency = matrix(0L, 2L, 3L))),
        "must be square")
})

test_that("CelltypeNetworkResult S4 class validation", {
    # Empty object
    net <- new("CelltypeNetworkResult")
    expect_s4_class(net, "CelltypeNetworkResult")

    # Invalid: non-square precision
    expect_error(
        validObject(new("CelltypeNetworkResult",
            precision_matrix = matrix(1:6, 2, 3),
            adjacency_matrix = matrix(0, 2, 3))),
        "square")
})

test_that("network-utils internal helpers work", {
    # .precision_to_adjacency
    prec <- matrix(c(1, 0.5, 0, 0.5, 1, 0.3, 0, 0.3, 1),
                   3, 3)
    adj <- exposomeSC:::.precision_to_adjacency(prec)
    expect_equal(adj[1, 2], 1)
    expect_equal(adj[1, 3], 0)
    expect_equal(adj[2, 3], 1)
    expect_equal(diag(adj), c(0, 0, 0))

    # .count_edge_types
    ni <- S4Vectors::DataFrame(
        feature = c("g1", "g2", "m1"),
        omic_layer = c("transcript", "transcript", "metabolite"),
        block = c(1L, 1L, 2L))
    counts <- exposomeSC:::.count_edge_types(adj, ni)
    expect_equal(counts$total, 2)
    expect_equal(counts$within_transcript, 1)
    expect_equal(counts$cross_omic, 1)
})

## Twenty donors with one cell type, thirty genes and three metabolites;
## gene G1 has no counts and metabolite M3 the same value in every donor
.constant_fixture <- function(seed = 7) {
    set.seed(seed)
    donor_ids <- paste0("D", seq_len(20))
    donors <- rep(donor_ids, each = 20)
    counts <- matrix(stats::rpois(30 * length(donors), 10), nrow = 30,
        dimnames = list(paste0("G", seq_len(30)),
            paste0("c", seq_along(donors))))
    counts[1, ] <- 0L
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = counts),
        colData = S4Vectors::DataFrame(donor_id = donors,
            cell_type = "Mono"))
    exposures <- matrix(stats::rnorm(40), ncol = 2,
        dimnames = list(donor_ids, c("PM2.5", "Pb")))
    metab <- matrix(stats::rnorm(60), nrow = 20,
        dimnames = list(donor_ids, paste0("M", 1:3)))
    metab[, "M3"] <- 1
    list(scee = build_scee(sce, exposures, sample_col = "donor_id"),
        metab = metab)
}

test_that("constant features are left out of a cell-type network", {
    skip_if_not_installed("glasso")
    fx <- .constant_fixture()
    msgs <- character()
    net <- withCallingHandlers(
        run_celltype_network(fx$scee, fx$metab, celltype = "Mono",
            method = "block_glasso", transform = "log_cpm"),
        message = function(m) {
            msgs <<- c(msgs, conditionMessage(m))
            invokeRestart("muffleMessage")
        })
    expect_true(any(grepl("1 constant transcript", msgs)))
    expect_true(any(grepl("1 constant metabolite", msgs)))
    ## the message names them, each after its count
    expect_true(any(grepl(paste("Leaving 1 constant transcript(s) (G1) and",
        "1 constant metabolite(s) (M3) out of the network for 'Mono'."),
        msgs, fixed = TRUE)))
    expect_identical(metadata(net)$constant_features, c("G1", "M3"))
    expect_false(any(c("G1", "M3") %in% networkNodes(net)$feature))
    expect_equal(metadata(net)$n_transcripts, 29L)
    expect_equal(metadata(net)$n_metabolites, 2L)
    expect_equal(nrow(networkNodes(net)), 31L)
    ## with every metabolite constant there is nothing to connect
    all_constant <- fx$metab
    all_constant[, ] <- 1
    expect_error(
        suppressMessages(run_celltype_network(fx$scee, all_constant,
            celltype = "Mono", method = "block_glasso", transform = "log_cpm")),
        "constant across")
})

test_that("constant features are found before the adjustment", {
    skip_if_not_installed("glasso")
    fx <- .constant_fixture()
    net <- suppressMessages(run_celltype_network(fx$scee, fx$metab,
        celltype = "Mono", method = "block_glasso", covariates = "Pb",
        transform = "log_cpm"))
    expect_identical(metadata(net)$constant_features, c("G1", "M3"))
    expect_false("M3" %in% networkNodes(net)$feature)
    all_constant <- fx$metab
    all_constant[, ] <- 1
    expect_error(
        suppressMessages(run_celltype_network(fx$scee, all_constant,
            celltype = "Mono", method = "block_glasso", covariates = "Pb",
            transform = "log_cpm")),
        "constant across")
    ## the selection stage of run_exposure_network() never picks a
    ## constant feature, with covariates as without
    selected <- suppressMessages(run_exposure_network(fx$scee, fx$metab,
        celltype = "Mono", exposure = "PM2.5", covariates = "Pb",
        selection_method = "univariate", fdr_threshold = 1,
        network_method = "block_glasso", transform = "log_cpm"))
    expect_false(any(c("G1", "M3") %in% networkNodes(selected)$feature))
})

test_that("unnamed and missing metabolite columns are reported clearly", {
    skip_if_not_installed("glasso")
    fx <- .constant_fixture()
    unnamed <- fx$metab[, c(3, 1, 2)]
    colnames(unnamed) <- NULL
    net <- suppressMessages(run_celltype_network(fx$scee, unnamed,
        celltype = "Mono", method = "block_glasso", transform = "log_cpm"))
    ## the constant first column is metab_1 in the metadata and in no node
    expect_identical(metadata(net)$constant_features, c("G1", "metab_1"))
    expect_identical(utils::tail(networkNodes(net)$feature, 2L),
        c("metab_2", "metab_3"))
    expect_error(
        suppressMessages(run_celltype_network(fx$scee,
            fx$metab[, 0, drop = FALSE], celltype = "Mono",
            method = "block_glasso", transform = "log_cpm")),
        "has no columns")
})

test_that("run_celltype_network says that it needs the metabolite matrix", {
    fx <- .constant_fixture()
    needed <- paste("'metabolites' must be a matrix with one row for each",
        "donor. It can be left out only when 'precomputed_network' is given.")
    ## its default is NULL, for a call that wraps a precomputed network
    expect_null(formals(run_celltype_network)$metabolites)
    expect_error(run_celltype_network(fx$scee, celltype = "Mono"), needed,
        fixed = TRUE)
    expect_error(
        run_celltype_network(fx$scee, as.data.frame(fx$metab),
            celltype = "Mono"),
        needed, fixed = TRUE)
    expect_error(
        run_celltype_network(fx$scee, fx$metab[, "M1"], celltype = "Mono"),
        needed, fixed = TRUE)
})

test_that("missing or infinite metabolite values are an error", {
    skip_if_not_installed("glasso")
    fx <- .constant_fixture()
    gaps <- fx$metab
    gaps["D3", "M2"] <- NA
    gaps["D5", "M1"] <- Inf
    expect_error(
        suppressMessages(run_celltype_network(fx$scee, gaps,
            celltype = "Mono", method = "block_glasso", transform = "log_cpm")),
        "for the analysed donors in: M1, M2.", fixed = TRUE)
    ## a NaN used to be replaced by the column mean without a word
    not_a_number <- fx$metab
    not_a_number["D7", "M1"] <- NaN
    expect_error(
        suppressMessages(run_celltype_network(fx$scee, not_a_number,
            celltype = "Mono", method = "block_glasso", transform = "log_cpm")),
        "for the analysed donors in: M1.", fixed = TRUE)
    expect_error(
        suppressMessages(run_celltype_network(fx$scee, gaps,
            celltype = "Mono", method = "block_glasso", exposure = "PM2.5",
            transform = "log_cpm")),
        "in: M1, M2.", fixed = TRUE)
    ## the selection of run_exposure_network() used to skip such a
    ## metabolite; a data frame is checked as well
    expect_error(
        suppressMessages(run_exposure_network(fx$scee, as.data.frame(gaps),
            celltype = "Mono", exposure = "PM2.5",
            selection_method = "univariate", fdr_threshold = 1,
            network_method = "block_glasso", transform = "log_cpm")),
        "in: M1, M2.", fixed = TRUE)
    ## a data frame of numbers is taken as its matrix, also when covariates
    ## are regressed out; one with a column of another kind is refused
    adjusted <- function(metabolites) {
        suppressMessages(run_exposure_network(fx$scee, metabolites,
            celltype = "Mono", exposure = "PM2.5", covariates = "Pb",
            selection_method = "univariate", fdr_threshold = 1,
            network_method = "block_glasso", transform = "log_cpm"))
    }
    from_matrix <- adjusted(fx$metab)
    from_frame <- adjusted(as.data.frame(fx$metab))
    expect_identical(networkNodes(from_frame), networkNodes(from_matrix))
    expect_identical(networkPrecision(from_frame),
        networkPrecision(from_matrix))
    expect_identical(metadata(from_matrix)$covariates, "Pb")
    words <- as.data.frame(fx$metab)
    words$M2 <- as.character(words$M2)
    expect_error(adjusted(words), "'metabolites' must hold numbers only.",
        fixed = TRUE)
    ## a missing value of a donor that is not analysed is no error: D99 has
    ## no cells, and D3 is left out for its missing covariate
    exposures <- exposureData(fx$scee)
    exposures["D3", "Pb"] <- NA
    scee <- build_scee(methods::as(fx$scee, "SingleCellExperiment"),
        exposures, sample_col = "donor_id")
    others <- rbind(fx$metab, D99 = c(NA, 1, 1))
    others["D3", "M2"] <- NA
    net <- suppressMessages(run_celltype_network(scee, others,
        celltype = "Mono", method = "block_glasso", covariates = "Pb",
        transform = "log_cpm"))
    expect_equal(metadata(net)$n_donors, 19L)
})

test_that("unknown exposure or covariate names are an error", {
    skip_if_not_installed("glasso")
    fx <- .constant_fixture()
    expect_error(
        suppressMessages(run_celltype_network(fx$scee, fx$metab,
            celltype = "Mono", method = "block_glasso",
            covariates = "not_a_column", transform = "log_cpm")),
        "not_a_column")
    expect_error(
        suppressMessages(run_celltype_network(fx$scee, fx$metab,
            celltype = "Mono", method = "block_glasso", exposure = "typo",
            transform = "log_cpm")),
        "typo")
    expect_error(
        suppressMessages(run_exposure_network(fx$scee, fx$metab,
            celltype = "Mono", exposure = "PM2.5", covariates = "typo",
            selection_method = "univariate",
            network_method = "block_glasso", transform = "log_cpm")),
        "typo")
    expect_error(
        suppressMessages(run_exposure_network(fx$scee, fx$metab,
            celltype = "Mono", exposure = "typo",
            selection_method = "univariate",
            network_method = "block_glasso", transform = "log_cpm")),
        "typo")
    ## one unknown name among valid ones, and several unknown names
    expect_error(
        suppressMessages(run_celltype_network(fx$scee, fx$metab,
            celltype = "Mono", method = "block_glasso",
            covariates = c("Pb", "typo"), transform = "log_cpm")),
        "exposureData: typo$")
    expect_error(
        suppressMessages(run_celltype_network(fx$scee, fx$metab,
            celltype = "Mono", method = "block_glasso", exposure = "x1",
            covariates = c("Pb", "x2"), transform = "log_cpm")),
        "exposureData: x1, x2$")
    ## a valid adjustment is recorded as before
    net <- suppressMessages(run_celltype_network(fx$scee, fx$metab,
        celltype = "Mono", method = "block_glasso", exposure = "PM2.5",
        covariates = "Pb", transform = "log_cpm"))
    expect_identical(metadata(net)$exposure, "PM2.5")
    expect_identical(metadata(net)$covariates, "Pb")
})

test_that("the network data are adjusted for variables of any name", {
    skip_if_not_installed("glasso")
    fx <- .constant_fixture()
    with_names <- function(names) {
        scee <- fx$scee
        exposures <- exposureData(scee)
        colnames(exposures) <- names
        exposureData(scee) <- exposures
        scee
    }
    network <- function(names, ...) {
        suppressMessages(run_celltype_network(with_names(names), fx$metab,
            celltype = "Mono", method = "block_glasso", transform = "log_cpm",
            ...))
    }
    reference <- network(c("PM2.5", "Pb"), exposure = "PM2.5",
        covariates = "Pb")
    ## a formula would read "..." and "..1" as arguments of the calling
    ## function; the adjustment builds none
    for (name in c("...", "..1", ".", "age-years")) {
        as_covariate <- network(c("PM2.5", name), exposure = "PM2.5",
            covariates = name)
        expect_identical(networkPrecision(as_covariate),
            networkPrecision(reference))
        expect_identical(metadata(as_covariate)$covariates, name)
        as_exposure <- network(c(name, "Pb"), exposure = name,
            covariates = "Pb")
        expect_identical(networkPrecision(as_exposure),
            networkPrecision(reference))
    }
    selected <- function(names, covariate) {
        suppressMessages(run_exposure_network(with_names(names), fx$metab,
            celltype = "Mono", exposure = "PM2.5", covariates = covariate,
            selection_method = "univariate", fdr_threshold = 1,
            network_method = "block_glasso", transform = "log_cpm"))
    }
    expect_identical(networkPrecision(selected(c("PM2.5", "..1"), "..1")),
        networkPrecision(selected(c("PM2.5", "Pb"), "Pb")))

    ## the residuals are those of lm() for each column, with their names,
    ## also for a single column
    residualize <- exposomeSC:::.residualize_matrix
    set.seed(1)
    values <- matrix(stats::rnorm(40), 20, 2,
        dimnames = list(paste0("D", 1:20), c("f1", "f2")))
    adjustment <- matrix(stats::rnorm(40), 20, 2,
        dimnames = list(NULL, c("...", "..1")))
    by_lm <- vapply(1:2, function(j) {
        unname(stats::residuals(stats::lm(values[, j] ~ adjustment)))
    }, numeric(20))
    residuals <- residualize(values, adjustment)
    expect_equal(unname(residuals), by_lm, tolerance = 1e-12)
    expect_identical(dimnames(residuals), dimnames(values))
    one <- residualize(values[, "f2", drop = FALSE], adjustment)
    expect_identical(dimnames(one), list(rownames(values), "f2"))
    expect_equal(unname(one[, 1]), by_lm[, 2], tolerance = 1e-12)
})

test_that("donors with a missing covariate are left out of the network", {
    skip_if_not_installed("glasso")
    fx <- .constant_fixture()
    exposures <- exposureData(fx$scee)
    exposures["D3", "Pb"] <- NA
    scee <- build_scee(methods::as(fx$scee, "SingleCellExperiment"),
        exposures, sample_col = "donor_id")
    msgs <- character()
    net <- withCallingHandlers(
        run_celltype_network(scee, fx$metab, celltype = "Mono",
            method = "block_glasso", covariates = "Pb", transform = "log_cpm"),
        message = function(m) {
            msgs <<- c(msgs, conditionMessage(m))
            invokeRestart("muffleMessage")
        })
    expect_true(any(grepl("1 donor(s) with a missing", msgs, fixed = TRUE)))
    expect_equal(metadata(net)$n_donors, 19L)
    expect_false("D3" %in% metadata(net)$donors)
    msgs <- character()
    selected <- withCallingHandlers(
        run_exposure_network(scee, fx$metab, celltype = "Mono",
            exposure = "PM2.5", covariates = "Pb",
            selection_method = "univariate", fdr_threshold = 1,
            network_method = "block_glasso", transform = "log_cpm"),
        message = function(m) {
            msgs <<- c(msgs, conditionMessage(m))
            invokeRestart("muffleMessage")
        })
    expect_equal(metadata(selected)$n_donors, 19L)
    ## only the covariate is adjusted for here
    expect_true(any(grepl("1 donor(s) with a missing covariate value",
        msgs, fixed = TRUE)))
    ## too few donors with values of the adjustment variables
    exposures[paste0("D", 1:16), "Pb"] <- NA
    sparse <- build_scee(methods::as(fx$scee, "SingleCellExperiment"),
        exposures, sample_col = "donor_id")
    expect_error(
        suppressMessages(run_celltype_network(sparse, fx$metab,
            celltype = "Mono", method = "block_glasso", covariates = "Pb",
            transform = "log_cpm")),
        "Only 4 of the donors shared between pseudobulk and metabolites")
})

test_that("the network functions take one cell type and one exposure name", {
    skip_if_not_installed("glasso")
    set.seed(4)
    donors <- sprintf("D%02d", 1:20)
    donor <- rep(donors, each = 20)
    counts <- matrix(stats::rpois(30 * length(donor), 20), nrow = 30,
        dimnames = list(paste0("G", 1:30), paste0("c", seq_along(donor))))
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = counts),
        colData = S4Vectors::DataFrame(donor_id = donor, cell_type = "Mono",
            cluster = 1L))
    exposures <- cbind(E1 = stats::rnorm(20), E2 = stats::rnorm(20))
    rownames(exposures) <- donors
    scee <- build_scee(sce, exposures, sample_col = "donor_id")
    metab <- matrix(stats::rnorm(60), nrow = 20,
        dimnames = list(donors, paste0("M", 1:3)))
    network <- function(...) {
        suppressMessages(run_celltype_network(scee, metab,
            method = "block_glasso", ..., transform = "log_cpm"))
    }
    expect_error(network(celltype = c("Mono", "NK")),
        "'celltype' must be a single name")
    expect_identical(
        networkCelltype(network(celltype = 1, celltype_col = "cluster")), "1")
    ## a factor exposure is its label, not the column of its integer code
    expect_equal(
        networkPrecision(network(celltype = "Mono", exposure = factor("E2"))),
        networkPrecision(network(celltype = "Mono", exposure = "E2")))
    expect_error(
        suppressMessages(run_exposure_network(scee, metab,
            celltype = c("Mono", "NK"), exposure = "E1",
            selection_method = "univariate", network_method = "block_glasso",
            transform = "log_cpm")),
        "'celltype' must be a single name")
})

## Twelve donors with twelve cells of one cell type each, two metabolites
## and `n_genes` genes with about 120 counts per donor
.transform_fixture <- function(n_genes) {
    set.seed(11)
    donor_ids <- sprintf("D%02d", 1:12)
    donors <- rep(donor_ids, each = 12)
    counts <- matrix(stats::rpois(n_genes * length(donors), 10),
        nrow = n_genes, dimnames = list(paste0("G", seq_len(n_genes)),
            paste0("c", seq_along(donors))))
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = counts),
        colData = S4Vectors::DataFrame(donor_id = donors,
            cell_type = "Mono"))
    exposures <- matrix(stats::rnorm(12), ncol = 1,
        dimnames = list(donor_ids, "PM2.5"))
    metab <- matrix(stats::rnorm(24), nrow = 12,
        dimnames = list(donor_ids, c("M1", "M2")))
    pseudobulk <- vapply(split(seq_along(donors), donors),
        function(cells) rowSums(counts[, cells]), numeric(n_genes))
    list(scee = build_scee(sce, exposures, sample_col = "donor_id"),
        metab = metab, pseudobulk = pseudobulk)
}

## The transcripts of a network of the six most variable genes
.transform_network <- function(fx, ...) {
    net <- suppressMessages(run_celltype_network(fx$scee, fx$metab,
        celltype = "Mono", method = "block_glasso", top_var_genes = 6L, ...))
    nodes <- networkNodes(net)
    list(net = net,
        transcripts = nodes$feature[nodes$omic_layer == "transcript"])
}

.most_variable <- function(values, n = 6L) {
    names(sort(apply(values, 1, stats::var), decreasing = TRUE))[seq_len(n)]
}

test_that("transform = \"log_cpm\" gives log counts per million", {
    skip_if_not_installed("glasso")
    fx <- .transform_fixture(30L)
    by_cpm <- .transform_network(fx, transform = "log_cpm")
    expect_identical(metadata(by_cpm$net)$transform, "log_cpm")
    log_cpm <- log2(t(t(fx$pseudobulk) / colSums(fx$pseudobulk) * 1e6) + 1)
    expect_setequal(by_cpm$transcripts, .most_variable(log_cpm))
    ## the values themselves, donors by genes
    donors <- colnames(fx$pseudobulk)
    expect_equal(.network_transform(fx$pseudobulk, donors, "log_cpm", "Mono"),
        t(log_cpm))
    ## and in the matrix of the network, where each column is scaled
    assembled <- .assemble_crossomic(fx$scee, fx$metab, "Mono",
        transform = "log_cpm", genes = c("G3", "G17"))
    expect_equal(assembled$data_matrix[, c("G3", "G17")],
        scale(t(log_cpm)[, c("G3", "G17")]), ignore_attr = TRUE)
    expect_error(.transform_network(fx, transform = "rlog"),
        "'arg' should be one of")
    selected <- suppressMessages(run_exposure_network(fx$scee, fx$metab,
        celltype = "Mono", exposure = "PM2.5",
        selection_method = "univariate", fdr_threshold = 1,
        network_method = "block_glasso", transform = "log_cpm"))
    expect_identical(metadata(selected)$transform, "log_cpm")
})

test_that("transform = \"vst\" is DESeq2's transformation or an error", {
    skip_if_not_installed("glasso")
    skip_if_not_installed("DESeq2")
    ## vst() needs 1000 genes with a mean normalised count above 5
    large <- .transform_fixture(1200L)
    by_vst <- .transform_network(large)
    expect_identical(metadata(by_vst$net)$transform, "vst")
    expect_identical(.transform_network(large, transform = "vst")$net,
        by_vst$net)
    vst_values <- suppressMessages({
        dds <- DESeq2::DESeqDataSetFromMatrix(large$pseudobulk,
            data.frame(donor = colnames(large$pseudobulk)), design = ~ 1)
        SummarizedExperiment::assay(DESeq2::vst(dds, blind = TRUE))
    })
    expect_setequal(by_vst$transcripts, .most_variable(vst_values))
    ## the values themselves, which are not those of the other transform
    donors <- colnames(large$pseudobulk)
    transformed <- suppressMessages(
        .network_transform(large$pseudobulk, donors, "vst", "Mono"))
    expect_equal(transformed, t(vst_values))
    expect_false(isTRUE(all.equal(transformed,
        .network_transform(large$pseudobulk, donors, "log_cpm", "Mono"))))

    ## with fewer genes the call stops; it does not use log-CPM values
    small <- .transform_fixture(30L)
    failed <- paste0("^DESeq2::vst\\(\\) failed for 'Mono': fewer than 1000 ",
        "genes have a mean normalised count above 5\\. Use ",
        "transform = \"log_cpm\" for counts that vst\\(\\) cannot ",
        "transform\\.$")
    expect_error(.transform_network(small), failed)
    expect_error(.transform_network(small, transform = "vst"), failed)
    expect_error(
        suppressMessages(run_exposure_network(small$scee, small$metab,
            celltype = "Mono", exposure = "PM2.5",
            selection_method = "univariate", fdr_threshold = 1,
            network_method = "block_glasso")),
        failed)
})

test_that("transform = \"vst\" stops when DESeq2 is not installed", {
    skip_if_not_installed("glasso")
    fx <- .transform_fixture(30L)
    testthat::local_mocked_bindings(
        .is_installed = function(package) package != "DESeq2",
        .package = "exposomeSC")
    missing_package <- paste0("Package 'DESeq2' is required for ",
        "transform = \"vst\". Install it with ",
        "BiocManager::install('DESeq2'), or use transform = \"log_cpm\".")
    expect_error(.transform_network(fx), missing_package, fixed = TRUE)
    expect_error(
        suppressMessages(run_exposure_network(fx$scee, fx$metab,
            celltype = "Mono", exposure = "PM2.5",
            selection_method = "univariate", fdr_threshold = 1,
            network_method = "block_glasso")),
        missing_package, fixed = TRUE)
    ## log-CPM values need no suggested package
    expect_identical(
        metadata(.transform_network(fx, transform = "log_cpm")$net)$transform,
        "log_cpm")
})

## Sixteen donors with twelve cells of each of two cell types, 30 genes and
## three metabolites
.genes_fixture <- function() {
    set.seed(21)
    donor_ids <- sprintf("D%02d", 1:16)
    donors <- rep(rep(donor_ids, each = 12), times = 2)
    types <- rep(c("Mono", "NK"), each = 16 * 12)
    counts <- matrix(stats::rpois(30 * length(donors), 10), nrow = 30,
        dimnames = list(paste0("G", 1:30), paste0("c", seq_along(donors))))
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = counts),
        colData = S4Vectors::DataFrame(donor_id = donors, cell_type = types))
    exposures <- matrix(stats::rnorm(16), ncol = 1,
        dimnames = list(donor_ids, "PM2.5"))
    metab <- matrix(stats::rnorm(48), nrow = 16,
        dimnames = list(donor_ids, paste0("M", 1:3)))
    list(scee = build_scee(sce, exposures, sample_col = "donor_id"),
        metab = metab)
}

test_that("genes names the transcripts of a network after the transform", {
    skip_if_not_installed("glasso")
    fx <- .genes_fixture()
    chosen <- c("G7", "G2", "G19", "G11")
    assemble <- function(x, ...) {
        .assemble_crossomic(x, fx$metab, "Mono", transform = "log_cpm", ...)
    }
    every_gene <- assemble(fx$scee)
    named <- assemble(fx$scee, genes = chosen)
    ## the genes in the order given, then the metabolites
    expect_identical(named$node_info$feature, c(chosen, colnames(fx$metab)))
    expect_identical(named$n_transcripts, 4L)
    ## with the values that they have when every gene is kept
    expect_identical(named$data_matrix[, chosen],
        every_gene$data_matrix[, chosen])
    ## which an object cut down to these genes does not give
    cut_down <- assemble(fx$scee[chosen, ])
    expect_false(isTRUE(all.equal(cut_down$data_matrix[, chosen],
        named$data_matrix[, chosen])))

    network <- function(celltype, ...) {
        suppressMessages(run_celltype_network(fx$scee, fx$metab,
            celltype = celltype, method = "block_glasso",
            transform = "log_cpm", ...))
    }
    transcripts <- function(net) {
        nodes <- networkNodes(net)
        nodes$feature[nodes$omic_layer == "transcript"]
    }
    mono <- network("Mono", genes = chosen)
    nk <- network("NK", genes = chosen)
    expect_identical(transcripts(mono), chosen)
    expect_identical(networkNodes(nk)$feature, networkNodes(mono)$feature)
    ## a factor names the genes by its labels
    expect_identical(transcripts(network("Mono", genes = factor(chosen))),
        chosen)
    ## two cell types with the same genes can be compared; the most variable
    ## genes of each are not the same genes
    expect_s4_class(run_comparative_network(mono, nk), "NetworkComparison")
    expect_error(
        run_comparative_network(network("Mono", top_var_genes = 4L),
            network("NK", top_var_genes = 4L)),
        "All networks must have identical features")

    ## top_var_genes keeps the most variable of the genes named
    pseudobulk <- SummarizedExperiment::assay(
        donor_pseudobulk(fx$scee)[["Mono"]], "counts")
    log_cpm <- log2(t(t(pseudobulk) / colSums(pseudobulk) * 1e6) + 1)
    variance <- apply(log_cpm[chosen, ], 1, stats::var)
    expect_setequal(
        transcripts(network("Mono", genes = chosen, top_var_genes = 2L)),
        names(sort(variance, decreasing = TRUE))[1:2])

    ## names that are not rows of x, or not distinct names, are errors
    expect_error(network("Mono", genes = c("G1", "nope", "G2", "also")),
        "'genes' not found in rownames(x): nope, also", fixed = TRUE)
    not_names <- paste("'genes' must be NULL or a vector of distinct row",
        "names of 'x'.")
    expect_error(network("Mono", genes = c("G1", "G1")), not_names,
        fixed = TRUE)
    expect_error(network("Mono", genes = character()), not_names,
        fixed = TRUE)
    expect_error(network("Mono", genes = 1:3), not_names, fixed = TRUE)
    expect_error(network("Mono", genes = c("G1", NA)), not_names,
        fixed = TRUE)

    ## run_exposure_network() selects among the genes named
    selection <- function(...) {
        suppressMessages(run_exposure_network(fx$scee, fx$metab,
            celltype = "Mono", exposure = "PM2.5",
            selection_method = "univariate", fdr_threshold = 1,
            network_method = "block_glasso", transform = "log_cpm", ...))
    }
    selected <- selection(genes = chosen)
    expect_identical(metadata(selected)$n_total_features, 7L)
    expect_true(all(networkNodes(selected)$feature %in%
        c(chosen, colnames(fx$metab))))
    expect_identical(metadata(selection())$n_total_features, 33L)
    expect_error(selection(genes = "nope"),
        "'genes' not found in rownames(x): nope", fixed = TRUE)
})

test_that("a named gene that is constant in one cell type is left out there", {
    skip_if_not_installed("glasso")
    fx <- .genes_fixture()
    chosen <- c("G7", "G2", "G19", "G11")
    ## G7 without a count in the NK cells
    sce <- methods::as(fx$scee, "SingleCellExperiment")
    counts <- SummarizedExperiment::assay(sce, "counts")
    counts["G7", sce$cell_type == "NK"] <- 0L
    SummarizedExperiment::assay(sce, "counts") <- counts
    scee <- build_scee(sce, exposureData(fx$scee), sample_col = "donor_id")
    network <- function(celltype, genes) {
        run_celltype_network(scee, fx$metab, celltype = celltype,
            method = "block_glasso", transform = "log_cpm", genes = genes)
    }
    mono <- suppressMessages(network("Mono", chosen))
    msgs <- character()
    nk <- withCallingHandlers(network("NK", chosen),
        message = function(m) {
            msgs <<- c(msgs, conditionMessage(m))
            invokeRestart("muffleMessage")
        })
    expect_true(any(grepl(paste("Leaving 1 constant transcript(s) (G7) and",
        "0 constant metabolite(s) out of the network for 'NK'."), msgs,
        fixed = TRUE)))
    expect_identical(metadata(nk)$constant_features, "G7")
    expect_identical(networkNodes(nk)$feature,
        c("G2", "G19", "G11", colnames(fx$metab)))
    ## the comparison says which feature the networks do not share
    expect_error(run_comparative_network(mono, nk),
        paste("Network for 'NK' has different features from 'Mono'.",
            "Only in 'Mono': G7."), fixed = TRUE)
    ## and names the network that has it when that is the second one
    expect_error(run_comparative_network(nk, mono),
        paste("Network for 'Mono' has different features from 'NK'.",
            "Only in 'Mono': G7."), fixed = TRUE)
    ## without that gene the two cell types can be compared
    expect_s4_class(
        suppressMessages(run_comparative_network(
            network("Mono", chosen[-1]), network("NK", chosen[-1]))),
        "NetworkComparison")
    ## the same features in another order are not the same networks either
    expect_error(
        suppressMessages(run_comparative_network(
            network("Mono", chosen[-1]), network("NK", rev(chosen[-1])))),
        "The features are the same, in another order.", fixed = TRUE)

    ## a name that is on two rows of x cannot name a gene
    rownames(sce)[2] <- "G1"
    twice <- build_scee(sce, exposureData(fx$scee), sample_col = "donor_id")
    expect_error(
        run_celltype_network(twice, fx$metab, celltype = "Mono",
            method = "block_glasso", transform = "log_cpm",
            genes = c("G1", "G3")),
        "'genes' that are the name of more than one row of 'x': G1",
        fixed = TRUE)
})

test_that("a note of DESeq2 on another dispersion fit is passed on", {
    quiet <- function() {
        message("converting counts to integer mode")
        1
    }
    expect_no_message(value <- .relay_vst_notes(quiet(), "Mono"))
    expect_identical(value, 1)
    noted <- function() {
        message("-- note: fitType='parametric', but the dispersion trend ",
            "was not well captured by the function: y = a/x + b, and a ",
            "local regression fit was automatically substituted.")
        2
    }
    ## every message that reaches the caller: the note of the package is
    ## replaced by one of ours, which is not repeated
    msgs <- character()
    value <- withCallingHandlers(.relay_vst_notes(noted(), "Mono"),
        message = function(m) {
            msgs <<- c(msgs, conditionMessage(m))
            invokeRestart("muffleMessage")
        })
    expect_length(msgs, 1L)
    expect_identical(msgs,
        paste0("[exposomeSC] DESeq2::vst() for 'Mono': the parametric fit of ",
            "the dispersion trend failed, and DESeq2 used a local fit in ",
            "its place.\n"))
    expect_false(any(grepl("fitType", msgs, fixed = TRUE)))
    expect_identical(value, 2)
})

.coglasso_needed <- paste0("Package 'coglasso' (>= 1.1.0) is required for ",
    "method = \"coglasso\". Install it with install.packages('coglasso'), ",
    "or use method = \"block_glasso\".")

test_that("coglasso is asked for when it is not installed", {
    testthat::local_mocked_bindings(
        .is_installed = function(package) package != "coglasso",
        .package = "exposomeSC")
    expect_error(.require_coglasso("method"), .coglasso_needed, fixed = TRUE)
    expect_error(.require_coglasso("network_method"),
        "or use network_method = \"block_glasso\".", fixed = TRUE)
})

test_that("coglasso older than version 1.1.0 is not used", {
    testthat::local_mocked_bindings(
        .is_installed = function(package) TRUE,
        .installed_version = function(package) package_version("1.0.2"),
        .package = "exposomeSC")
    expect_error(.require_coglasso("method"), .coglasso_needed, fixed = TRUE)
})

test_that("coglasso of version 1.1.0 is accepted", {
    testthat::local_mocked_bindings(
        .is_installed = function(package) TRUE,
        .installed_version = function(package) package_version("1.1.0"),
        .package = "exposomeSC")
    expect_true(.require_coglasso("method"))
})

test_that("run_celltype_network asks for coglasso before using the data", {
    fx <- .constant_fixture()
    testthat::local_mocked_bindings(
        .is_installed = function(package) package != "coglasso",
        .assemble_crossomic = function(...) stop("the data were assembled"),
        .package = "exposomeSC")
    expect_error(
        run_celltype_network(fx$scee, fx$metab, celltype = "Mono",
            method = "coglasso", transform = "log_cpm"),
        .coglasso_needed, fixed = TRUE)
})

test_that("a shortened list of names says how many are left out", {
    expect_identical(.some_names(letters[1:7]), "a, b, c, d, e and 2 more")
    expect_identical(.some_names(letters[1:3]), "a, b, c")
})

test_that("covariates given as a factor are taken by their labels", {
    skip_if_not_installed("glasso")
    fx <- .constant_fixture()
    network <- function(covariates) {
        suppressMessages(run_celltype_network(fx$scee, fx$metab,
            celltype = "Mono", method = "block_glasso",
            covariates = covariates, transform = "log_cpm"))
    }
    by_name <- network("Pb")
    expect_identical(network(factor("Pb")), by_name)
    expect_identical(metadata(by_name)$covariates, "Pb")
    selection <- function(covariates) {
        suppressMessages(run_exposure_network(fx$scee, fx$metab,
            celltype = "Mono", exposure = "PM2.5", covariates = covariates,
            selection_method = "univariate", fdr_threshold = 1,
            network_method = "block_glasso", transform = "log_cpm"))
    }
    selected <- selection("Pb")
    expect_identical(selection(factor("Pb")), selected)
    expect_identical(metadata(selected)$covariates, "Pb")
})

test_that("a name on two rows of x cannot name a node", {
    skip_if_not_installed("glasso")
    fx <- .constant_fixture()
    sce <- methods::as(fx$scee, "SingleCellExperiment")
    rownames(sce)[2] <- "G3"
    twice <- build_scee(sce, exposureData(fx$scee), sample_col = "donor_id")
    repeated <- paste("Two nodes of a network cannot have the same name:",
        "G3. Give the rows of 'x' and the columns of 'metabolites' distinct",
        "names, or name the genes to use with 'genes'.")
    expect_error(
        suppressMessages(run_celltype_network(twice, fx$metab,
            celltype = "Mono", method = "block_glasso", transform = "log_cpm",
            genes = NULL)),
        repeated, fixed = TRUE)
    expect_error(
        suppressMessages(run_exposure_network(twice, fx$metab,
            celltype = "Mono", exposure = "PM2.5", genes = NULL,
            selection_method = "univariate", fdr_threshold = 1,
            network_method = "block_glasso", transform = "log_cpm")),
        repeated, fixed = TRUE)
    ## genes that leave the name out give a network
    net <- suppressMessages(run_celltype_network(twice, fx$metab,
        celltype = "Mono", method = "block_glasso", transform = "log_cpm",
        genes = paste0("G", 4:12)))
    expect_identical(metadata(net)$n_transcripts, 9L)
    ## two metabolites of one name, and a metabolite named as a gene
    same <- function(to, from) {
        metab <- fx$metab
        colnames(metab)[colnames(metab) == from] <- to
        suppressMessages(run_celltype_network(fx$scee, metab,
            celltype = "Mono", method = "block_glasso",
            transform = "log_cpm"))
    }
    expect_error(same("M1", "M2"), "cannot have the same name: M1.",
        fixed = TRUE)
    expect_error(same("G4", "M2"), "cannot have the same name: G4.",
        fixed = TRUE)
})

test_that("a precomputed network is wrapped before the exposure is read", {
    fx <- .constant_fixture()
    adj <- matrix(0L, 3L, 3L)
    adj[1L, 2L] <- adj[2L, 1L] <- 1L
    precomputed <- list(adjacency = adj,
        feature_names = c("G2", "G3", "M1"),
        feature_blocks = c("transcript", "transcript", "metabolite"))
    wrapped <- suppressMessages(run_celltype_network(fx$scee,
        celltype = "Mono", precomputed_network = precomputed))
    expect_identical(
        suppressMessages(run_celltype_network(fx$scee, celltype = "Mono",
            precomputed_network = precomputed, exposure = c("a", "b"))),
        wrapped)
})
