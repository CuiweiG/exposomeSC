# tests/testthat/test-network-build.R

test_that("run_celltype_network works with block_glasso", {
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
        sample_col = "donor_id",
        method = "block_glasso",
        stability = FALSE,
        top_var_genes = 10)

    expect_s4_class(net, "CelltypeNetworkResult")
    expect_equal(net@celltype, "Mono")
    expect_equal(net@method, "block_glasso")
    expect_true(nrow(net@precision_matrix) > 0)
    expect_equal(nrow(net@precision_matrix),
                 ncol(net@precision_matrix))
    expect_equal(nrow(net@node_info),
                 nrow(net@precision_matrix))

    # Check node_info
    expect_true("omic_layer" %in% colnames(net@node_info))
    expect_true(any(net@node_info$omic_layer == "transcript"))
    expect_true(any(net@node_info$omic_layer == "metabolite"))

    # Show method
    expect_output(show(net), "CelltypeNetworkResult")
})

test_that("coglasso networks return edge-level selection frequencies", {
    skip_if_not_installed("coglasso")
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
        celltype = "Mono", sample_col = "donor_id", method = "coglasso",
        nlambda_w = 3L, nlambda_b = 3L, top_var_genes = 6L, ...)

    net <- build(rep_num = 4L, subsample_ratio = 0.7)
    stab <- net@stability_scores
    p <- nrow(net@precision_matrix)
    ## one frequency per pair of features, not a single summary number
    expect_equal(dim(stab), c(p, p))
    expect_identical(dimnames(stab), dimnames(net@precision_matrix))
    expect_true(isSymmetric(unname(stab)))
    expect_true(all(stab >= 0 & stab <= 1))
    ## each frequency is a count out of rep_num subsamples
    expect_true(all(abs(stab * 4 - round(stab * 4)) < 1e-8))
    expect_equal(net@metadata$selection$rep_num, 4L)
    expect_equal(net@metadata$selection$subsample_ratio, 0.7)
    expect_length(net@metadata$selection$stars_variability, 1L)

    expect_equal(nrow(build(rep_num = 4L, stability = FALSE)@stability_scores), 0L)
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
    expect_equal(net@celltype, "Mono")
    expect_equal(net@method, "precomputed")
    expect_equal(net@metadata$source, "list")
    expect_equal(net@metadata$n_features, 3L)
    ## two undirected edges, counted once each
    expect_equal(net@metadata$n_edges, 2L)
    expect_identical(net@node_info$feature, c("G1", "G2", "M1"))
    expect_identical(net@node_info$omic,
                     c("transcript", "transcript", "metabolite"))
    expect_identical(dimnames(net@adjacency_matrix),
                     list(c("G1", "G2", "M1"), c("G1", "G2", "M1")))
    expect_identical(dimnames(net@precision_matrix),
                     dimnames(net@adjacency_matrix))
    ## metabolites and the estimation settings are not consulted
    expect_identical(
        suppressMessages(run_celltype_network(scee, celltype = "Mono",
            precomputed_network = precomputed, method = "coglasso")),
        net)

    ## optional elements may be omitted
    minimal <- suppressMessages(run_celltype_network(scee, celltype = "Mono",
        precomputed_network = list(adjacency = adj)))
    expect_equal(nrow(minimal@stability_scores), 0L)
    expect_true(all(is.na(minimal@precision_matrix)))
    expect_identical(minimal@node_info$feature, paste0("F", seq_len(3L)))
    expect_identical(minimal@node_info$omic, rep("unknown", 3L))

    expect_error(
        run_celltype_network(scee, celltype = "Mono",
                             precomputed_network = list(precision = prec)),
        "\\$adjacency")
    expect_error(
        run_celltype_network(scee, celltype = "Mono",
                             precomputed_network = seq_len(3L)),
        "must be a list")
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
