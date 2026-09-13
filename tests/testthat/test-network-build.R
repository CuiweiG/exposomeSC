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
