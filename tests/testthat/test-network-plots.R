# tests/testthat/test-network-plots.R

.make_network_comparison_for_plot <- function(categories = character()) {
    empty_edges <- S4Vectors::DataFrame(
        node_i = character(),
        node_j = character(),
        category = character(),
        max_diff = numeric(),
        p_value = numeric(),
        p_adjusted = numeric()
    )
    make_edges <- function(selected) {
        if (!length(selected)) {
            return(empty_edges)
        }
        S4Vectors::DataFrame(
            node_i = paste0("G", seq_along(selected)),
            node_j = rep("M1", length(selected)),
            category = selected,
            max_diff = rep(0.2, length(selected)),
            p_value = rep(0.1, length(selected)),
            p_adjusted = rep(0.1, length(selected))
        )
    }

    shared <- categories[categories == "shared"]
    differential <- categories[categories != "shared"]
    methods::new(
        "NetworkComparison",
        diff_edges = make_edges(differential),
        shared_edges = make_edges(shared),
        celltype_specific = list(),
        summary = list(
            celltypes = c("A", "B"),
            n_edges = c(A = 2L, B = 1L),
            n_shared = length(shared),
            n_differential = length(differential),
            jaccard = list(A_vs_B = 0.5)
        )
    )
}

test_that("plot functions handle empty networks gracefully", {
    skip_if_not_installed("ggplot2")
    library(S4Vectors)

    net <- new("CelltypeNetworkResult",
        precision_matrix = matrix(c(1, 0, 0, 1), 2, 2),
        adjacency_matrix = matrix(c(0, 0, 0, 0), 2, 2),
        stability_scores = matrix(nrow = 0, ncol = 0),
        node_info = DataFrame(
            feature = c("g1", "m1"),
            omic_layer = c("transcript", "metabolite"),
            block = c(1L, 2L)),
        celltype = "Test",
        method = "block_glasso",
        metadata = list(n_donors = 10L))

    # Should handle gracefully (no edges to plot)
    p <- plot_celltype_network(net)
    expect_s3_class(p, "gg")
})

test_that("plot_stability_surface handles no stability", {
    skip_if_not_installed("ggplot2")
    library(S4Vectors)

    net <- new("CelltypeNetworkResult",
        precision_matrix = matrix(1, 1, 1),
        adjacency_matrix = matrix(0, 1, 1),
        stability_scores = matrix(nrow = 0, ncol = 0),
        node_info = DataFrame(
            feature = "g1", omic_layer = "transcript",
            block = 1L),
        celltype = "Test", method = "test",
        metadata = list())

    p <- plot_stability_surface(net)
    expect_s3_class(p, "gg")
})

test_that("plot_temporal_dynamics handles empty dynamics", {
    skip_if_not_installed("ggplot2")

    tn <- new("TemporalNetwork")
    p <- plot_temporal_dynamics(tn)
    expect_s3_class(p, "gg")
})

test_that("plot_network_comparison handles every edge-category boundary", {
    skip_if_not_installed("ggplot2")

    cases <- list(
        empty = character(),
        shared_only = "shared",
        differential_only = "differential",
        unique_only = "unique_A",
        mixed = c("shared", "differential", "unique_A")
    )
    for (case_name in names(cases)) {
        comparison <- .make_network_comparison_for_plot(cases[[case_name]])
        plot <- plot_network_comparison(comparison, highlight = "all")
        expect_true(inherits(plot, "gg"), info = case_name)
    }
})
