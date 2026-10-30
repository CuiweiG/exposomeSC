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

test_that("plot_celltype_network draws visible, correctly labelled edges", {
    skip_if_not_installed("ggplot2")
    skip_if_not_installed("ggraph")
    skip_if_not_installed("igraph")
    skip_if_not_installed("tidygraph")

    features <- c("G1", "G2", "M1")
    precision <- diag(3)
    ## partial precisions of the size real networks produce
    precision[1, 3] <- precision[3, 1] <- -0.04
    precision[2, 3] <- precision[3, 2] <- 0.03
    adjacency <- (precision != 0) * 1L
    diag(adjacency) <- 0L
    dimnames(precision) <- dimnames(adjacency) <- list(features, features)
    net <- new("CelltypeNetworkResult",
        precision_matrix = precision,
        adjacency_matrix = adjacency,
        stability_scores = matrix(nrow = 0, ncol = 0),
        node_info = S4Vectors::DataFrame(feature = features,
            omic_layer = c("transcript", "transcript", "metabolite")),
        celltype = "Test", method = "precomputed",
        metadata = list(n_donors = 10L))

    plot <- plot_celltype_network(net, layout = "circle")
    built <- ggplot2::ggplot_build(plot)
    edge_layer <- built$data[[1]]
    expect_true(all(edge_layer$edge_alpha >= 0.35 - 1e-8))
    ## both edges are cross-omic, so that is the only label shown
    colour_scale <- built$plot$scales$get_scales("edge_colour")
    expect_identical(as.character(colour_scale$get_labels()), "Cross-omic")
    ## the subtitle counts the edges drawn, even without stored counts
    expect_match(plot$labels$subtitle, "2 edges (2 cross-omic)", fixed = TRUE)
})

test_that("plot_edge_stability handles no stability", {
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

    p <- plot_edge_stability(net)
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

.stability_network <- function(layers) {
    features <- c("G1", "G2", "M1")
    stability <- matrix(0, 3, 3, dimnames = list(features, features))
    stability[1, 3] <- stability[3, 1] <- 0.85
    stability[1, 2] <- stability[2, 1] <- 0.4
    methods::new(
        "CelltypeNetworkResult",
        precision_matrix = diag(3),
        adjacency_matrix = (stability > 0.5) * 1L,
        stability_scores = stability,
        node_info = S4Vectors::DataFrame(
            feature = features,
            omic_layer = layers
        ),
        celltype = "Monocyte",
        method = "precomputed",
        metadata = list()
    )
}

test_that("plot_edge_stability labels a single edge type correctly", {
    skip_if_not_installed("ggplot2")
    ## every edge joins the two layers
    net <- .stability_network(c("transcript", "metabolite", "metabolite"))
    p <- plot_edge_stability(net)
    guide <- ggplot2::get_guide_data(p, "fill")
    expect_identical(as.character(guide$.label), "Cross-omic")
    bars <- ggplot2::layer_data(p, 1L)
    expect_true(all(bars$fill == "#E69F00"))

    ## every edge stays within one layer
    net <- .stability_network(rep("transcript", 3))
    p <- plot_edge_stability(net)
    guide <- ggplot2::get_guide_data(p, "fill")
    expect_identical(as.character(guide$.label), "Within-omic")
    expect_true(all(ggplot2::layer_data(p, 1L)$fill == "#0072B2"))
})

test_that("plot_edge_stability matches colours to edge types", {
    skip_if_not_installed("ggplot2")
    net <- .stability_network(c("transcript", "transcript", "metabolite"))
    p <- plot_edge_stability(net)
    bars <- ggplot2::layer_data(p, 1L)
    ## bars are drawn from the top edge (G1 - M1, cross-omic, 0.85)
    expect_identical(bars$fill[bars$x == 0.85], "#E69F00")
    expect_identical(bars$fill[bars$x == 0.4], "#0072B2")
    guide <- ggplot2::get_guide_data(p, "fill")
    expect_identical(as.character(guide$.label),
        c("Within-omic", "Cross-omic"))
})

test_that("plot_edge_stability draws a reference line only on request", {
    skip_if_not_installed("ggplot2")
    net <- .stability_network(c("transcript", "transcript", "metabolite"))
    geoms <- function(p) {
        vapply(p$layers, function(l) class(l$geom)[1], character(1))
    }
    expect_false("GeomVline" %in% geoms(plot_edge_stability(net)))
    p <- plot_edge_stability(net, threshold = 0.8)
    expect_true("GeomVline" %in% geoms(p))
    vline <- ggplot2::layer_data(p, which(geoms(p) == "GeomVline"))
    expect_equal(vline$xintercept, 0.8)
    expect_error(plot_edge_stability(net, threshold = 1.5), "threshold")
    expect_error(plot_edge_stability(net, threshold = c(0.5, 0.8)),
        "threshold")
})

test_that("plot_temporal_dynamics labels presence by name", {
    skip_if_not_installed("ggplot2")
    edge_dynamics <- S4Vectors::DataFrame(
        node_i = c("G1", "G2"),
        node_j = c("M1", "M1"),
        cross_omic = c(TRUE, TRUE),
        status = c("persistent", "persistent"),
        present_baseline = c(TRUE, TRUE),
        present_follow_up = c(TRUE, TRUE)
    )
    temporal <- methods::new(
        "TemporalNetwork",
        networks = list(),
        edge_dynamics = edge_dynamics,
        timepoints = c("baseline", "follow_up")
    )
    p <- plot_temporal_dynamics(temporal)
    tiles <- ggplot2::layer_data(p, 1L)
    expect_true(all(tiles$fill == "#0072B2"))
    guide <- ggplot2::get_guide_data(p, "fill")
    expect_identical(as.character(guide$.label), c("Absent", "Present"))
    expect_identical(guide$fill, c("#E5E5E5", "#0072B2"))
})
