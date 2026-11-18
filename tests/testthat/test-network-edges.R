# tests/testthat/test-network-edges.R
# The edge summaries are written without loops over pairs of nodes. The
# loops they replaced are kept here as the reference: every table must be
# the same, row for row.

.loop_count_edge_types <- function(adjacency, node_info) {
    layers <- node_info$omic_layer
    p <- nrow(adjacency)
    total <- within_tx <- within_met <- cross <- 0L
    for (i in seq_len(p - 1L)) {
        for (j in (i + 1L):p) {
            if (adjacency[i, j] != 0) {
                total <- total + 1L
                if (layers[i] == layers[j]) {
                    if (layers[i] == "transcript")
                        within_tx <- within_tx + 1L
                    else
                        within_met <- within_met + 1L
                } else {
                    cross <- cross + 1L
                }
            }
        }
    }
    list(total = total, within_transcript = within_tx,
        within_metabolite = within_met, cross_omic = cross)
}

.loop_extract_edges <- function(adjacency, node_info, precision = NULL,
                                stability = NULL) {
    p <- nrow(adjacency)
    edges <- list()
    idx <- 0L
    for (i in seq_len(p - 1L)) {
        for (j in (i + 1L):p) {
            if (adjacency[i, j] != 0) {
                idx <- idx + 1L
                edges[[idx]] <- list(
                    node_i = node_info$feature[i],
                    node_j = node_info$feature[j],
                    weight = if (!is.null(precision))
                        precision[i, j] else NA_real_,
                    stability = if (!is.null(stability))
                        stability[i, j] else NA_real_,
                    layer_i = node_info$omic_layer[i],
                    layer_j = node_info$omic_layer[j],
                    cross_omic = node_info$omic_layer[i] !=
                        node_info$omic_layer[j])
            }
        }
    }
    if (length(edges) == 0) {
        return(S4Vectors::DataFrame(node_i = character(),
            node_j = character(), weight = numeric(), stability = numeric(),
            layer_i = character(), layer_j = character(),
            cross_omic = logical()))
    }
    S4Vectors::DataFrame(do.call(rbind.data.frame,
        lapply(edges, as.data.frame, stringsAsFactors = FALSE)))
}

## the edge table of run_comparative_network(), by the loop it used
.loop_comparison_edges <- function(networks) {
    nodes <- networkNodes(networks[[1L]])
    celltypes <- vapply(networks, networkCelltype, character(1))
    partial <- lapply(networks, function(network) {
        theta <- networkPrecision(network)
        scale <- sqrt(diag(theta))
        -theta / outer(scale, scale)
    })
    adjacencies <- lapply(networks, networkAdjacency)
    union_adj <- Reduce(`+`, adjacencies) > 0
    p <- nrow(nodes)
    rows <- list()
    for (i in seq_len(p - 1L)) {
        for (j in (i + 1L):p) {
            if (!union_adj[i, j]) next
            values <- vapply(partial, function(m) m[i, j], numeric(1))
            present <- vapply(adjacencies, function(m) m[i, j] != 0,
                logical(1))
            category <- if (all(present)) {
                "shared"
            } else if (sum(present) == 1L) {
                paste0("unique_", celltypes[which(present)])
            } else {
                "differential"
            }
            rows[[length(rows) + 1L]] <- data.frame(
                node_i = nodes$feature[i], node_j = nodes$feature[j],
                layer_i = nodes$omic_layer[i], layer_j = nodes$omic_layer[j],
                cross_omic = nodes$omic_layer[i] != nodes$omic_layer[j],
                category = category, max_diff = max(values) - min(values),
                stringsAsFactors = FALSE)
        }
    }
    S4Vectors::DataFrame(do.call(rbind, rows))
}

.random_network <- function(seed, p, density, celltype, n_transcripts) {
    set.seed(seed)
    ## paste0() of an empty sequence gives one string, hence the subset
    features <- c(paste0("g", seq_len(n_transcripts)),
        paste0("m", seq_len(p - n_transcripts)))[seq_len(p)]
    upper <- upper.tri(matrix(0, p, p))
    adjacency <- matrix(0L, p, p)
    adjacency[upper] <- stats::rbinom(sum(upper), 1L, density)
    adjacency <- adjacency + t(adjacency)
    weight <- matrix(0, p, p)
    weight[upper] <- stats::runif(sum(upper), -0.3, 0.3)
    precision <- adjacency * (weight + t(weight))
    diag(precision) <- stats::runif(p, 1, 2)
    stability <- matrix(0, p, p)
    stability[upper] <- stats::runif(sum(upper))
    stability <- stability + t(stability)
    ## the entries below the diagonal are halved, so that a table read from
    ## [j, i] in place of [i, j] would differ from that of the loops
    lower <- lower.tri(precision)
    precision[lower] <- precision[lower] / 2
    stability[lower] <- stability[lower] / 2
    dimnames(precision) <- dimnames(adjacency) <- dimnames(stability) <-
        list(features, features)
    CelltypeNetworkResult(precision = precision, adjacency = adjacency,
        node_info = S4Vectors::DataFrame(feature = features,
            omic_layer = rep(c("transcript", "metabolite"),
                c(n_transcripts, p - n_transcripts))),
        celltype = celltype, stability = stability)
}

test_that("edges are listed by row and then by column", {
    upper_edges <- exposomeSC:::.upper_edges
    ## the smallest case in which the two orders differ: which() alone
    ## lists (2, 3) before (1, 4)
    adjacency <- matrix(0L, 4, 4)
    adjacency[1, 4] <- adjacency[4, 1] <- 1L
    adjacency[2, 3] <- adjacency[3, 2] <- 1L
    expect_identical(upper_edges(adjacency),
        matrix(c(1L, 2L, 4L, 3L), 2L))
    by_column <- which(upper.tri(adjacency) & adjacency != 0,
        arr.ind = TRUE, useNames = FALSE)
    expect_identical(by_column, matrix(c(2L, 1L, 3L, 4L), 2L))
    ## only entries above the diagonal are edges, whatever their value
    odd <- matrix(0, 3, 3)
    odd[1, 3] <- -2
    odd[3, 2] <- 1
    diag(odd) <- 5
    expect_identical(upper_edges(odd), matrix(c(1L, 3L), 1L))
    expect_identical(nrow(upper_edges(matrix(0, 3, 3))), 0L)
    expect_identical(dim(upper_edges(matrix(0, 1, 1))), c(0L, 2L))
    expect_identical(dim(upper_edges(matrix(0, 0, 0))), c(0L, 2L))
    ## a missing value is an error above the diagonal and unread below it
    odd[2, 1] <- NA
    expect_identical(upper_edges(odd), matrix(c(1L, 3L), 1L))
    odd[1, 2] <- NA
    expect_error(upper_edges(odd),
        "The adjacency matrix has a missing value above its diagonal.",
        fixed = TRUE)
})

test_that("edge counts and edge tables equal those of the loops", {
    count <- exposomeSC:::.count_edge_types
    extract <- exposomeSC:::.extract_edges
    for (case in list(c(p = 4, density = 0.5, tx = 2),
        c(p = 12, density = 0.2, tx = 7), c(p = 12, density = 0.9, tx = 5),
        c(p = 12, density = 0, tx = 6), c(p = 9, density = 0.4, tx = 9),
        c(p = 2, density = 1, tx = 1), c(p = 1, density = 0, tx = 1))) {
        network <- .random_network(2026L, case[["p"]], case[["density"]],
            "A", case[["tx"]])
        adjacency <- networkAdjacency(network)
        nodes <- networkNodes(network)
        expect_identical(count(adjacency, nodes),
            .loop_count_edge_types(adjacency, nodes))
        expect_identical(extract(adjacency, nodes),
            .loop_extract_edges(adjacency, nodes))
        expect_identical(
            extract(adjacency, nodes, networkPrecision(network),
                networkStability(network)),
            .loop_extract_edges(adjacency, nodes, networkPrecision(network),
                networkStability(network)))
    }
    ## the order in a table: by the first node, then by the second
    adjacency <- matrix(0L, 4, 4)
    adjacency[1, 4] <- adjacency[4, 1] <- 1L
    adjacency[2, 3] <- adjacency[3, 2] <- 1L
    nodes <- S4Vectors::DataFrame(feature = c("g1", "g2", "m1", "m2"),
        omic_layer = rep(c("transcript", "metabolite"), each = 2L))
    edges <- extract(adjacency, nodes, precision = matrix(1:16 / 10, 4, 4))
    expect_identical(edges$node_i, c("g1", "g2"))
    expect_identical(edges$node_j, c("m2", "m1"))
    expect_identical(edges$weight, c(1.3, 1))
    expect_identical(edges$stability, c(NA_real_, NA_real_))
    expect_identical(edges$cross_omic, c(TRUE, TRUE))
    expect_null(rownames(edges))
    ## the stability scores are those above the diagonal as well
    scored <- extract(adjacency, nodes, stability = matrix(1:16 / 100, 4, 4))
    expect_identical(scored$stability, c(13, 10) / 100)
    expect_identical(count(adjacency, nodes), list(total = 2L,
        within_transcript = 0L, within_metabolite = 0L, cross_omic = 2L))
    ## a same-layer edge counts as within_metabolite unless its layer is
    ## "transcript", and a node with an edge needs a layer
    nodes$omic_layer <- c("transcript", "protein", "protein", "transcript")
    expect_identical(count(adjacency, nodes), list(total = 2L,
        within_transcript = 1L, within_metabolite = 1L, cross_omic = 0L))
    nodes$omic_layer[2L] <- NA
    expect_error(count(adjacency, nodes),
        "A node with an edge has no omic layer.", fixed = TRUE)
    expect_error(count(adjacency, nodes[, "feature", drop = FALSE]),
        "A node with an edge has no omic layer.", fixed = TRUE)
    ## no loop over pairs of nodes is left
    expect_false("for" %in% all.names(body(count)))
    expect_false("for" %in% all.names(body(extract)))
})

test_that("the comparison lists its edges as the loop did", {
    for (density in c(0.15, 0.6)) {
        networks <- lapply(1:3, function(k) {
            .random_network(100L + k, 10L, density, paste0("ct", k), 6L)
        })
        comparison <- do.call(run_comparative_network, networks)
        edges <- .loop_comparison_edges(networks)
        shared <- edges$category == "shared"
        expect_identical(sharedEdges(comparison), edges[shared, ])
        expect_identical(diffEdges(comparison), edges[!shared, ])
        for (celltype in names(celltypeSpecificEdges(comparison))) {
            expect_identical(celltypeSpecificEdges(comparison)[[celltype]],
                edges[edges$category == paste0("unique_", celltype), ])
        }
        expect_identical(comparisonSummary(comparison)$n_total_edges,
            nrow(edges))
    }
    ## all three categories occur in the denser comparison
    networks <- lapply(1:3, function(k) {
        .random_network(100L + k, 10L, 0.6, paste0("ct", k), 6L)
    })
    categories <- .loop_comparison_edges(networks)$category
    expect_true(all(c("shared", "differential", "unique_ct1") %in%
        categories))
    ## the comparison finds its edges with the function the summaries use
    expect_true(".upper_edges" %in%
        all.names(body(run_comparative_network)))
})

test_that("names on the columns of a node table do not reach the edge tables", {
    features <- c("g1", "g2", "m1", "m2")
    layer_of <- c(g1 = "transcript", g2 = "transcript", m1 = "metabolite",
        m2 = "metabolite")
    make <- function(celltype, edges, layers, feature = features) {
        precision <- diag(4)
        adjacency <- matrix(0L, 4, 4)
        for (edge in edges) {
            precision[edge[1], edge[2]] <- precision[edge[2], edge[1]] <- -0.3
            adjacency[edge[1], edge[2]] <- adjacency[edge[2], edge[1]] <- 1L
        }
        dimnames(precision) <- dimnames(adjacency) <- list(features, features)
        CelltypeNetworkResult(precision, adjacency,
            node_info = S4Vectors::DataFrame(feature = feature,
                omic_layer = layers),
            celltype = celltype)
    }
    compare <- function(layers, feature = features) {
        run_comparative_network(
            make("A", list(c(1, 2), c(1, 3), c(1, 4)), layers, feature),
            make("B", list(c(1, 2), c(1, 3)), layers, feature))
    }
    plain <- compare(unname(layer_of))
    expect_null(rownames(sharedEdges(plain)))
    expect_null(rownames(diffEdges(plain)))
    ## a lookup by feature, layer_of[features], leaves the names of the
    ## features on the column; a feature it does not know gets the name NA
    unknown <- stats::setNames(layer_of, c("g1", "g2", "m1", NA))
    expect_identical(compare(layer_of[features]), plain)
    expect_identical(compare(unknown), plain)
    expect_identical(
        compare(unname(layer_of), stats::setNames(features, features)), plain)
    network <- make("A", list(c(1, 2), c(1, 3), c(1, 4)), unknown)
    edges <- exposomeSC:::.extract_edges(networkAdjacency(network),
        networkNodes(network), networkPrecision(network))
    expect_null(rownames(edges))
    expect_identical(edges$node_j, c("g2", "m1", "m2"))
    expect_identical(edges$layer_j,
        c("transcript", "metabolite", "metabolite"))
})

test_that("networks without nodes are compared and plotted", {
    empty <- function(celltype) {
        CelltypeNetworkResult(matrix(numeric(), 0L, 0L),
            matrix(integer(), 0L, 0L),
            node_info = S4Vectors::DataFrame(feature = character(),
                omic_layer = character()),
            celltype = celltype)
    }
    comparison <- run_comparative_network(empty("A"), empty("B"))
    expect_s4_class(comparison, "NetworkComparison")
    expect_identical(nrow(sharedEdges(comparison)), 0L)
    expect_identical(nrow(diffEdges(comparison)), 0L)
    expect_identical(comparisonSummary(comparison)$celltypes, c("A", "B"))
    expect_identical(colnames(sharedEdges(comparison)),
        c("node_i", "node_j", "layer_i", "layer_j", "cross_omic", "category",
            "max_diff"))
    skip_if_not_installed("ggplot2")
    skip_if_not_installed("ggraph")
    skip_if_not_installed("igraph")
    skip_if_not_installed("tidygraph")
    expect_message(plot <- plot_celltype_network(empty("A")),
        "No edges to plot", fixed = TRUE)
    expect_identical(plot$labels$title, "A: no edges detected")
})

test_that("show() counts the cross-omic edges without a loop", {
    adjacency <- matrix(0L, 4, 4)
    adjacency[1, 4] <- adjacency[4, 1] <- 1L
    adjacency[2, 3] <- adjacency[3, 2] <- 1L
    adjacency[1, 2] <- adjacency[2, 1] <- 1L
    features <- c("g1", "g2", "m1", "m2")
    dimnames(adjacency) <- list(features, features)
    precision <- adjacency + diag(4)
    network <- CelltypeNetworkResult(precision = precision,
        adjacency = adjacency,
        node_info = S4Vectors::DataFrame(feature = features,
            omic_layer = rep(c("transcript", "metabolite"), each = 2L)),
        celltype = "Mono", metadata = list(n_donors = 12L))
    expect_identical(utils::capture.output(show(network)), c(
        "CelltypeNetworkResult for 'Mono'",
        "  Method: precomputed",
        "  Nodes: 4 (2 transcripts, 2 metabolites)",
        "  Edges: 3 (2 cross-omic)",
        "  Donors: 12"))
    show_method <- methods::getMethod("show", "CelltypeNetworkResult")
    expect_false("for" %in% all.names(body(show_method)))
})
