# R/network-plots.R
# Visualisation for cross-omic network inference

#' @include AllClasses.R
#' @include network-utils.R
#' @importFrom S4Vectors DataFrame
NULL

# Suppress R CMD check NOTEs for ggplot2 .data pronoun
utils::globalVariables(".data")

#' Plot cell-type-specific cross-omic network
#'
#' Renders a cell-type network graph with nodes coloured by
#' omic layer and edges weighted by partial correlation
#' strength or stability.
#'
#' @param network \code{\linkS4class{CelltypeNetworkResult}}.
#' @param layout Character; graph layout algorithm.
#'   \code{"fr"} (Fruchterman-Reingold, default),
#'   \code{"circle"}, \code{"grid"}, or \code{"bipartite"}.
#' @param color_by Character; node coloring scheme.
#'   \code{"omic_layer"} (default), \code{"stability"},
#'   or \code{"community"}.
#' @param highlight_cross_omic Logical; emphasize
#'   transcript-metabolite edges. Default TRUE.
#' @param label_top Integer; label top N nodes by degree.
#'   Default 10.
#' @param edge_alpha_by Character; \code{"weight"} (default)
#'   or \code{"stability"}.
#'
#' @return A \code{ggplot} object (requires \pkg{ggraph}
#'   and \pkg{igraph}).
#'
#' @export
#' @examples
#' features <- c("G1", "G2", "M1")
#' precision <- diag(3)
#' precision[1, 3] <- precision[3, 1] <- -0.4
#' adjacency <- (precision != 0) * 1L
#' diag(adjacency) <- 0L
#' dimnames(precision) <- dimnames(adjacency) <- list(features, features)
#' network <- methods::new(
#'     "CelltypeNetworkResult",
#'     precision_matrix = precision,
#'     adjacency_matrix = adjacency,
#'     stability_scores = matrix(nrow = 0, ncol = 0),
#'     node_info = S4Vectors::DataFrame(
#'         feature = features,
#'         omic_layer = c("transcript", "transcript", "metabolite")
#'     ),
#'     celltype = "Monocyte",
#'     method = "precomputed",
#'     metadata = list(n_donors = 20L)
#' )
#' plot_packages <- c("igraph", "ggraph", "ggplot2", "tidygraph")
#' if (all(vapply(plot_packages, requireNamespace, logical(1), quietly = TRUE))) {
#'     network_plot <- plot_celltype_network(
#'         network,
#'         layout = "circle",
#'         label_top = 3L
#'     )
#'     inherits(network_plot, "ggplot")
#' }
plot_celltype_network <- function(network,
                                   layout = c("fr", "circle",
                                              "grid",
                                              "bipartite"),
                                   color_by = c("omic_layer",
                                                "stability",
                                                "community"),
                                   highlight_cross_omic = TRUE,
                                   label_top = 10L,
                                   edge_alpha_by = c("weight",
                                                     "stability")) {

    layout <- match.arg(layout)
    color_by <- match.arg(color_by)
    edge_alpha_by <- match.arg(edge_alpha_by)

    if (!requireNamespace("igraph", quietly = TRUE) ||
        !requireNamespace("ggraph", quietly = TRUE) ||
        !requireNamespace("ggplot2", quietly = TRUE))
        stop("Packages 'igraph', 'ggraph', and 'ggplot2' ",
             "required for plot_celltype_network()")

    if (!requireNamespace("tidygraph", quietly = TRUE))
        stop("Package 'tidygraph' required for ",
             "plot_celltype_network()")

    stopifnot(is(network, "CelltypeNetworkResult"))

    ## Build edge list
    edges <- .extract_edges(
        network@adjacency_matrix,
        network@node_info,
        precision = network@precision_matrix,
        stability = if (nrow(network@stability_scores) > 0)
            network@stability_scores else NULL)

    if (nrow(edges) == 0) {
        message("[exposomeSC] No edges to plot.")
        return(ggplot2::ggplot() +
            ggplot2::theme_void() +
            ggplot2::ggtitle(paste0(
                network@celltype, ": no edges detected")))
    }

    ## Build igraph
    g <- igraph::graph_from_data_frame(
        d = data.frame(
            from = as.character(edges$node_i),
            to = as.character(edges$node_j),
            weight = abs(as.numeric(edges$weight)),
            stability = as.numeric(edges$stability),
            cross_omic = as.logical(edges$cross_omic),
            stringsAsFactors = FALSE),
        directed = FALSE,
        vertices = data.frame(
            name = as.character(network@node_info$feature),
            omic_layer = as.character(
                network@node_info$omic_layer),
            stringsAsFactors = FALSE)
    )

    ## Node degree
    deg <- igraph::degree(g)
    igraph::V(g)$degree <- deg

    ## Community detection for color_by = "community"
    if (color_by == "community") {
        comm <- igraph::cluster_louvain(g)
        igraph::V(g)$community <- as.character(
            igraph::membership(comm))
    }

    ## Labels: top N by degree
    top_nodes <- names(sort(deg, decreasing = TRUE))[
        seq_len(min(label_top, length(deg)))]
    igraph::V(g)$label <- ifelse(
        igraph::V(g)$name %in% top_nodes,
        igraph::V(g)$name, "")

    ## Convert to tidygraph
    tg <- tidygraph::as_tbl_graph(g)

    ## Build ggraph plot
    layout_algo <- switch(layout,
        fr        = "fr",
        circle    = "circle",
        grid      = "grid",
        bipartite = "bipartite"
    )

    ## Handle bipartite layout
    if (layout == "bipartite") {
        igraph::V(g)$type <- igraph::V(g)$omic_layer ==
            "metabolite"
        tg <- tidygraph::as_tbl_graph(g)
    }

    p <- ggraph::ggraph(tg, layout = layout_algo)

    ## Edge aesthetics
    if (highlight_cross_omic) {
        p <- p + ggraph::geom_edge_link(
            ggplot2::aes(
                alpha = if (edge_alpha_by == "stability" &&
                            any(!is.na(
                                igraph::E(g)$stability)))
                    .data$stability else .data$weight,
                color = .data$cross_omic),
            show.legend = TRUE) +
            ggraph::scale_edge_color_manual(
                values = c("FALSE" = "grey60",
                           "TRUE" = "#E74C3C"),
                labels = c("Within-omic", "Cross-omic"),
                name = "Edge type")
    } else {
        p <- p + ggraph::geom_edge_link(
            ggplot2::aes(
                alpha = .data$weight),
            color = "grey50",
            show.legend = FALSE)
    }

    ## Node aesthetics
    fill_aes <- switch(color_by,
        omic_layer = ggplot2::aes(fill = .data$omic_layer),
        community  = ggplot2::aes(fill = .data$community),
        stability  = ggplot2::aes(fill = .data$omic_layer)
    )

    p <- p +
        ggraph::geom_node_point(
            fill_aes,
            shape = 21, size = 4, color = "white") +
        ggraph::geom_node_text(
            ggplot2::aes(label = .data$label),
            repel = TRUE, size = 3)

    ## Colour scales
    if (color_by == "omic_layer") {
        p <- p + ggplot2::scale_fill_manual(
            values = c(transcript = "#3498DB",
                       metabolite = "#E67E22"),
            name = "Omic layer")
    }

    ## Theme
    ec <- network@metadata$edge_counts
    subtitle <- sprintf(
        "%d edges (%d cross-omic) | %d donors",
        if (!is.null(ec)) ec$total else sum(
            network@adjacency_matrix[
                upper.tri(network@adjacency_matrix)] != 0),
        if (!is.null(ec)) ec$cross_omic else 0L,
        if (!is.null(network@metadata$n_donors))
            network@metadata$n_donors else 0L)

    p <- p +
        ggraph::theme_graph() +
        ggplot2::ggtitle(
            paste0("Cross-omic network: ", network@celltype),
            subtitle = subtitle)

    p
}


#' Side-by-side network comparison across cell types
#'
#' Visualises differential network edges between cell types.
#'
#' @param comparison \code{\linkS4class{NetworkComparison}}.
#' @param highlight Character; \code{"differential"} (default),
#'   \code{"shared"}, or \code{"all"}.
#'
#' @return A \code{ggplot} object.
#'
#' @export
#' @examples
#' shared_edge <- S4Vectors::DataFrame(
#'     node_i = "G1", node_j = "M1", category = "shared",
#'     max_diff = 0.1, p_value = 0.4, p_adjusted = 0.4
#' )
#' differential_edge <- S4Vectors::DataFrame(
#'     node_i = "G2", node_j = "M1", category = "differential",
#'     max_diff = 0.5, p_value = 0.01, p_adjusted = 0.02
#' )
#' comparison <- methods::new(
#'     "NetworkComparison",
#'     diff_edges = differential_edge,
#'     shared_edges = shared_edge,
#'     celltype_specific = list(),
#'     summary = list(
#'         celltypes = c("Monocyte", "B cell"),
#'         n_edges = c(Monocyte = 2L, `B cell` = 1L),
#'         n_shared = 1L,
#'         n_differential = 1L,
#'         jaccard = list(Monocyte_vs_B_cell = 0.5)
#'     )
#' )
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'     comparison_plot <- plot_network_comparison(
#'         comparison,
#'         highlight = "differential"
#'     )
#'     inherits(comparison_plot, "ggplot")
#' }
plot_network_comparison <- function(comparison,
                                     highlight = c(
                                         "differential",
                                         "shared", "all")) {
    highlight <- match.arg(highlight)

    if (!requireNamespace("ggplot2", quietly = TRUE))
        stop("Package 'ggplot2' required")

    stopifnot(is(comparison, "NetworkComparison"))

    ## Build edge-by-celltype heatmap
    edge_tables <- list(
        comparison@shared_edges,
        comparison@diff_edges
    )
    non_empty <- vapply(edge_tables, nrow, integer(1)) > 0L
    edge_tables <- edge_tables[non_empty]
    all_edges <- if (!length(edge_tables)) {
        NULL
    } else if (length(edge_tables) == 1L) {
        edge_tables[[1L]]
    } else {
        do.call(rbind, edge_tables)
    }

    if (is.null(all_edges) || nrow(all_edges) == 0) {
        return(ggplot2::ggplot() +
            ggplot2::theme_void() +
            ggplot2::ggtitle("No edges to compare"))
    }

    ## Summary bar chart: edges per cell type
    ct_names <- comparison@summary$celltypes
    n_edges <- comparison@summary$n_edges

    if (length(n_edges) > 0) {
        df <- data.frame(
            celltype = names(n_edges),
            n_edges = as.integer(n_edges),
            stringsAsFactors = FALSE)

        p <- ggplot2::ggplot(df, ggplot2::aes(
            x = stats::reorder(.data$celltype, -.data$n_edges),
            y = .data$n_edges)) +
            ggplot2::geom_col(
                fill = "#3498DB", alpha = 0.8) +
            ggplot2::geom_text(
                ggplot2::aes(label = .data$n_edges),
                vjust = -0.5, size = 4) +
            ggplot2::labs(
                title = "Network comparison across cell types",
                subtitle = sprintf(
                    "%d shared, %d differential edges",
                    comparison@summary$n_shared,
                    comparison@summary$n_differential),
                x = "Cell type",
                y = "Number of edges") +
            ggplot2::theme_minimal(base_size = 14)

        ## Add Jaccard annotations
        jac <- comparison@summary$jaccard
        if (length(jac) > 0) {
            jac_text <- paste(
                vapply(names(jac), function(nm) {
                    sprintf("%s: J=%.2f", nm, jac[[nm]])
                }, character(1)),
                collapse = "\n")
            p <- p + ggplot2::annotate(
                "text", x = Inf, y = Inf,
                label = jac_text,
                hjust = 1.1, vjust = 1.1,
                size = 3, color = "grey40")
        }

        return(p)
    }

    ggplot2::ggplot() + ggplot2::theme_void()
}


#' Stability calibration surface
#'
#' Displays the lambda-pi calibration heatmap for
#' hyperparameter transparency, reproducing the StARS
#' stability surface.
#'
#' @param network \code{\linkS4class{CelltypeNetworkResult}}
#'   with \code{stability=TRUE}.
#'
#' @return A \code{ggplot} object.
#'
#' @export
#' @examples
#' features <- c("G1", "G2", "M1")
#' stability <- matrix(0, 3, 3, dimnames = list(features, features))
#' stability[1, 3] <- stability[3, 1] <- 0.85
#' precision <- diag(3)
#' dimnames(precision) <- list(features, features)
#' network <- methods::new(
#'     "CelltypeNetworkResult",
#'     precision_matrix = precision,
#'     adjacency_matrix = (stability > 0) * 1L,
#'     stability_scores = stability,
#'     node_info = S4Vectors::DataFrame(
#'         feature = features,
#'         omic_layer = c("transcript", "transcript", "metabolite")
#'     ),
#'     celltype = "Monocyte",
#'     method = "precomputed",
#'     metadata = list(n_donors = 20L)
#' )
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'     stability_plot <- plot_stability_surface(network)
#'     inherits(stability_plot, "ggplot")
#' }
plot_stability_surface <- function(network) {

    if (!requireNamespace("ggplot2", quietly = TRUE))
        stop("Package 'ggplot2' required")

    stopifnot(is(network, "CelltypeNetworkResult"))

    stab <- network@stability_scores
    if (nrow(stab) == 0) {
        message("[exposomeSC] No stability scores available. ",
                "Re-run with stability=TRUE.")
        return(ggplot2::ggplot() +
            ggplot2::theme_void() +
            ggplot2::ggtitle("No stability data available"))
    }

    ## Extract upper triangle stability scores
    ut <- upper.tri(stab)
    stab_vals <- stab[ut]
    ni <- network@node_info

    ## Build heatmap data
    pairs <- which(ut, arr.ind = TRUE)
    df <- data.frame(
        feature_i = ni$feature[pairs[, 1]],
        feature_j = ni$feature[pairs[, 2]],
        stability = stab_vals,
        cross_omic = ni$omic_layer[pairs[, 1]] !=
            ni$omic_layer[pairs[, 2]],
        stringsAsFactors = FALSE
    )

    ## Filter to edges with non-zero stability
    df <- df[df$stability > 0, ]

    if (nrow(df) == 0) {
        return(ggplot2::ggplot() +
            ggplot2::theme_void() +
            ggplot2::ggtitle("No stable edges found"))
    }

    ## Sort by stability
    df <- df[order(-df$stability), ]

    ## Plot top edges as bar chart with stability threshold
    top_n <- min(50L, nrow(df))
    df_top <- df[seq_len(top_n), ]
    df_top$edge_label <- paste0(df_top$feature_i, " - ",
                                 df_top$feature_j)
    df_top$edge_label <- factor(df_top$edge_label,
        levels = rev(df_top$edge_label))

    ggplot2::ggplot(df_top, ggplot2::aes(
        x = .data$stability,
        y = .data$edge_label,
        fill = .data$cross_omic)) +
        ggplot2::geom_col(alpha = 0.8) +
        ggplot2::geom_vline(
            xintercept = 0.6,
            linetype = "dashed", color = "red") +
        ggplot2::scale_fill_manual(
            values = c("FALSE" = "#3498DB",
                       "TRUE" = "#E74C3C"),
            labels = c("Within-omic", "Cross-omic"),
            name = "Edge type") +
        ggplot2::labs(
            title = paste0("Edge stability: ",
                           network@celltype),
            subtitle = "Dashed line: pi=0.6 threshold",
            x = "Stability selection probability",
            y = NULL) +
        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::theme(
            axis.text.y = ggplot2::element_text(size = 7))
}


#' Temporal network dynamics panel
#'
#' Tracks edge presence across time points, showing
#' emerging, disappearing, and persistent edges.
#'
#' @param temporal \code{\linkS4class{TemporalNetwork}}.
#' @param track_node Character or NULL; highlight a specific
#'   node's edges across time. Default NULL.
#'
#' @return A \code{ggplot} object.
#'
#' @export
#' @examples
#' edge_dynamics <- S4Vectors::DataFrame(
#'     node_i = c("G1", "G2"),
#'     node_j = c("M1", "M1"),
#'     cross_omic = c(TRUE, TRUE),
#'     status = c("persistent", "emerging"),
#'     present_baseline = c(TRUE, FALSE),
#'     present_follow_up = c(TRUE, TRUE)
#' )
#' temporal_network <- methods::new(
#'     "TemporalNetwork",
#'     networks = list(),
#'     edge_dynamics = edge_dynamics,
#'     timepoints = c("baseline", "follow_up")
#' )
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'     temporal_plot <- plot_temporal_dynamics(temporal_network)
#'     inherits(temporal_plot, "ggplot")
#' }
plot_temporal_dynamics <- function(temporal,
                                    track_node = NULL) {

    if (!requireNamespace("ggplot2", quietly = TRUE))
        stop("Package 'ggplot2' required")

    stopifnot(is(temporal, "TemporalNetwork"))

    ed <- temporal@edge_dynamics
    if (nrow(ed) == 0) {
        return(ggplot2::ggplot() +
            ggplot2::theme_void() +
            ggplot2::ggtitle("No edge dynamics to plot"))
    }

    tps <- temporal@timepoints

    ## Filter to specific node if requested
    if (!is.null(track_node)) {
        mask <- ed$node_i == track_node |
            ed$node_j == track_node
        ed <- ed[mask, ]
        if (nrow(ed) == 0) {
            return(ggplot2::ggplot() +
                ggplot2::theme_void() +
                ggplot2::ggtitle(paste0(
                    "No edges for node: ", track_node)))
        }
    }

    ## Build long-form data for plotting
    rows <- list()
    idx <- 0L
    for (r in seq_len(nrow(ed))) {
        edge_label <- paste0(ed$node_i[r], " - ",
                              ed$node_j[r])
        for (tp in tps) {
            col <- paste0("present_", tp)
            present <- if (col %in% colnames(ed))
                as.logical(ed[[col]][r]) else FALSE
            idx <- idx + 1L
            rows[[idx]] <- data.frame(
                edge = edge_label,
                timepoint = tp,
                present = present,
                status = as.character(ed$status[r]),
                cross_omic = as.logical(ed$cross_omic[r]),
                stringsAsFactors = FALSE)
        }
    }

    plot_df <- do.call(rbind, rows)
    plot_df$timepoint <- factor(plot_df$timepoint,
                                 levels = tps)

    ## Sort edges by status for visual grouping
    edge_order <- unique(plot_df$edge[
        order(plot_df$status)])
    plot_df$edge <- factor(plot_df$edge,
                            levels = rev(edge_order))

    status_colors <- c(
        persistent   = "#2ECC71",
        emerging     = "#3498DB",
        disappearing = "#E74C3C",
        transient    = "#F39C12"
    )

    title <- "Temporal edge dynamics"
    if (!is.null(track_node))
        title <- paste0(title, ": ", track_node)

    ggplot2::ggplot(plot_df, ggplot2::aes(
        x = .data$timepoint,
        y = .data$edge,
        fill = .data$present)) +
        ggplot2::geom_tile(color = "white", linewidth = 0.5) +
        ggplot2::scale_fill_manual(
            values = c("TRUE" = "#2C3E50",
                       "FALSE" = "#ECF0F1"),
            labels = c("Absent", "Present"),
            name = "Edge") +
        ggplot2::facet_grid(
            .data$status ~ .,
            scales = "free_y",
            space = "free_y") +
        ggplot2::labs(
            title = title,
            x = "Time point",
            y = NULL) +
        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::theme(
            axis.text.y = ggplot2::element_text(size = 7),
            strip.text.y = ggplot2::element_text(
                angle = 0, hjust = 0))
}
