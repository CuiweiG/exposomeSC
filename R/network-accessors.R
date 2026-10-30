# R/network-accessors.R
# Constructors and accessors for the network result classes

#' @include AllClasses.R
#' @include AllGenerics.R
#' @importFrom S4Vectors DataFrame metadata
NULL

#' Create a CelltypeNetworkResult
#'
#' Wraps a network estimated elsewhere, or built by hand, in a
#' \code{\linkS4class{CelltypeNetworkResult}} so that it can be compared
#' with \code{\link{run_comparative_network}} and plotted. The networks
#' returned by \code{\link{run_celltype_network}} and
#' \code{\link{run_exposure_network}} are already of this class.
#'
#' @param precision Numeric square matrix: the precision (inverse
#'   covariance) matrix, one row and column per feature.
#' @param adjacency Square matrix of edge indicators (1 = edge, 0 = no
#'   edge) with the dimensions of \code{precision}.
#' @param node_info \code{DataFrame} or \code{data.frame} with one row per
#'   feature, in matrix order, and at least a \code{feature} column. The
#'   plotting and comparison functions also use an \code{omic_layer} column
#'   (\code{"transcript"} or \code{"metabolite"}).
#' @param celltype Character scalar: the cell type.
#' @param method Character scalar: how the network was estimated. Default
#'   \code{"precomputed"}.
#' @param stability Optional matrix of edge selection frequencies with the
#'   dimensions of \code{adjacency}; \code{NULL} (default) stores none.
#' @param metadata List of further information, for example
#'   \code{n_donors}, which \code{run_comparative_network()} uses for its
#'   Fisher-Z mode.
#'
#' @return A valid \code{CelltypeNetworkResult}.
#'
#' @seealso \code{\link{network-accessors}}
#' @export
#' @examples
#' features <- c("G1", "G2", "M1")
#' precision <- diag(3)
#' precision[1, 3] <- precision[3, 1] <- -0.4
#' adjacency <- (precision != 0) * 1L
#' diag(adjacency) <- 0L
#' dimnames(precision) <- dimnames(adjacency) <- list(features, features)
#' network <- CelltypeNetworkResult(
#'     precision = precision,
#'     adjacency = adjacency,
#'     node_info = S4Vectors::DataFrame(
#'         feature = features,
#'         omic_layer = c("transcript", "transcript", "metabolite")),
#'     celltype = "Monocyte",
#'     metadata = list(n_donors = 20L))
#' network
CelltypeNetworkResult <- function(precision, adjacency, node_info, celltype,
                                  method = "precomputed", stability = NULL,
                                  metadata = list()) {
    if (!is(node_info, "DataFrame")) {
        node_info <- S4Vectors::DataFrame(node_info, check.names = FALSE)
    }
    if (is.null(stability)) {
        stability <- matrix(numeric(), nrow = 0L, ncol = 0L)
    }
    new("CelltypeNetworkResult",
        precision_matrix = as.matrix(precision),
        adjacency_matrix = as.matrix(adjacency),
        stability_scores = as.matrix(stability),
        node_info        = node_info,
        celltype         = celltype,
        method           = method,
        metadata         = metadata)
}

#' Create a NetworkComparison
#'
#' Builds a \code{\linkS4class{NetworkComparison}} from edge tables. Most
#' users obtain one from \code{\link{run_comparative_network}} instead.
#'
#' @param diff_edges,shared_edges \code{DataFrame}s with one row per edge
#'   and at least the columns \code{node_i} and \code{node_j}.
#'   \code{plot_network_comparison()} also uses \code{category} and
#'   \code{max_diff}. Edges in \code{shared_edges} have category
#'   \code{"shared"}.
#' @param celltype_specific Named list of \code{DataFrame}s, one per cell
#'   type, with the edges found only in that cell type.
#' @param summary List with, for example, \code{celltypes} (the compared
#'   cell types), \code{n_edges} (named edge counts) and \code{jaccard}
#'   (named pairwise overlaps).
#'
#' @return A valid \code{NetworkComparison}.
#'
#' @seealso \code{\link{network-accessors}}
#' @export
#' @examples
#' shared <- S4Vectors::DataFrame(node_i = "G1", node_j = "M1",
#'     category = "shared", max_diff = 0.1)
#' differential <- S4Vectors::DataFrame(node_i = "G2", node_j = "M1",
#'     category = "differential", max_diff = 0.5)
#' NetworkComparison(
#'     diff_edges = differential,
#'     shared_edges = shared,
#'     summary = list(celltypes = c("Monocyte", "B cell")))
NetworkComparison <- function(diff_edges = S4Vectors::DataFrame(),
                              shared_edges = S4Vectors::DataFrame(),
                              celltype_specific = list(),
                              summary = list()) {
    new("NetworkComparison",
        diff_edges        = diff_edges,
        shared_edges      = shared_edges,
        celltype_specific = celltype_specific,
        summary           = summary)
}

#' @rdname network-accessors
#' @export
setMethod("networkPrecision", "CelltypeNetworkResult",
    function(x) x@precision_matrix)

#' @rdname network-accessors
#' @export
setMethod("networkAdjacency", "CelltypeNetworkResult",
    function(x) x@adjacency_matrix)

#' @rdname network-accessors
#' @export
setMethod("networkStability", "CelltypeNetworkResult",
    function(x) x@stability_scores)

#' @rdname network-accessors
#' @export
setMethod("networkNodes", "CelltypeNetworkResult",
    function(x) x@node_info)

#' @rdname network-accessors
#' @export
setMethod("networkCelltype", "CelltypeNetworkResult",
    function(x) x@celltype)

#' @rdname network-accessors
#' @export
setMethod("networkMethod", "CelltypeNetworkResult",
    function(x) x@method)

#' @rdname network-accessors
#' @param ... Ignored by \code{metadata()}.
#' @export
setMethod("metadata", "CelltypeNetworkResult",
    function(x, ...) x@metadata)

#' @rdname network-accessors
#' @export
setMethod("diffEdges", "NetworkComparison",
    function(x) x@diff_edges)

#' @rdname network-accessors
#' @export
setMethod("sharedEdges", "NetworkComparison",
    function(x) x@shared_edges)

#' @rdname network-accessors
#' @export
setMethod("celltypeSpecificEdges", "NetworkComparison",
    function(x) x@celltype_specific)

#' @rdname network-accessors
#' @export
setMethod("comparisonSummary", "NetworkComparison",
    function(x) x@summary)
