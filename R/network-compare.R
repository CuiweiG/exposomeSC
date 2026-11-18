# R/network-compare.R
# Compare networks across cell types or conditions

#' @include AllClasses.R
#' @include network-utils.R
#' @importFrom S4Vectors DataFrame
NULL

#' Compare networks across cell types
#'
#' Classifies the union of edges from two or more cell-type-specific
#' networks estimated on the same features as shared, unique to one
#' cell type or differential, and reports pairwise Jaccard overlap.
#' The comparison is descriptive: no edge is tested and no p-value is
#' reported.
#'
#' @param ... Two or more \code{\linkS4class{CelltypeNetworkResult}}
#'   objects. All must share the same feature set (same genes
#'   and metabolites).
#'
#' @return A \code{\linkS4class{NetworkComparison}} object. Its edge
#'   tables contain \code{max_diff}, the range of the edge's partial
#'   correlation across networks.
#'
#' @details
#' Partial correlations are computed from each precision matrix
#' \eqn{\Theta} as \eqn{-\Theta_{ij} / \sqrt{\Theta_{ii}\Theta_{jj}}}.
#'
#' This is the single-cell analogue of the temporal edge comparison in
#' Cheng et al. (ES&T 2024).
#'
#' @references
#' Cheng SL et al. (2024). Multiomic signatures of traffic-related
#'   air pollution in London reveal potential short-term perturbations
#'   in gut microbiome-related pathways. \emph{Environ Sci Technol}
#'   58:8771-8782. \doi{10.1021/acs.est.3c09148}
#'
#' @export
#' @examples
#' features <- c("G1", "G2", "M1")
#' make_network <- function(celltype, cross_weight) {
#'     precision <- diag(3)
#'     precision[1, 3] <- precision[3, 1] <- cross_weight
#'     adjacency <- (precision != 0) * 1L
#'     diag(adjacency) <- 0L
#'     dimnames(precision) <- dimnames(adjacency) <- list(features, features)
#'     CelltypeNetworkResult(
#'         precision = precision,
#'         adjacency = adjacency,
#'         node_info = S4Vectors::DataFrame(
#'             feature = features,
#'             omic_layer = c("transcript", "transcript", "metabolite")
#'         ),
#'         celltype = celltype,
#'         metadata = list(n_donors = 20L)
#'     )
#' }
#' monocyte_network <- make_network("Monocyte", -0.45)
#' b_cell_network <- make_network("B cell", -0.15)
#' comparison <- run_comparative_network(monocyte_network, b_cell_network)
#' comparison
run_comparative_network <- function(...) {
    nets <- list(...)

    if (length(nets) < 2L)
        stop("At least 2 CelltypeNetworkResult objects required")

    ## Validate all are CelltypeNetworkResult; an argument given by name is
    ## named in the message, as every argument is taken to be a network
    arg_names <- names(nets)
    for (i in seq_along(nets)) {
        if (!is(nets[[i]], "CelltypeNetworkResult")) {
            given_as <- if (!is.null(arg_names) && nzchar(arg_names[[i]])) {
                paste0(" ('", arg_names[[i]], "')")
            }
            stop("Argument ", i, given_as,
                 " is not a CelltypeNetworkResult")
        }
    }

    ## Get cell type names
    ct_names <- vapply(nets, networkCelltype, character(1))
    if (any(duplicated(ct_names)))
        stop("Duplicate cell types found: ",
             paste(ct_names[duplicated(ct_names)],
                   collapse = ", "))

    ## Check feature alignment
    ref_nodes <- networkNodes(nets[[1]])
    ref_features <- ref_nodes$feature
    p <- length(ref_features)

    for (i in seq_along(nets)[-1]) {
        fi <- networkNodes(nets[[i]])$feature
        if (length(fi) != p || !all(fi == ref_features)) {
            only_first <- setdiff(ref_features, fi)
            only_this <- setdiff(fi, ref_features)
            detail <- c(
                if (length(only_first)) paste0(" Only in '", ct_names[1],
                    "': ", .some_names(only_first), "."),
                if (length(only_this)) paste0(" Only in '", ct_names[i],
                    "': ", .some_names(only_this), "."),
                if (!length(only_first) && !length(only_this) &&
                    length(fi) == p)
                    " The features are the same, in another order.")
            stop("All networks must have identical features. ",
                 "Network for '", ct_names[i],
                 "' has different features from '",
                 ct_names[1], "'.", detail)
        }
    }

    ## --- Partial correlations from the precision matrices ---
    precisions <- lapply(nets, function(n) {
        theta <- networkPrecision(n)
        root_diag <- sqrt(diag(theta))
        partial <- -theta / outer(root_diag, root_diag)
        diag(partial) <- 1
        partial
    })
    adjacencies <- lapply(nets, networkAdjacency)

    ## --- Identify all edges across networks ---
    union_adj <- Reduce(`+`, adjacencies)
    union_adj <- (union_adj > 0) * 1L

    ## --- Classify each edge ---
    ij <- .upper_edges(union_adj)
    n_union <- nrow(ij)

    if (n_union == 0) {
        all_edges <- S4Vectors::DataFrame(
            node_i = character(), node_j = character(),
            layer_i = character(), layer_j = character(),
            cross_omic = logical(), category = character(),
            max_diff = numeric())

        return(new("NetworkComparison",
            diff_edges        = all_edges,
            shared_edges      = all_edges,
            celltype_specific = list(),
            summary = list(
                celltypes = ct_names,
                n_edges = integer(0),
                jaccard = 0)))
    }

    ## One row per edge and one column per network: the partial correlation
    ## of the edge and whether the network has it
    pcor <- matrix(
        vapply(precisions, function(pm) pm[ij], numeric(n_union)),
        nrow = n_union
    )
    present <- matrix(
        vapply(adjacencies, function(am) am[ij] != 0, logical(n_union)),
        nrow = n_union
    )
    n_present <- rowSums(present)
    category <- rep("differential", n_union)
    category[n_present == length(nets)] <- "shared"
    ## for an edge found in one network, the sum of the positions of the
    ## networks that have it is the position of that network
    in_one <- n_present == 1L
    category[in_one] <- paste0(
        "unique_",
        ct_names[drop(present %*% seq_along(nets))[in_one]]
    )

    ## row.names = NULL: data.frame() would otherwise take the row names
    ## from names that a column of the node table may carry
    layers <- ref_nodes$omic_layer
    all_df <- data.frame(
        node_i     = ref_features[ij[, 1L]],
        node_j     = ref_features[ij[, 2L]],
        layer_i    = layers[ij[, 1L]],
        layer_j    = layers[ij[, 2L]],
        cross_omic = layers[ij[, 1L]] != layers[ij[, 2L]],
        category   = category,
        max_diff   = apply(pcor, 1L, max) - apply(pcor, 1L, min),
        row.names = NULL,
        stringsAsFactors = FALSE
    )
    all_edges <- S4Vectors::DataFrame(all_df)

    ## Partition edges
    shared <- all_edges[all_edges$category == "shared", ]
    diff <- all_edges[grepl("unique_|differential",
                             all_edges$category), ]

    ct_specific <- list()
    for (ct in ct_names) {
        pattern <- paste0("unique_", ct)
        ct_rows <- all_edges$category == pattern
        if (any(ct_rows))
            ct_specific[[ct]] <- all_edges[ct_rows, ]
    }

    ## Jaccard similarity (pairwise)
    n_edges_per <- vapply(adjacencies, function(am) {
        sum(am[upper.tri(am)] != 0)
    }, integer(1))
    names(n_edges_per) <- ct_names

    jaccard_pairs <- list()
    for (a in seq_along(nets)[-length(nets)]) {
        for (b in (a + 1):length(nets)) {
            am_a <- adjacencies[[a]][upper.tri(adjacencies[[a]])]
            am_b <- adjacencies[[b]][upper.tri(adjacencies[[b]])]
            inter <- sum(am_a != 0 & am_b != 0)
            union <- sum(am_a != 0 | am_b != 0)
            j <- if (union > 0) inter / union else 0
            jaccard_pairs[[paste(ct_names[a], ct_names[b],
                                  sep = "_vs_")]] <- j
        }
    }

    new("NetworkComparison",
        diff_edges        = diff,
        shared_edges      = shared,
        celltype_specific = ct_specific,
        summary = list(
            celltypes     = ct_names,
            n_edges       = n_edges_per,
            jaccard       = jaccard_pairs,
            n_total_edges = nrow(all_edges),
            n_shared      = nrow(shared),
            n_differential = nrow(diff)
        )
    )
}
