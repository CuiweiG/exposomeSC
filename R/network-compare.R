# R/network-compare.R
# Compare networks across cell types or conditions

#' @include AllClasses.R
#' @include network-utils.R
#' @importFrom S4Vectors DataFrame
#' @importFrom stats p.adjust pnorm
NULL

#' Compare networks across cell types
#'
#' Classifies the union of edges from two or more cell-type-specific
#' networks estimated on the same features as shared, unique to one
#' cell type or differential, and reports pairwise Jaccard overlap.
#' The default \code{method = "descriptive"} reports no p-values.
#'
#' \code{method = "fisher_z"} also compares the partial correlations
#' of each edge with Fisher's Z test. This mode is experimental and
#' warns once per session: networks estimated from the same donors are
#' not independent samples, and penalised partial correlations do not
#' follow the sampling distribution the test assumes, so its p-values
#' are not validated.
#'
#' @param ... Two or more \code{\linkS4class{CelltypeNetworkResult}}
#'   objects. All must share the same feature set (same genes
#'   and metabolites).
#' @param method Character; \code{"descriptive"} (default) or
#'   \code{"fisher_z"}.
#' @param adjust Character; p-value adjustment method used with
#'   \code{method = "fisher_z"}. Default \code{"BH"}.
#'
#' @return A \code{\linkS4class{NetworkComparison}} object. Its edge
#'   tables contain \code{max_diff}, the range of the edge's partial
#'   correlation across networks, and \code{p_value} and
#'   \code{p_adjusted}, which are \code{NA} for
#'   \code{method = "descriptive"}.
#'
#' @details
#' Partial correlations are computed from each precision matrix
#' \eqn{\Theta} as \eqn{-\Theta_{ij} / \sqrt{\Theta_{ii}\Theta_{jj}}}.
#' For two networks, Fisher's Z statistic is
#' \code{(atanh(r1) - atanh(r2)) / sqrt(1/(n1 - 3) + 1/(n2 - 3))},
#' where \code{n1} and \code{n2} are the donor counts recorded in the
#' networks' metadata; with more than two networks the largest pairwise
#' statistic is used.
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
#'     methods::new(
#'         "CelltypeNetworkResult",
#'         precision_matrix = precision,
#'         adjacency_matrix = adjacency,
#'         stability_scores = matrix(nrow = 0, ncol = 0),
#'         node_info = S4Vectors::DataFrame(
#'             feature = features,
#'             omic_layer = c("transcript", "transcript", "metabolite")
#'         ),
#'         celltype = celltype,
#'         method = "precomputed",
#'         metadata = list(n_donors = 20L)
#'     )
#' }
#' monocyte_network <- make_network("Monocyte", -0.45)
#' b_cell_network <- make_network("B cell", -0.15)
#' comparison <- run_comparative_network(monocyte_network, b_cell_network)
#' comparison
run_comparative_network <- function(...,
                                     method = c("descriptive",
                                                "fisher_z"),
                                     adjust = "BH") {
    method <- match.arg(method)
    nets <- list(...)

    if (length(nets) < 2L)
        stop("At least 2 CelltypeNetworkResult objects required")

    ## Validate all are CelltypeNetworkResult
    for (i in seq_along(nets)) {
        if (!is(nets[[i]], "CelltypeNetworkResult"))
            stop("Argument ", i,
                 " is not a CelltypeNetworkResult")
    }

    ## Get cell type names
    ct_names <- vapply(nets, function(n) n@celltype, character(1))
    if (any(duplicated(ct_names)))
        stop("Duplicate cell types found: ",
             paste(ct_names[duplicated(ct_names)],
                   collapse = ", "))

    ## Check feature alignment
    ref_features <- nets[[1]]@node_info$feature
    p <- length(ref_features)

    for (i in seq_along(nets)[-1]) {
        fi <- nets[[i]]@node_info$feature
        if (length(fi) != p || !all(fi == ref_features))
            stop("All networks must have identical features. ",
                 "Network for '", ct_names[i],
                 "' has different features from '",
                 ct_names[1], "'.")
    }

    if (method == "fisher_z") {
        .warn_experimental("run_comparative_network_fisher_z", paste0(
            "run_comparative_network(method = \"fisher_z\") is ",
            "experimental: networks from the same donors are not ",
            "independent samples and penalised partial correlations do ",
            "not follow the distribution Fisher's Z test assumes, so its ",
            "p-values are not validated. This warning is shown once per ",
            "session."))
    }

    ## --- Partial correlations from the precision matrices ---
    precisions <- lapply(nets, function(n) {
        theta <- n@precision_matrix
        root_diag <- sqrt(diag(theta))
        partial <- -theta / outer(root_diag, root_diag)
        diag(partial) <- 1
        partial
    })
    adjacencies <- lapply(nets, function(n) n@adjacency_matrix)
    n_donors <- vapply(nets, function(n) {
        md <- n@metadata
        if (!is.null(md$n_donors)) md$n_donors else NA_integer_
    }, integer(1))

    ## --- Identify all edges across networks ---
    union_adj <- Reduce(`+`, adjacencies)
    union_adj <- (union_adj > 0) * 1L

    ## --- For each edge, test differential ---
    edge_list <- list()
    idx <- 0L

    for (i in seq_len(p - 1L)) {
        for (j in (i + 1L):p) {
            if (union_adj[i, j] == 0) next

            idx <- idx + 1L
            pcor_vals <- vapply(precisions, function(pm) {
                pm[i, j]
            }, numeric(1))
            names(pcor_vals) <- ct_names

            present_in <- vapply(adjacencies, function(am) {
                am[i, j] != 0
            }, logical(1))

            ## Determine edge category
            if (all(present_in)) {
                category <- "shared"
            } else if (sum(present_in) == 1L) {
                category <- paste0("unique_",
                    ct_names[which(present_in)])
            } else {
                category <- "differential"
            }

            ## Statistical test
            pval <- NA_real_
            if (method == "fisher_z" && length(nets) == 2L &&
                all(is.finite(pcor_vals))) {
                r1 <- pcor_vals[1]
                r2 <- pcor_vals[2]
                n1 <- n_donors[1]
                n2 <- n_donors[2]
                if (!is.na(n1) && !is.na(n2) &&
                    n1 > 3 && n2 > 3) {
                    z1 <- atanh(min(max(r1, -0.999), 0.999))
                    z2 <- atanh(min(max(r2, -0.999), 0.999))
                    se <- sqrt(1 / (n1 - 3) + 1 / (n2 - 3))
                    z_stat <- abs(z1 - z2) / se
                    pval <- 2 * pnorm(-z_stat)
                }
            } else if (method == "fisher_z" &&
                       length(nets) > 2L && all(is.finite(pcor_vals))) {
                ## Multi-group: pairwise max test
                z_max <- 0
                for (a in seq_along(nets)[-length(nets)]) {
                    for (b in (a + 1L):length(nets)) {
                        r_a <- pcor_vals[a]
                        r_b <- pcor_vals[b]
                        n_a <- n_donors[a]
                        n_b <- n_donors[b]
                        if (is.na(n_a) || is.na(n_b) ||
                            n_a <= 3 || n_b <= 3) next
                        z_a <- atanh(min(max(r_a, -0.999),
                                         0.999))
                        z_b <- atanh(min(max(r_b, -0.999),
                                         0.999))
                        se <- sqrt(1/(n_a - 3) + 1/(n_b - 3))
                        z_max <- max(z_max,
                            abs(z_a - z_b) / se)
                    }
                }
                if (z_max > 0) pval <- 2 * pnorm(-z_max)
            }

            edge_list[[idx]] <- data.frame(
                node_i     = ref_features[i],
                node_j     = ref_features[j],
                layer_i    = nets[[1]]@node_info$omic_layer[i],
                layer_j    = nets[[1]]@node_info$omic_layer[j],
                cross_omic = nets[[1]]@node_info$omic_layer[i] !=
                    nets[[1]]@node_info$omic_layer[j],
                category   = category,
                max_diff   = max(pcor_vals) - min(pcor_vals),
                p_value    = pval,
                stringsAsFactors = FALSE
            )
        }
    }

    if (length(edge_list) == 0) {
        all_edges <- S4Vectors::DataFrame(
            node_i = character(), node_j = character(),
            layer_i = character(), layer_j = character(),
            cross_omic = logical(), category = character(),
            max_diff = numeric(), p_value = numeric(),
            p_adjusted = numeric())

        return(new("NetworkComparison",
            diff_edges        = all_edges,
            shared_edges      = all_edges,
            celltype_specific = list(),
            summary = list(
                celltypes = ct_names,
                n_edges = integer(0),
                jaccard = 0,
                method = method)))
    }

    all_df <- do.call(rbind, edge_list)
    all_df$p_adjusted <- p.adjust(all_df$p_value,
                                   method = adjust)
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
            method        = method,
            n_total_edges = nrow(all_edges),
            n_shared      = nrow(shared),
            n_differential = nrow(diff)
        )
    )
}
