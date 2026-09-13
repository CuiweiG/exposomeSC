# R/network-temporal.R
# Temporal cross-omic network evolution

#' @include AllClasses.R
#' @include network-build.R
#' @include network-utils.R
#' @importFrom S4Vectors DataFrame
NULL

#' Temporal cross-omic network evolution
#'
#' For longitudinal designs (e.g., pre/post exposure), builds
#' networks at each time point and tracks edge dynamics:
#' which edges emerge, disappear, or persist.
#'
#' @param scee_list Named list of
#'   \code{\linkS4class{SingleCellExposomeExperiment}} objects,
#'   one per time point.
#' @param metabolites_list Named list of metabolite matrices,
#'   one per time point. Names must match \code{scee_list}.
#' @param celltype Character; cell type to analyse.
#' @param celltype_col Character; column in \code{colData}
#'   for cell types. Default \code{"cell_type"}.
#' @param sample_col Character or NULL; donor ID column.
#' @param timepoints Character vector; ordered time point
#'   labels. Default: \code{names(scee_list)}.
#' @param ... Additional arguments passed to
#'   \code{\link{run_celltype_network}}.
#'
#' @return A \code{\linkS4class{TemporalNetwork}} object
#'   containing:
#'   \describe{
#'     \item{networks}{List of CelltypeNetworkResult per
#'       time point.}
#'     \item{edge_dynamics}{DataFrame tracking each edge's
#'       presence/absence at each time point, with
#'       categorization as emerging, disappearing, persistent,
#'       transient or intermittent.}
#'     \item{timepoints}{Ordered time point labels.}
#'   }
#'
#' @details
#' Builds a separate network at each time point using
#' \code{run_celltype_network()}, then compares edge sets
#' across time. Edge dynamics are categorized as:
#' \describe{
#'   \item{persistent}{Present at all time points.}
#'   \item{emerging}{Absent at baseline, present at later
#'     time points.}
#'   \item{disappearing}{Present at baseline, absent later.}
#'   \item{transient}{Present at intermediate time points
#'     only.}
#'   \item{intermittent}{Present at the first and last time
#'     points but absent in between.}
#' }
#'
#' This mirrors the temporal network approach in Cheng et al.
#' (ES&T 2024), adapted for cell-type-resolved analysis.
#'
#' @references
#' Cheng SL et al. (2024). Multiomic signatures of traffic-related
#'   air pollution in London reveal potential short-term perturbations
#'   in gut microbiome-related pathways. \emph{Environ Sci Technol}
#'   58:8771-8782. \doi{10.1021/acs.est.3c09148}
#'
#' @export
#' @examples
#' donor_ids <- paste0("D", seq_len(8L))
#' cell_donor <- rep(donor_ids, each = 10L)
#' cell_ids <- paste0("cell", seq_along(cell_donor))
#' counts <- outer(
#'     seq_len(4L),
#'     seq_along(cell_ids),
#'     function(gene, cell) 2L + ((5L * gene + 3L * cell) %% 9L)
#' )
#' storage.mode(counts) <- "integer"
#' dimnames(counts) <- list(paste0("G", seq_len(4L)), cell_ids)
#' sce <- SingleCellExperiment::SingleCellExperiment(
#'     assays = list(counts = counts),
#'     colData = S4Vectors::DataFrame(
#'         cell_id = cell_ids,
#'         donor_id = cell_donor,
#'         cell_type = "Monocyte"
#'     )
#' )
#' exposure <- matrix(seq_len(8L), ncol = 1L,
#'     dimnames = list(donor_ids, "exposure"))
#' scee <- build_scee(sce, exposure, sample_col = "donor_id")
#' baseline_metabolites <- outer(
#'     seq_len(8L),
#'     seq_len(2L),
#'     function(donor, feature) (donor + feature)^2 / 10
#' )
#' dimnames(baseline_metabolites) <- list(donor_ids, c("M1", "M2"))
#' follow_up_metabolites <- baseline_metabolites
#' follow_up_metabolites[, "M1"] <- follow_up_metabolites[, "M1"] +
#'     seq_len(8L) / 5
#' temporal <- run_temporal_network(
#'     scee_list = list(baseline = scee, follow_up = scee),
#'     metabolites_list = list(
#'         baseline = baseline_metabolites,
#'         follow_up = follow_up_metabolites
#'     ),
#'     celltype = "Monocyte",
#'     sample_col = "donor_id",
#'     method = "block_glasso",
#'     stability = FALSE,
#'     min_cells = 5L,
#'     top_var_genes = 3L
#' )
#' temporal
run_temporal_network <- function(scee_list, metabolites_list,
                                  celltype,
                                  celltype_col = "cell_type",
                                  sample_col = NULL,
                                  timepoints = names(scee_list),
                                  ...) {

    stopifnot(is.list(scee_list))
    stopifnot(is.list(metabolites_list))
    stopifnot(length(scee_list) >= 2L)

    if (is.null(timepoints))
        timepoints <- paste0("T", seq_along(scee_list))

    if (length(timepoints) != length(scee_list))
        stop("Length of timepoints must match scee_list")

    if (length(metabolites_list) != length(scee_list))
        stop("metabolites_list must have same length ",
             "as scee_list")

    names(scee_list) <- timepoints
    names(metabolites_list) <- timepoints

    ## --- Build network at each time point ---
    networks <- list()
    for (tp in timepoints) {
        message(sprintf(
            "[exposomeSC] Building temporal network: %s", tp))
        networks[[tp]] <- run_celltype_network(
            scee_list[[tp]],
            metabolites_list[[tp]],
            celltype      = celltype,
            celltype_col  = celltype_col,
            sample_col    = sample_col,
            ...)
    }

    ## --- Track edge dynamics ---
    ## Use union of all features
    all_features <- Reduce(union, lapply(networks, function(n) {
        n@node_info$feature
    }))

    ## Build presence matrix: edges x timepoints
    edge_key_list <- list()

    for (tp in timepoints) {
        net <- networks[[tp]]
        am <- net@adjacency_matrix
        ni <- net@node_info
        p <- nrow(am)

        for (i in seq_len(p - 1L)) {
            for (j in (i + 1L):p) {
                if (am[i, j] != 0) {
                    key <- paste(sort(c(ni$feature[i],
                                         ni$feature[j])),
                                  collapse = "||")
                    if (is.null(edge_key_list[[key]])) {
                        edge_key_list[[key]] <- list(
                            node_i = ni$feature[min(i, j)],
                            node_j = ni$feature[max(i, j)],
                            layer_i = ni$omic_layer[min(i, j)],
                            layer_j = ni$omic_layer[max(i, j)],
                            cross_omic = ni$omic_layer[i] !=
                                ni$omic_layer[j],
                            present = setNames(
                                rep(FALSE, length(timepoints)),
                                timepoints)
                        )
                    }
                    edge_key_list[[key]]$present[tp] <- TRUE
                }
            }
        }
    }

    ## Build edge_dynamics DataFrame
    if (length(edge_key_list) > 0) {
        rows <- lapply(edge_key_list, function(e) {
            baseline <- e$present[1]
            final <- e$present[length(timepoints)]
            all_present <- all(e$present)
            any_present <- any(e$present)

            status <- if (all_present) {
                "persistent"
            } else if (!baseline && final) {
                "emerging"
            } else if (baseline && !final) {
                "disappearing"
            } else if (baseline && final) {
                "intermittent"
            } else {
                "transient"
            }

            row <- data.frame(
                node_i     = e$node_i,
                node_j     = e$node_j,
                layer_i    = e$layer_i,
                layer_j    = e$layer_j,
                cross_omic = e$cross_omic,
                status     = status,
                n_timepoints = sum(e$present),
                stringsAsFactors = FALSE
            )
            ## Add presence columns for each timepoint
            for (tp in timepoints) {
                row[[paste0("present_", tp)]] <- e$present[tp]
            }
            row
        })

        edge_dynamics <- S4Vectors::DataFrame(
            do.call(rbind, rows))
    } else {
        edge_dynamics <- S4Vectors::DataFrame(
            node_i       = character(),
            node_j       = character(),
            layer_i      = character(),
            layer_j      = character(),
            cross_omic   = logical(),
            status       = character(),
            n_timepoints = integer()
        )
    }

    ## Summary
    if (nrow(edge_dynamics) > 0) {
        status_tab <- table(edge_dynamics$status)
        message(sprintf(
            "[exposomeSC] Temporal dynamics: %d edges total",
            nrow(edge_dynamics)))
        for (s in names(status_tab)) {
            message(sprintf("  %s: %d", s, status_tab[s]))
        }
    }

    new("TemporalNetwork",
        networks      = networks,
        edge_dynamics = edge_dynamics,
        timepoints    = timepoints
    )
}
