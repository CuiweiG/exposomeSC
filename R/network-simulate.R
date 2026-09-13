# R/network-simulate.R
# Simulation engine for cross-omic network benchmarking

#' @include AllClasses.R
#' @include network-utils.R
#' @importFrom S4Vectors DataFrame
#' @importFrom stats rnorm rpois runif rbinom cov var
NULL

#' Simulate cross-omic network data with known structure
#'
#' Generates donor-level metabolite data and cell-type-specific transcript
#' pseudobulk counts whose latent log-expression follows a known Gaussian
#' graphical model, with cell-type heterogeneity and optional
#' exposure-dependent composition. Intended for benchmarking edge recovery
#' and false-positive control.
#'
#' @param n_donors Integer; number of donors. Default 100.
#' @param n_celltypes Integer; number of cell types. Default 4.
#' @param n_transcripts Integer; number of transcript
#'   features. Default 30.
#' @param n_metabolites Integer; number of metabolite
#'   features. Default 10.
#' @param n_cells_per_donor Integer; cells per donor, allocated to cell types
#'   by the composition. Pseudobulk depth scales with the number of cells.
#'   Default 200.
#' @param edge_density Numeric; proportion of possible edges
#'   present in the base network. Default 0.1.
#' @param shared_edge_frac Numeric; fraction of base-network edges kept by
#'   cell types 1 and 2, each of which also gains the same number of new
#'   edges. Default 0.6.
#' @param cross_omic_frac Numeric; fraction of edges that
#'   are cross-omic. Default 0.2.
#' @param exposure_effect Numeric; effect of the exposure on the mean of
#'   affected features and on composition. Default 0.5.
#' @param composition_confounding Logical; make the proportions of cell
#'   types 1 and 2 depend on the exposure. Default TRUE.
#' @param seed Integer or NULL; random seed. The caller's random number
#'   generator state is restored on return. Default NULL.
#'
#' @return List with:
#'   \describe{
#'     \item{pseudobulk}{List of count matrices (donors x transcripts), one
#'       per cell type.}
#'     \item{metabolites}{Metabolite matrix (donors x metabolites), shared by
#'       all cell types.}
#'     \item{exposure}{Named numeric vector of exposure values.}
#'     \item{composition}{Donor-by-cell-type matrix of proportions.}
#'     \item{n_cells}{Donor-by-cell-type matrix of simulated cell numbers.}
#'     \item{true_networks}{List of precision matrices, one per cell type,
#'       of the latent transcript log-expression and the metabolites
#'       conditional on the exposure.}
#'     \item{true_adjacency}{List of the corresponding adjacency matrices.}
#'     \item{celltype_names}{Character vector of cell type names.}
#'     \item{node_info}{DataFrame of feature metadata.}
#'     \item{params}{List of simulation parameters.}
#'   }
#'
#' @details
#' Simulation design:
#' \enumerate{
#'   \item A sparse positive-definite base precision matrix
#'     \eqn{\Theta} is generated over transcripts and metabolites with
#'     within-omic and cross-omic edges.
#'   \item Cell types 1 and 2 keep a fraction \code{shared_edge_frac} of
#'     the base edges and gain as many new edges; cell type 3 has an
#'     independently generated precision matrix; cell types 4 and above
#'     have no transcript-transcript or cross-omic edges.
#'   \item Metabolites are drawn once per donor from a multivariate normal
#'     distribution whose precision is the metabolite block of the base
#'     matrix, \eqn{P_m}.
#'   \item For each cell type, latent transcript log-expression is drawn
#'     from its conditional distribution given the metabolites,
#'     \eqn{N(-\Theta_{tt}^{-1}\Theta_{tm} m, \Theta_{tt}^{-1})}. The joint
#'     precision of transcripts and metabolites is then
#'     \eqn{\Theta_{tt}}, \eqn{\Theta_{tm}} and
#'     \eqn{P_m + \Theta_{mt}\Theta_{tt}^{-1}\Theta_{tm}}; it is returned in
#'     \code{true_networks}, so its metabolite block can contain edges
#'     induced by the cross-omic edges.
#'   \item The exposure shifts the mean of the first 30\% of metabolites
#'     and, in cell types 1 and 2, of the first 20\% of transcripts; it
#'     does not change any network.
#'   \item Pseudobulk counts are Poisson with mean equal to the number of
#'     cells times \eqn{\exp(\text{latent} + 1)}, which makes them
#'     overdispersed relative to a Poisson with a fixed mean.
#'   \item If \code{composition_confounding = TRUE}, the proportions of
#'     cell types 1 and 2 increase with the exposure, so pseudobulk depth
#'     varies with the exposure.
#' }
#'
#' @export
#' @examples
#' sim <- simulate_crossomic_network(
#'     n_donors = 50, n_transcripts = 20,
#'     n_metabolites = 5, seed = 42)
#' names(sim)
#' sim$true_adjacency[[1]][1:5, 1:5]
simulate_crossomic_network <- function(
        n_donors = 100L,
        n_celltypes = 4L,
        n_transcripts = 30L,
        n_metabolites = 10L,
        n_cells_per_donor = 200L,
        edge_density = 0.1,
        shared_edge_frac = 0.6,
        cross_omic_frac = 0.2,
        exposure_effect = 0.5,
        composition_confounding = TRUE,
        seed = NULL) {

    if (!is.null(seed)) .local_rng_scope(seed)
    stopifnot(n_transcripts >= 1L, n_metabolites >= 1L, n_celltypes >= 1L,
              n_donors >= 2L, n_cells_per_donor >= n_celltypes)

    p <- n_transcripts + n_metabolites
    tx <- seq_len(n_transcripts)
    met <- n_transcripts + seq_len(n_metabolites)
    ct_names <- paste0("CellType_", LETTERS[seq_len(n_celltypes)])
    donor_ids <- paste0("D", seq_len(n_donors))
    feature_names <- c(paste0("gene_", tx), paste0("metab_", seq_len(n_metabolites)))

    ## --- Generate exposure ---
    exposure <- setNames(rnorm(n_donors), donor_ids)

    ## --- Generate base precision matrix ---
    base_prec <- .generate_sparse_precision(
        p, edge_density, n_transcripts,
        n_metabolites, cross_omic_frac)

    ## --- Cell-type-specific precision matrices ---
    ct_prec <- list()
    for (ct_idx in seq_len(n_celltypes)) {
        ct_prec[[ct_names[ct_idx]]] <- if (ct_idx <= 2) {
            .perturb_precision(
                base_prec, frac_change = 1 - shared_edge_frac,
                seed = if (is.null(seed)) NULL else seed + ct_idx)
        } else if (ct_idx == 3) {
            .generate_sparse_precision(
                p, edge_density, n_transcripts,
                n_metabolites, cross_omic_frac)
        } else {
            diag(1, p)
        }
    }

    ## --- Donor-level metabolites with precision P_m ---
    met_precision <- base_prec[met, met, drop = FALSE]
    metabolites <- matrix(rnorm(n_donors * n_metabolites),
                          n_donors, n_metabolites) %*%
        chol(solve(met_precision))
    n_affected_met <- max(1L, floor(n_metabolites * 0.3))
    metabolites[, seq_len(n_affected_met)] <-
        metabolites[, seq_len(n_affected_met)] + exposure_effect * exposure
    dimnames(metabolites) <- list(donor_ids, feature_names[met])

    ## --- Composition and cell numbers ---
    comp_mat <- matrix(
        1 / n_celltypes, n_donors, n_celltypes,
        dimnames = list(donor_ids, ct_names))
    if (composition_confounding) {
        for (ct_idx in seq_len(min(2, n_celltypes))) {
            comp_mat[, ct_idx] <- comp_mat[, ct_idx] +
                exposure_effect * 0.1 * exposure
        }
        comp_mat[comp_mat < 0.01] <- 0.01
        comp_mat <- comp_mat / rowSums(comp_mat)
    }
    n_cells <- round(n_cells_per_donor * comp_mat)
    n_cells[n_cells < 1] <- 1
    storage.mode(n_cells) <- "integer"

    ## --- Transcripts conditional on metabolites, per cell type ---
    pseudobulk_list <- list()
    true_prec <- list()
    true_adj <- list()
    for (ct_idx in seq_len(n_celltypes)) {
        ct <- ct_names[ct_idx]
        theta <- ct_prec[[ct]]
        theta_tt_inv <- solve(theta[tx, tx, drop = FALSE])
        theta_tm <- theta[tx, met, drop = FALSE]
        conditional_mean <- -metabolites %*% t(theta_tm) %*% theta_tt_inv
        log_expr <- conditional_mean +
            matrix(rnorm(n_donors * n_transcripts),
                   n_donors, n_transcripts) %*% chol(theta_tt_inv)
        if (ct_idx <= 2) {
            n_affected <- max(1L, floor(n_transcripts * 0.2))
            log_expr[, seq_len(n_affected)] <-
                log_expr[, seq_len(n_affected)] + exposure_effect * exposure
        }
        mu <- n_cells[, ct] * exp(log_expr + 1)
        counts <- matrix(rpois(length(mu), mu), n_donors, n_transcripts,
                         dimnames = list(donor_ids, feature_names[tx]))
        pseudobulk_list[[ct]] <- counts

        ## Joint precision of latent transcripts and metabolites
        joint <- matrix(0, p, p, dimnames = list(feature_names, feature_names))
        joint[tx, tx] <- theta[tx, tx]
        joint[tx, met] <- theta_tm
        joint[met, tx] <- t(theta_tm)
        joint[met, met] <- met_precision + t(theta_tm) %*% theta_tt_inv %*% theta_tm
        joint[abs(joint) < 1e-12] <- 0
        true_prec[[ct]] <- joint
        true_adj[[ct]] <- .precision_to_adjacency(joint)
    }

    ## --- Node info ---
    node_info <- S4Vectors::DataFrame(
        feature = feature_names,
        omic_layer = c(rep("transcript", n_transcripts),
                       rep("metabolite", n_metabolites)),
        block = c(rep(1L, n_transcripts),
                  rep(2L, n_metabolites))
    )

    list(
        pseudobulk     = pseudobulk_list,
        metabolites    = metabolites,
        exposure       = exposure,
        composition    = comp_mat,
        n_cells        = n_cells,
        true_networks  = true_prec,
        true_adjacency = true_adj,
        celltype_names = ct_names,
        node_info      = node_info,
        params         = list(
            n_donors             = n_donors,
            n_celltypes          = n_celltypes,
            n_transcripts        = n_transcripts,
            n_metabolites        = n_metabolites,
            n_cells_per_donor    = n_cells_per_donor,
            edge_density         = edge_density,
            shared_edge_frac     = shared_edge_frac,
            cross_omic_frac      = cross_omic_frac,
            exposure_effect      = exposure_effect,
            composition_confounding = composition_confounding
        )
    )
}

# -------------------------------------------------------
# Internal helpers for simulation
# -------------------------------------------------------

#' Generate a sparse positive-definite precision matrix
#'   with block structure
#'
#' @keywords internal
#' @noRd
.generate_sparse_precision <- function(p, density,
                                        n_tx, n_met,
                                        cross_frac) {
    ## Start with identity
    prec <- diag(1, p)

    ## Total possible off-diagonal edges
    n_possible <- p * (p - 1) / 2
    n_edges <- max(1L, floor(n_possible * density))

    ## Allocate edges: cross-omic vs within-omic
    n_cross <- max(0L, floor(n_edges * cross_frac))
    n_within <- n_edges - n_cross

    ## Within-omic edges (transcript-transcript + metabolite-metabolite)
    if (n_within > 0) {
        within_pairs <- list()
        idx <- 0L
        ## Transcript-transcript
        for (i in seq_len(max(n_tx - 1L, 0L))) {
            for (j in (i + 1L):n_tx) {
                idx <- idx + 1L
                within_pairs[[idx]] <- c(i, j)
            }
        }
        ## Metabolite-metabolite
        if (n_met >= 2L) {
            for (i in (n_tx + 1L):(p - 1L)) {
                for (j in (i + 1L):p) {
                    idx <- idx + 1L
                    within_pairs[[idx]] <- c(i, j)
                }
            }
        }
        if (length(within_pairs) > 0) {
            sel <- sample(length(within_pairs),
                           min(n_within, length(within_pairs)))
            for (s in sel) {
                pair <- within_pairs[[s]]
                val <- runif(1, 0.2, 0.5) *
                    sample(c(-1, 1), 1)
                prec[pair[1], pair[2]] <- val
                prec[pair[2], pair[1]] <- val
            }
        }
    }

    ## Cross-omic edges
    if (n_cross > 0 && n_tx > 0 && n_met > 0) {
        cross_pairs <- list()
        idx <- 0L
        for (i in seq_len(n_tx)) {
            for (j in (n_tx + 1L):p) {
                idx <- idx + 1L
                cross_pairs[[idx]] <- c(i, j)
            }
        }
        sel <- sample(length(cross_pairs),
                       min(n_cross, length(cross_pairs)))
        for (s in sel) {
            pair <- cross_pairs[[s]]
            val <- runif(1, 0.15, 0.4) *
                sample(c(-1, 1), 1)
            prec[pair[1], pair[2]] <- val
            prec[pair[2], pair[1]] <- val
        }
    }

    ## Ensure positive definiteness by adding to diagonal
    eig <- eigen(prec, symmetric = TRUE, only.values = TRUE)
    if (min(eig$values) < 0.1) {
        prec <- prec + diag(abs(min(eig$values)) + 0.1, p)
    }

    prec
}

#' Perturb a precision matrix by adding/removing edges
#'
#' @keywords internal
#' @noRd
.perturb_precision <- function(prec, frac_change = 0.4,
                                seed = NULL) {
    ## The caller keeps drawing from this stream after we return, so the
    ## scope ends when the caller exits rather than here.
    if (!is.null(seed)) .local_rng_scope(seed, envir = parent.frame())
    p <- nrow(prec)

    ## Find existing edges
    ut <- upper.tri(prec)
    edges <- which(ut & prec != 0, arr.ind = TRUE)
    non_edges <- which(ut & prec == 0, arr.ind = TRUE)

    n_change <- max(1L, floor(nrow(edges) * frac_change))

    ## Remove some edges
    if (nrow(edges) > 0 && n_change > 0) {
        to_remove <- sample(nrow(edges),
                             min(n_change, nrow(edges)))
        for (r in to_remove) {
            i <- edges[r, 1]
            j <- edges[r, 2]
            prec[i, j] <- 0
            prec[j, i] <- 0
        }
    }

    ## Add some new edges
    if (nrow(non_edges) > 0 && n_change > 0) {
        to_add <- sample(nrow(non_edges),
                          min(n_change, nrow(non_edges)))
        for (a in to_add) {
            i <- non_edges[a, 1]
            j <- non_edges[a, 2]
            val <- runif(1, 0.15, 0.4) *
                sample(c(-1, 1), 1)
            prec[i, j] <- val
            prec[j, i] <- val
        }
    }

    ## Re-ensure positive definiteness
    eig <- eigen(prec, symmetric = TRUE, only.values = TRUE)
    if (min(eig$values) < 0.1) {
        prec <- prec + diag(abs(min(eig$values)) + 0.1,
                             nrow(prec))
    }

    prec
}
