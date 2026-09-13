# R/network-simulate.R
# Simulation engine for cross-omic network benchmarking

#' @include AllClasses.R
#' @include network-utils.R
#' @importFrom S4Vectors DataFrame
#' @importFrom stats rnorm rpois runif rbinom cov var
NULL

#' Simulate cross-omic network data with known structure
#'
#' Generates realistic multi-omics data with known network
#' structure, cell-type heterogeneity, and composition
#' confounding. Used for benchmarking edge recovery and
#' false positive control.
#'
#' @param n_donors Integer; number of donors. Default 100.
#' @param n_celltypes Integer; number of cell types. Default 4.
#' @param n_transcripts Integer; number of transcript
#'   features. Default 30.
#' @param n_metabolites Integer; number of metabolite
#'   features. Default 10.
#' @param n_cells_per_donor Integer; average cells per donor
#'   (before cell type allocation). Default 200.
#' @param edge_density Numeric; proportion of possible edges
#'   present in the true network. Default 0.1.
#' @param shared_edge_frac Numeric; fraction of edges shared
#'   across all cell types. Default 0.6.
#' @param cross_omic_frac Numeric; fraction of edges that
#'   are cross-omic. Default 0.2.
#' @param exposure_effect Numeric; effect size of exposure
#'   on cell composition and gene expression. Default 0.5.
#' @param composition_confounding Logical; simulate
#'   exposure-driven composition changes. Default TRUE.
#' @param seed Integer or NULL; random seed. Default NULL.
#'
#' @return List with:
#'   \describe{
#'     \item{scee_list}{List of mock SCEEs per cell type
#'       (simplified for benchmarking).}
#'     \item{metabolites}{Metabolite matrix (donors x metabolites).}
#'     \item{exposure}{Named numeric vector of exposure values.}
#'     \item{true_networks}{List of true precision matrices
#'       per cell type.}
#'     \item{true_adjacency}{List of true adjacency matrices
#'       per cell type.}
#'     \item{celltype_names}{Character vector of cell type names.}
#'     \item{node_info}{DataFrame of feature metadata.}
#'     \item{params}{List of simulation parameters.}
#'   }
#'
#' @details
#' Simulation design:
#' \enumerate{
#'   \item Generate a shared "base" precision matrix with
#'     block structure (transcript-transcript,
#'     metabolite-metabolite, cross-omic).
#'   \item For each cell type, create a variant:
#'     cell types 1 and 2 share 60\% of edges; cell type 3
#'     has unique edges; cell type 4 is a null (no
#'     exposure effect on network).
#'   \item Generate pseudobulk counts from negative binomial
#'     (mimicking real scRNA-seq aggregation).
#'   \item Generate metabolite values from multivariate
#'     normal with the cell-type-specific covariance.
#'   \item If \code{composition_confounding=TRUE}, make
#'     cell proportions dependent on exposure, creating
#'     spurious bulk-level correlations.
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

    p <- n_transcripts + n_metabolites
    ct_names <- paste0("CellType_", LETTERS[seq_len(n_celltypes)])
    donor_ids <- paste0("D", seq_len(n_donors))

    ## --- Generate exposure ---
    exposure <- setNames(rnorm(n_donors), donor_ids)

    ## --- Generate base precision matrix ---
    base_prec <- .generate_sparse_precision(
        p, edge_density, n_transcripts,
        n_metabolites, cross_omic_frac)

    ## --- Cell-type-specific variants ---
    true_prec <- list()
    true_adj <- list()

    for (ct_idx in seq_len(n_celltypes)) {
        if (ct_idx <= 2) {
            ## Cell types 1-2: share base with minor perturbation
            prec_ct <- .perturb_precision(
                base_prec, frac_change = 1 - shared_edge_frac,
                seed = seed + ct_idx)
        } else if (ct_idx == 3) {
            ## Cell type 3: unique edges
            prec_ct <- .generate_sparse_precision(
                p, edge_density, n_transcripts,
                n_metabolites, cross_omic_frac)
        } else {
            ## Cell type 4+: null (diagonal)
            prec_ct <- diag(1, p)
        }

        true_prec[[ct_names[ct_idx]]] <- prec_ct
        true_adj[[ct_names[ct_idx]]] <-
            .precision_to_adjacency(prec_ct)
    }

    ## --- Generate metabolite data ---
    ## Use base precision to generate correlated metabolites
    met_cov <- tryCatch(
        solve(base_prec[
            (n_transcripts + 1):p,
            (n_transcripts + 1):p]),
        error = function(e) diag(1, n_metabolites)
    )

    ## Ensure positive definiteness
    eig <- eigen(met_cov, symmetric = TRUE)
    eig$values <- pmax(eig$values, 0.01)
    met_cov <- eig$vectors %*% diag(eig$values) %*%
        t(eig$vectors)

    L_met <- chol(met_cov)
    metabolites <- matrix(rnorm(n_donors * n_metabolites),
                           n_donors, n_metabolites) %*% L_met
    ## Add exposure effect on some metabolites
    n_affected_met <- max(1L,
        floor(n_metabolites * 0.3))
    for (j in seq_len(n_affected_met)) {
        metabolites[, j] <- metabolites[, j] +
            exposure_effect * exposure
    }
    dimnames(metabolites) <- list(donor_ids,
        paste0("metab_", seq_len(n_metabolites)))

    ## --- Generate pseudobulk counts ---
    ## For simplicity, create count matrices per cell type
    ## based on the true covariance structure

    ## Cell composition (exposure-dependent if confounding)
    comp_mat <- matrix(
        1 / n_celltypes, n_donors, n_celltypes,
        dimnames = list(donor_ids, ct_names))
    if (composition_confounding) {
        for (ct_idx in seq_len(min(2, n_celltypes))) {
            comp_mat[, ct_idx] <- comp_mat[, ct_idx] +
                exposure_effect * 0.1 * exposure
        }
        ## Normalise to sum to 1
        comp_mat <- comp_mat / rowSums(comp_mat)
        comp_mat[comp_mat < 0.01] <- 0.01
        comp_mat <- comp_mat / rowSums(comp_mat)
    }

    ## Generate count data per cell type
    pseudobulk_list <- list()
    for (ct in ct_names) {
        prec_ct <- true_prec[[ct]]
        tx_prec <- prec_ct[seq_len(n_transcripts),
                           seq_len(n_transcripts)]
        tx_cov <- tryCatch(
            solve(tx_prec),
            error = function(e) diag(1, n_transcripts)
        )
        eig2 <- eigen(tx_cov, symmetric = TRUE)
        eig2$values <- pmax(eig2$values, 0.01)
        tx_cov <- eig2$vectors %*% diag(eig2$values) %*%
            t(eig2$vectors)

        L_tx <- chol(tx_cov)
        log_expr <- matrix(rnorm(n_donors * n_transcripts),
                            n_donors, n_transcripts) %*% L_tx
        ## Add exposure effect for cell types 1-2
        ct_idx <- match(ct, ct_names)
        if (ct_idx <= 2) {
            n_affected <- max(1, floor(n_transcripts * 0.2))
            for (j in seq_len(n_affected)) {
                log_expr[, j] <- log_expr[, j] +
                    exposure_effect * exposure
            }
        }
        ## Convert to counts (negative binomial)
        mu <- exp(log_expr + 5)  # base expression ~150
        counts <- matrix(
            vapply(as.vector(mu), function(m) {
                rpois(1, m)
            }, integer(1)),
            n_donors, n_transcripts)
        dimnames(counts) <- list(donor_ids,
            paste0("gene_", seq_len(n_transcripts)))
        pseudobulk_list[[ct]] <- counts
    }

    ## --- Node info ---
    node_info <- S4Vectors::DataFrame(
        feature = c(paste0("gene_", seq_len(n_transcripts)),
                    paste0("metab_", seq_len(n_metabolites))),
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
        true_networks  = true_prec,
        true_adjacency = true_adj,
        celltype_names = ct_names,
        node_info      = node_info,
        params         = list(
            n_donors             = n_donors,
            n_celltypes          = n_celltypes,
            n_transcripts        = n_transcripts,
            n_metabolites        = n_metabolites,
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
        for (i in seq_len(n_tx - 1L)) {
            for (j in (i + 1L):n_tx) {
                idx <- idx + 1L
                within_pairs[[idx]] <- c(i, j)
            }
        }
        ## Metabolite-metabolite
        for (i in (n_tx + 1L):(p - 1L)) {
            for (j in (i + 1L):p) {
                idx <- idx + 1L
                within_pairs[[idx]] <- c(i, j)
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
