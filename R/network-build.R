# R/network-build.R
# Core: build cross-omic networks per cell type

#' @include AllClasses.R
#' @include network-utils.R
#' @importFrom S4Vectors DataFrame
#' @importFrom stats cov cor var
NULL

#' Build cell-type-resolved cross-omic network
#'
#' Given a SCEE with pseudobulk transcriptomics and donor-level
#' metabolomics, constructs a block-calibrated conditional
#' independence network for a specified cell type.
#'
#' @param scee \code{\linkS4class{SingleCellExposomeExperiment}}.
#' @param metabolites Numeric matrix; rows = donors,
#'   cols = metabolites. Row names must be donor IDs.
#' @param celltype Character; which cell type to build
#'   the network for.
#' @param celltype_col Character; column in \code{colData}
#'   identifying cell types. Default \code{"cell_type"}.
#' @param sample_col Character or NULL. Column for donor IDs.
#'   If NULL, inferred from \code{sampleMap}.
#' @param exposure Character or NULL; exposure variable to
#'   condition on (regressed out). Default NULL.
#' @param covariates Character vector or NULL; adjustment
#'   covariates. Default NULL.
#' @param method Character; \code{"coglasso"} (default) or
#'   \code{"block_glasso"}. The coglasso method uses
#'   collaborative graphical lasso with block-calibrated
#'   penalties (Albanese et al. 2024).
#' @param stability Logical; use stability selection via
#'   XStARS for robust edge selection. Default TRUE.
#' @param nlambda_w Integer; within-omic lambda grid size.
#'   Default 15.
#' @param nlambda_b Integer; between-omic lambda grid size.
#'   Default 15.
#' @param subsample_ratio Numeric; proportion of samples for
#'   stability selection subsamples. Default 0.8.
#' @param min_cells Integer; minimum cells per pseudobulk
#'   sample. Default 10.
#' @param top_var_genes Integer or NULL; restrict to top N
#'   most variable genes for dimensionality reduction.
#'   Default NULL (all genes). Recommended: 50-200.
#' @param lambda_w Numeric; within-block penalty for
#'   block_glasso method. Default 0.3.
#' @param lambda_b Numeric; between-block penalty for
#'   block_glasso method. Default 0.5.
#' @param precomputed_network Optional; a pre-estimated network
#'   to wrap into a CelltypeNetworkResult without re-estimation.
#'   Accepts: (1) a tempoNet \code{StableNetwork} object,
#'   (2) a named list with \code{adjacency} (binary matrix),
#'   \code{stability} (numeric matrix, optional),
#'   \code{precision} (numeric matrix, optional), and
#'   \code{feature_names} (character vector).
#'   When provided, \code{metabolites}, \code{method},
#'   \code{stability}, and lambda parameters are ignored.
#'   The \code{scee} is still required for metadata extraction.
#'
#' @return A \code{\linkS4class{CelltypeNetworkResult}} object.
#'
#' @details
#' The workflow proceeds as follows:
#' \enumerate{
#'   \item Pseudobulk aggregation for the specified cell type,
#'     reusing existing \pkg{exposomeSC} infrastructure.
#'   \item Variance-stabilising transformation (DESeq2 VST)
#'     on counts, or log-CPM fallback.
#'   \item Residualization: if \code{exposure} or
#'     \code{covariates} are specified, regresses them out
#'     from both transcriptomic and metabolomic matrices.
#'     The resulting network captures the exposure-adjusted
#'     conditional dependency structure.
#'   \item Column-binding: \code{[VST residuals |
#'     metabolite residuals]}.
#'   \item Block-calibrated graphical lasso estimation.
#'     With \code{coglasso}: collaborative graphical lasso
#'     with separate within/between-omic penalties and
#'     XStARS stability selection. With \code{block_glasso}:
#'     thresholded block-penalized inverse covariance.
#'   \item Returns a \code{CelltypeNetworkResult} object.
#' }
#'
#' @references
#' Albanese A, Kohlen W, Behrouzi P (2024). Collaborative
#'   graphical lasso. \emph{arXiv} 2403.18602.
#'
#' Liu H, Roeder K, Wasserman L (2010). Stability approach
#'   to regularization selection (StARS) for high dimensional
#'   graphical models. \emph{Advances in Neural Information
#'   Processing Systems}, 23.
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
#' metabolites <- outer(
#'     seq_len(8L),
#'     seq_len(2L),
#'     function(donor, feature) (donor + feature)^2 / 10
#' )
#' dimnames(metabolites) <- list(donor_ids, c("M1", "M2"))
#' network <- run_celltype_network(
#'     scee,
#'     metabolites,
#'     celltype = "Monocyte",
#'     sample_col = "donor_id",
#'     method = "block_glasso",
#'     stability = FALSE,
#'     min_cells = 5L,
#'     top_var_genes = 3L
#' )
#' network
run_celltype_network <- function(scee, metabolites = NULL, celltype,
                                  celltype_col = "cell_type",
                                  sample_col = NULL,
                                  exposure = NULL,
                                  covariates = NULL,
                                  method = c("coglasso",
                                             "block_glasso"),
                                  stability = TRUE,
                                  nlambda_w = 15L,
                                  nlambda_b = 15L,
                                  subsample_ratio = 0.8,
                                  min_cells = 10L,
                                  top_var_genes = NULL,
                                  lambda_w = 0.3,
                                  lambda_b = 0.5,
                                  precomputed_network = NULL) {

    stopifnot(is(scee, "SingleCellExposomeExperiment"))

    ## --- Handle precomputed network (tempoNet / list) ---
    if (!is.null(precomputed_network)) {
        return(.wrap_precomputed(precomputed_network, celltype))
    }

    method <- match.arg(method)
    stopifnot(is.matrix(metabolites))

    ## --- Assemble cross-omic data ---
    assembled <- .assemble_crossomic(
        scee, metabolites, celltype,
        celltype_col = celltype_col,
        sample_col   = sample_col,
        exposure     = exposure,
        covariates   = covariates,
        min_cells    = min_cells,
        top_var_genes = top_var_genes)

    X <- assembled$data_matrix
    n <- nrow(X)
    p <- ncol(X)
    blocks <- .make_block_indicator(
        assembled$n_transcripts,
        assembled$n_metabolites)

    message(sprintf(
        "[exposomeSC] Building %s network for '%s': ",
        method, celltype),
        sprintf("n=%d donors, p=%d features ",
            n, p),
        sprintf("(%d transcripts + %d metabolites)",
            assembled$n_transcripts,
            assembled$n_metabolites))

    ## --- Estimate network ---
    if (method == "coglasso") {
        if (!requireNamespace("coglasso", quietly = TRUE))
            stop("Package 'coglasso' required for method='coglasso'. ",
                 "Install with: install.packages('coglasso')")

        ## coglasso requires data matrix and block indicator
        cg_result <- coglasso::bs(
            X,
            p = assembled$n_transcripts,
            nlambda_w = nlambda_w,
            nlambda_b = nlambda_b
        )

        ## Select best model via XStARS
        cg_stars <- coglasso::xstars(
            cg_result,
            rep_num = ceiling(1 / (1 - subsample_ratio)),
            stars_thresh = 0.1
        )

        ## Extract precision and adjacency
        ## coglasso >= 1.1.0 uses sel_icov / sel_adj
        precision <- as.matrix(cg_stars$sel_icov)
        adjacency <- as.matrix(cg_stars$sel_adj)
        diag(adjacency) <- 0L

        ## Stability scores (variability from XStARS)
        stab_mat <- if (stability && !is.null(cg_stars$sel_variability)) {
            tryCatch(
                as.matrix(cg_stars$sel_variability),
                error = function(e) matrix(nrow = 0, ncol = 0)
            )
        } else {
            matrix(nrow = 0, ncol = 0)
        }

    } else {
        ## block_glasso fallback
        S <- cov(X)
        bg <- .block_glasso(S, blocks,
                             lambda_w = lambda_w,
                             lambda_b = lambda_b)
        precision <- bg$wi
        adjacency <- bg$adj
        stab_mat <- matrix(nrow = 0, ncol = 0)
    }

    ## Name dimensions
    feat_names <- assembled$node_info$feature
    dimnames(precision) <- list(feat_names, feat_names)
    dimnames(adjacency) <- list(feat_names, feat_names)
    if (nrow(stab_mat) > 0 &&
        nrow(stab_mat) == length(feat_names) &&
        ncol(stab_mat) == length(feat_names))
        dimnames(stab_mat) <- list(feat_names, feat_names)

    ## --- Build result object ---
    edge_counts <- .count_edge_types(adjacency,
                                      assembled$node_info)

    new("CelltypeNetworkResult",
        precision_matrix = precision,
        adjacency_matrix = adjacency,
        stability_scores = stab_mat,
        node_info        = assembled$node_info,
        celltype         = celltype,
        method           = method,
        metadata         = list(
            n_donors       = n,
            n_features     = p,
            n_transcripts  = assembled$n_transcripts,
            n_metabolites  = assembled$n_metabolites,
            exposure       = exposure,
            covariates     = covariates,
            edge_counts    = edge_counts,
            donors         = assembled$donors,
            nlambda_w      = nlambda_w,
            nlambda_b      = nlambda_b
        )
    )
}


#' Exposure-driven network: exposure-associated features only
#'
#' Pre-selects features associated with exposure via stability
#' selection lasso or univariate testing, then builds a
#' cross-omic network on the selected features only. This
#' mirrors the two-stage approach in Cheng et al. (2024).
#'
#' @param scee \code{\linkS4class{SingleCellExposomeExperiment}}.
#' @param metabolites Numeric matrix; rows = donors,
#'   cols = metabolites.
#' @param celltype Character; cell type to analyse.
#' @param exposure Character; exposure variable name.
#' @param celltype_col Character; cell type column.
#'   Default \code{"cell_type"}.
#' @param sample_col Character or NULL; donor ID column.
#' @param covariates Character vector or NULL.
#' @param selection_method Character; \code{"stability_lasso"}
#'   (default) or \code{"univariate"}.
#' @param fdr_threshold Numeric; FDR threshold for univariate
#'   selection. Default 0.05.
#' @param pi_threshold Numeric; stability threshold for lasso
#'   selection. Default 0.6 (Meinshausen & Buhlmann 2010).
#' @param min_cells Integer; minimum cells per pseudobulk.
#' @param ... Additional arguments passed to
#'   \code{run_celltype_network}.
#'
#' @return A \code{\linkS4class{CelltypeNetworkResult}}
#'   restricted to exposure-associated features.
#'
#' @details
#' Stage 1: Feature selection. For each transcript and
#' metabolite, tests association with exposure.
#' \code{"univariate"} uses correlation tests with BH
#' correction. \code{"stability_lasso"} uses repeated
#' subsampling with lasso regression, selecting features
#' with selection probability above \code{pi_threshold}.
#'
#' Stage 2: Network estimation on selected features only,
#' WITHOUT residualizing exposure (to capture
#' exposure-driven rewiring).
#'
#' @references
#' Cheng SL et al. (2024). Multiomic signatures of traffic-related
#'   air pollution in London reveal potential short-term perturbations
#'   in gut microbiome-related pathways. \emph{Environ Sci Technol}
#'   58:8771-8782. \doi{10.1021/acs.est.3c09148}
#'
#' Meinshausen N, Buhlmann P (2010). Stability selection.
#'   \emph{J R Stat Soc Series B}, 72(4):417-473.
#'
#' @export
#' @examples
#' \dontrun{
#' net <- run_exposure_network(scee, metab,
#'     celltype = "Monocyte",
#'     exposure = "PM2.5",
#'     selection_method = "univariate")
#' }
run_exposure_network <- function(scee, metabolites, celltype,
                                  exposure,
                                  celltype_col = "cell_type",
                                  sample_col = NULL,
                                  covariates = NULL,
                                  selection_method = c(
                                      "stability_lasso",
                                      "univariate"),
                                  fdr_threshold = 0.05,
                                  pi_threshold = 0.6,
                                  min_cells = 10L,
                                  ...) {

    selection_method <- match.arg(selection_method)

    ## Assemble data WITHOUT residualization (we want exposure signal)
    assembled <- .assemble_crossomic(
        scee, metabolites, celltype,
        celltype_col = celltype_col,
        sample_col   = sample_col,
        exposure     = NULL,  # No residualization
        covariates   = covariates,
        min_cells    = min_cells)

    X <- assembled$data_matrix
    donors <- assembled$donors

    ## Get exposure values for shared donors
    exp_data <- slot(scee, "exposureData")
    if (!exposure %in% colnames(exp_data))
        stop("Exposure '", exposure, "' not in exposureData")
    y <- exp_data[donors, exposure]

    ## --- Stage 1: Feature selection ---
    p <- ncol(X)
    selected <- rep(FALSE, p)

    if (selection_method == "univariate") {
        ## Univariate correlation test per feature
        pvals <- vapply(seq_len(p), function(j) {
            ct <- cor.test(X[, j], y, method = "spearman")
            ct$p.value
        }, numeric(1))
        padj <- p.adjust(pvals, method = "BH")
        selected <- padj < fdr_threshold

    } else {
        ## Stability selection with lasso
        n_sub <- max(floor(nrow(X) * 0.8), 5L)
        n_boot <- 100L
        sel_freq <- numeric(p)

        for (b in seq_len(n_boot)) {
            idx <- sample(nrow(X), n_sub, replace = FALSE)
            X_sub <- X[idx, , drop = FALSE]
            y_sub <- y[idx]

            ## Simple lasso via coordinate descent
            ## Use scaled cross-validation lambda
            sel_freq <- sel_freq + tryCatch({
                fit <- lm(y_sub ~ X_sub)
                coefs <- coef(fit)[-1]  # drop intercept
                ## Select features with |t| > 2
                se <- summary(fit)$coefficients[-1, 2]
                se[se == 0] <- Inf
                t_stat <- abs(coefs / se)
                as.numeric(t_stat > 2)
            }, error = function(e) rep(0, p))
        }

        sel_prob <- sel_freq / n_boot
        selected <- sel_prob >= pi_threshold
    }

    n_selected <- sum(selected)
    if (n_selected < 3L)
        stop("Only ", n_selected, " features selected. ",
             "Try relaxing thresholds ",
             "(fdr_threshold or pi_threshold).")

    message(sprintf(
        "[exposomeSC] Exposure network: %d/%d features selected",
        n_selected, p))

    ## --- Stage 2: Build network on selected features ---
    ## Subset metabolites to only selected metabolite features
    tx_idx <- which(assembled$node_info$omic_layer == "transcript")
    met_idx <- which(assembled$node_info$omic_layer == "metabolite")

    sel_tx <- intersect(which(selected), tx_idx)
    sel_met <- intersect(which(selected), met_idx)

    ## Get original feature names
    sel_tx_names <- assembled$node_info$feature[sel_tx]
    sel_met_names <- assembled$node_info$feature[sel_met]

    ## Build subset matrices
    sel_X <- X[, selected, drop = FALSE]
    sel_blocks <- assembled$node_info$block[selected]
    sel_node_info <- assembled$node_info[selected, ]

    ## Estimate network on subset
    S <- cov(sel_X)
    n_tx_sel <- sum(sel_node_info$omic_layer == "transcript")
    n_met_sel <- sum(sel_node_info$omic_layer == "metabolite")

    if (requireNamespace("coglasso", quietly = TRUE) &&
        n_tx_sel > 0 && n_met_sel > 0) {
        cg_result <- coglasso::bs(
            sel_X,
            p = n_tx_sel,
            nlambda_w = 10L,
            nlambda_b = 10L
        )
        cg_stars <- coglasso::xstars(cg_result)
        precision <- as.matrix(cg_stars$sel_icov)
        adjacency <- as.matrix(cg_stars$sel_adj)
        diag(adjacency) <- 0L
        stab_mat <- matrix(nrow = 0, ncol = 0)
    } else {
        bg <- .block_glasso(S, sel_blocks)
        precision <- bg$wi
        adjacency <- bg$adj
        stab_mat <- matrix(nrow = 0, ncol = 0)
    }

    feat_names <- sel_node_info$feature
    dimnames(precision) <- list(feat_names, feat_names)
    dimnames(adjacency) <- list(feat_names, feat_names)

    edge_counts <- .count_edge_types(adjacency, sel_node_info)

    new("CelltypeNetworkResult",
        precision_matrix = precision,
        adjacency_matrix = adjacency,
        stability_scores = stab_mat,
        node_info        = sel_node_info,
        celltype         = celltype,
        method           = "exposure_driven",
        metadata         = list(
            n_donors          = length(donors),
            n_features        = n_selected,
            n_transcripts     = n_tx_sel,
            n_metabolites     = n_met_sel,
            exposure          = exposure,
            selection_method  = selection_method,
            n_total_features  = p,
            edge_counts       = edge_counts
        )
    )
}

# ============================================================
# Wrap a precomputed network into CelltypeNetworkResult
# ============================================================

#' @keywords internal
.wrap_precomputed <- function(net, celltype) {
    ## Accept tempoNet StableNetwork or plain list
    if (is(net, "StableNetwork")) {
        adj   <- as.matrix(net@adjacency)
        prec  <- as.matrix(net@precision)
        stab  <- as.matrix(net@stability_scores)
        fnames <- colnames(adj)
        if (is.null(fnames)) fnames <- paste0("F", seq_len(ncol(adj)))
        ## Infer block structure if available
        if (!is.null(net@block_sizes) && length(net@block_sizes) > 0) {
            blk <- rep(names(net@block_sizes), net@block_sizes)
        } else {
            blk <- rep("unknown", length(fnames))
        }
    } else if (is.list(net)) {
        if (is.null(net$adjacency))
            stop("precomputed_network list must have '$adjacency'.",
                 call. = FALSE)
        adj  <- as.matrix(net$adjacency)
        prec <- if (!is.null(net$precision))
            as.matrix(net$precision) else matrix(NA_real_, nrow(adj), ncol(adj))
        stab <- if (!is.null(net$stability))
            as.matrix(net$stability) else matrix(nrow = 0, ncol = 0)
        fnames <- if (!is.null(net$feature_names))
            net$feature_names else colnames(adj)
        if (is.null(fnames)) fnames <- paste0("F", seq_len(ncol(adj)))
        ## Block info
        if (!is.null(net$block_sizes)) {
            blk <- rep(names(net$block_sizes), net$block_sizes)
        } else if (!is.null(net$feature_blocks)) {
            blk <- net$feature_blocks
        } else {
            blk <- rep("unknown", length(fnames))
        }
    } else {
        stop("precomputed_network must be a StableNetwork or a list ",
             "with $adjacency.", call. = FALSE)
    }

    p <- length(fnames)
    dimnames(adj)  <- list(fnames, fnames)
    dimnames(prec) <- list(fnames, fnames)
    if (nrow(stab) == p && ncol(stab) == p)
        dimnames(stab) <- list(fnames, fnames)

    node_info <- S4Vectors::DataFrame(
        feature = fnames,
        omic    = blk,
        index   = seq_len(p)
    )

    n_edges <- sum(adj[upper.tri(adj)] == 1, na.rm = TRUE)
    message(sprintf(
        "[exposomeSC] Wrapping precomputed network for '%s': ",
        celltype),
        sprintf("p=%d features, %d edges", p, n_edges))

    new("CelltypeNetworkResult",
        precision_matrix = prec,
        adjacency_matrix = adj,
        stability_scores = stab,
        node_info        = node_info,
        celltype         = celltype,
        method           = "precomputed",
        metadata         = list(
            n_features = p,
            n_edges    = n_edges,
            source     = if (is(net, "StableNetwork"))
                "tempoNet" else "list"
        )
    )
}
