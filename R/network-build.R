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
            transform      = assembled$transform,
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
#' Pre-selects features associated with exposure by stability
#' selection or univariate testing, then builds a
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
#' @param selection_method Character; \code{"stability"}
#'   (default) or \code{"univariate"}.
#' @param fdr_threshold Numeric; FDR threshold for univariate
#'   selection. Default 0.05.
#' @param pi_threshold Numeric in (0.5, 1]; minimum proportion of
#'   subsamples in which a feature must be selected under
#'   \code{selection_method = "stability"}. Default 0.6.
#' @param stability_q Integer or \code{NULL}; number of features
#'   selected in each subsample under
#'   \code{selection_method = "stability"}. The default \code{NULL}
#'   uses \eqn{\max(1, \lfloor\sqrt{(2\pi_{thr} - 1)p}\rfloor)}, the
#'   largest value for which the Meinshausen-Buhlmann bound on the
#'   expected number of falsely selected features is at most 1 (when
#'   \eqn{(2\pi_{thr} - 1)p \ge 1}).
#' @param min_cells Integer; minimum cells per pseudobulk.
#' @param network_method Character; \code{"coglasso"} (default) or
#'   \code{"block_glasso"}, the estimator for the network on the
#'   selected features. \code{"coglasso"} requires the \pkg{coglasso}
#'   package and at least one selected transcript and metabolite;
#'   otherwise the block graphical lasso is used.
#'
#' @return A \code{\linkS4class{CelltypeNetworkResult}}
#'   restricted to exposure-associated features. Its metadata record
#'   the expression transform, the selection settings and, for
#'   stability selection, \code{pfer_bound}.
#'
#' @details
#' Stage 1: Feature selection. For each transcript and
#' metabolite, tests association with exposure.
#' \code{"univariate"} uses Spearman correlation tests (asymptotic
#' \emph{t} approximation) with BH correction; features without
#' variation across donors are not selected.
#' \code{"stability"} applies stability selection (Meinshausen and
#' Buhlmann 2010) to a marginal screen: in each of 100 random
#' subsamples of half the donors with an observed exposure, the
#' \code{stability_q} features with the largest absolute Spearman
#' correlation with the exposure are selected, and features selected
#' in a proportion of at least \code{pi_threshold} of subsamples are
#' retained. The screen remains defined when features outnumber
#' donors. Under the exchangeability and better-than-random-guessing
#' assumptions of Meinshausen and Buhlmann, the expected number of
#' falsely selected features is at most
#' \eqn{q^2 / ((2\pi_{thr} - 1)p)}, reported as \code{pfer_bound}.
#' At least three features must be selected.
#'
#' Transcript values are log-CPM or VST values and therefore relative
#' to each donor's library size. When the exposure strongly changes
#' genes that carry a large share of the counts, the remaining genes
#' shift in the opposite direction and can pass the screen; with few
#' measured genes this can dominate the univariate screen.
#'
#' Stage 2: Network estimation on selected features only,
#' without residualising the exposure (to capture
#' exposure-driven rewiring), with collaborative graphical lasso
#' model selection by XStARS (\code{"coglasso"}) or the block
#' graphical lasso (\code{"block_glasso"}).
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
#' set.seed(3)
#' donors <- sprintf("D%02d", 1:20)
#' donor <- rep(donors, each = 40)
#' exposure <- stats::setNames(stats::rnorm(20), donors)
#' counts <- matrix(stats::rpois(20 * length(donor), 20), nrow = 20,
#'     dimnames = list(paste0("G", 1:20), paste0("c", seq_along(donor))))
#' ## two genes rise and two fall with the exposure
#' effect <- c(0.8, 0.8, -0.8, -0.8)
#' for (g in 1:4)
#'     counts[g, ] <- stats::rpois(length(donor),
#'         20 * exp(effect[g] * exposure[donor]))
#' sce <- SingleCellExperiment::SingleCellExperiment(
#'     assays = list(counts = counts),
#'     colData = S4Vectors::DataFrame(donor_id = donor, cell_type = "Mono"))
#' scee <- build_scee(sce, matrix(exposure, ncol = 1,
#'     dimnames = list(donors, "PM2.5")), sample_col = "donor_id")
#' metabolites <- matrix(stats::rnorm(20 * 3), nrow = 20,
#'     dimnames = list(donors, paste0("M", 1:3)))
#' metabolites[, "M1"] <- metabolites[, "M1"] + 1.5 * exposure
#' net <- run_exposure_network(scee, metabolites, celltype = "Mono",
#'     exposure = "PM2.5", sample_col = "donor_id",
#'     selection_method = "univariate", network_method = "block_glasso")
#' net
run_exposure_network <- function(scee, metabolites, celltype,
                                  exposure,
                                  celltype_col = "cell_type",
                                  sample_col = NULL,
                                  covariates = NULL,
                                  selection_method = c(
                                      "stability",
                                      "univariate"),
                                  fdr_threshold = 0.05,
                                  pi_threshold = 0.6,
                                  stability_q = NULL,
                                  min_cells = 10L,
                                  network_method = c("coglasso",
                                                     "block_glasso")) {

    selection_method <- match.arg(selection_method)
    network_method <- match.arg(network_method)

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

    pfer_bound <- NA_real_
    varies <- apply(X, 2L, function(v) isTRUE(stats::sd(v) > 0))

    if (selection_method == "univariate") {
        ## Univariate Spearman test per feature; constant features are
        ## not tested and not selected
        pvals <- rep(NA_real_, p)
        pvals[varies] <- vapply(which(varies), function(j) {
            cor.test(X[, j], y, method = "spearman",
                     exact = FALSE)$p.value
        }, numeric(1))
        padj <- p.adjust(pvals, method = "BH")
        selected <- !is.na(padj) & padj < fdr_threshold

    } else {
        ## Stability selection (Meinshausen and Buhlmann 2010) over a
        ## marginal screen, which stays defined when p exceeds n
        if (!is.numeric(pi_threshold) || length(pi_threshold) != 1L ||
            is.na(pi_threshold) || pi_threshold <= 0.5 ||
            pi_threshold > 1)
            stop("pi_threshold must be a single number in (0.5, 1].")
        observed <- which(is.finite(y))
        n_half <- floor(length(observed) / 2)
        if (n_half < 5L)
            stop("Stability selection needs at least 10 donors with an ",
                 "observed exposure; use selection_method = ",
                 "'univariate'.")
        if (is.null(stability_q))
            stability_q <- max(1, floor(sqrt((2 * pi_threshold - 1) * p)))
        if (!is.numeric(stability_q) || length(stability_q) != 1L ||
            is.na(stability_q) || stability_q < 1 || stability_q > p ||
            stability_q != round(stability_q))
            stop("stability_q must be a single integer between 1 and ",
                 "the number of features (", p, ").")
        stability_q <- as.integer(stability_q)

        n_subsamples <- 100L
        sel_freq <- numeric(p)
        for (b in seq_len(n_subsamples)) {
            idx <- sample(observed, n_half)
            sel_freq <- sel_freq + .top_marginal_features(
                X[idx, , drop = FALSE], y[idx], stability_q)
        }

        sel_prob <- sel_freq / n_subsamples
        selected <- sel_prob >= pi_threshold
        pfer_bound <- stability_q^2 / ((2 * pi_threshold - 1) * p)
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

    use_coglasso <- network_method == "coglasso" &&
        n_tx_sel > 0 && n_met_sel > 0
    if (use_coglasso && !requireNamespace("coglasso", quietly = TRUE)) {
        message("[exposomeSC] Package 'coglasso' is not installed; ",
                "using the block graphical lasso.")
        use_coglasso <- FALSE
    }
    if (use_coglasso) {
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
            fdr_threshold     = if (selection_method == "univariate")
                fdr_threshold else NA_real_,
            pi_threshold      = if (selection_method == "stability")
                pi_threshold else NA_real_,
            stability_q       = if (is.null(stability_q))
                NA_integer_ else stability_q,
            pfer_bound        = pfer_bound,
            network_method    = if (use_coglasso) "coglasso" else
                "block_glasso",
            transform         = assembled$transform,
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
