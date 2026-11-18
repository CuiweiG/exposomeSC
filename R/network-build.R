# R/network-build.R
# Core: build cross-omic networks per cell type

#' @include AllClasses.R
#' @include network-utils.R
#' @importFrom S4Vectors DataFrame
#' @importFrom stats cov cor cor.test var
NULL

#' Build cell-type-resolved cross-omic network
#'
#' Given a SCEE with pseudobulk transcriptomics and donor-level
#' metabolomics, constructs a block-calibrated conditional
#' independence network for a specified cell type.
#'
#' @param x \code{\linkS4class{SingleCellExposomeExperiment}}. Each
#'   cell's donor is read from the \code{colData} column named by
#'   \code{\link{sampleColumn}(x)}.
#' @param metabolites Numeric matrix; rows = donors,
#'   cols = metabolites. Row names must be donor IDs. A matrix without
#'   columns, or a missing or infinite value for an analysed donor, is an
#'   error, and so is a column with the name of another column or of a
#'   gene that becomes a node: two nodes cannot have the same name.
#'   The matrix is needed to estimate a network; the default,
#'   \code{NULL}, is for a call with \code{precomputed_network}, which
#'   does not use it.
#' @param celltype Character; which cell type to build
#'   the network for.
#' @param celltype_col Character; column in \code{colData}
#'   identifying cell types. Default \code{"cell_type"}.
#' @param exposure Character or NULL; exposure variable to
#'   condition on (regressed out), a column of \code{exposureData(x)};
#'   another name is an error. Default NULL.
#' @param covariates Character vector or NULL; adjustment
#'   covariates, columns of \code{exposureData(x)}; another name is an
#'   error. Donors with a missing value of \code{exposure} or of a
#'   covariate are left out, with a message. Default NULL.
#' @param method Character; \code{"coglasso"} (default) or
#'   \code{"block_glasso"}. The coglasso method uses
#'   collaborative graphical lasso with block-calibrated
#'   penalties (Albanese et al. 2026) and needs the \pkg{coglasso}
#'   package, version 1.1.0 or later.
#' @param stability Logical; for \code{method = "coglasso"}, whether
#'   to return the edge-level selection frequencies from XStARS in
#'   \code{stability_scores}. The coglasso penalties are chosen by
#'   XStARS either way. \code{block_glasso} uses the fixed penalties
#'   \code{lambda_w} and \code{lambda_b}, performs no stability
#'   selection and ignores this argument. Default TRUE.
#' @param nlambda_w Integer; within-omic lambda grid size.
#'   Default 15.
#' @param nlambda_b Integer; between-omic lambda grid size.
#'   Default 15.
#' @param subsample_ratio Numeric strictly between 0 and 1;
#'   proportion of donors drawn in each XStARS subsample, passed to
#'   \code{coglasso::xstars()} as \code{stars_subsample_ratio}.
#'   Default 0.8.
#' @param rep_num Integer, at least 2; number of XStARS subsamples.
#'   Default 20, as in \code{coglasso::xstars()}.
#' @param min_cells Integer; minimum cells per pseudobulk
#'   sample. Default 10.
#' @param top_var_genes Integer or NULL; restrict to top N
#'   most variable genes for dimensionality reduction.
#'   Default NULL (all genes). Recommended: 50-200.
#' @param lambda_w Numeric; graphical-lasso penalty for pairs of
#'   features in the same omic layer, on the correlation scale, for
#'   the block_glasso method. Default 0.3.
#' @param lambda_b Numeric; graphical-lasso penalty for pairs of
#'   features in different omic layers, on the correlation scale,
#'   for the block_glasso method. Default 0.5.
#' @param precomputed_network Optional named list holding a
#'   pre-estimated network to wrap into a CelltypeNetworkResult
#'   without re-estimation, with elements \code{adjacency} (symmetric
#'   binary matrix, required), \code{precision} (numeric matrix, optional),
#'   \code{stability} (numeric matrix, optional),
#'   \code{feature_names} (character vector of unique, non-missing and
#'   non-empty names, optional; taken from the column names of
#'   \code{adjacency} when absent), and either
#'   \code{block_sizes} (a named vector giving the size of each omic
#'   block) or \code{feature_blocks} (a character vector labelling
#'   each feature), both optional. The block labels become the
#'   \code{omic_layer} column of the node table, so use
#'   \code{"transcript"} and \code{"metabolite"} for the summaries and
#'   plots to recognise the layers. A logical adjacency is stored as 0/1
#'   and its diagonal is set to zero.
#'   When provided, nothing is estimated, and of the other arguments only
#'   \code{celltype} is used. \code{x} must still be a
#'   \code{SingleCellExposomeExperiment}.
#' @param transform Character; the transform of the pseudobulk counts.
#'   \code{"vst"} (default) is the variance-stabilising transformation of
#'   \code{DESeq2::vst()}, blind to the design, and needs \pkg{DESeq2} and
#'   at least 1000 genes with a mean normalised count above 5;
#'   \code{"log_cpm"} is the base-2 logarithm of the counts per million
#'   plus one. The call stops when \code{"vst"} cannot be computed; it does
#'   not switch to the other transform. One choice stays with DESeq2: when
#'   the parametric fit of the dispersion trend fails, \code{vst()} uses a
#'   local fit, and a message says so.
#' @param genes Character vector or \code{NULL}; row names of \code{x}, the
#'   genes that become the transcript nodes of the network, in the order
#'   given. The transform is computed from all genes of \code{x} before
#'   these are taken, so the values of a gene do not depend on which genes
#'   are kept, as they would after subsetting \code{x}. When
#'   \code{top_var_genes} is smaller than their number, the most variable
#'   of these genes are kept, in the order of their variances. Networks of
#'   several cell types built
#'   with the same \code{genes} have the same transcripts, which
#'   \code{\link{run_comparative_network}} needs, unless one of the genes
#'   is constant in a cell type and is left out there (see Details). A name
#'   that is not a row name of \code{x}, or that is the name of more than
#'   one row, is an error. Default \code{NULL} (all genes).
#'
#' @return A \code{\linkS4class{CelltypeNetworkResult}} object.
#'
#' @details
#' The workflow proceeds as follows:
#' \enumerate{
#'   \item Pseudobulk aggregation for the specified cell type,
#'     reusing existing \pkg{exposomeSC} infrastructure.
#'   \item Transformation of the pseudobulk counts, chosen by
#'     \code{transform}: DESeq2's variance-stabilising transformation or
#'     log-CPM. The transform is recorded in
#'     \code{metadata(result)$transform}.
#'   \item Constant features: a transcript or metabolite that takes one
#'     value in every donor (for example a gene with no counts in the
#'     cell type) is left out, with a message that names up to five
#'     transcripts and five metabolites; all are recorded
#'     in \code{metadata(result)$constant_features}; the call stops when
#'     every transcript or every metabolite is constant. Networks of
#'     different cell types can therefore differ in their features, which
#'     \code{\link{run_comparative_network}} does not accept. For a
#'     comparison, name the genes with \code{genes} and leave a gene that
#'     is constant in one of the cell types out of \code{genes} for all of
#'     them. \code{top_var_genes} on its own selects different genes in
#'     different cell types.
#'   \item Adjustment: if \code{exposure} or
#'     \code{covariates} are specified, regresses them out
#'     from both transcriptomic and metabolomic matrices.
#'     The resulting network captures the exposure-adjusted
#'     conditional dependency structure.
#'   \item Column-binding of the transcript values and the metabolite
#'     values (the residuals of both, where something was regressed out),
#'     each column scaled to unit variance.
#'   \item Block-calibrated graphical lasso estimation.
#'     With \code{coglasso}: collaborative graphical lasso
#'     with separate within/between-omic penalties and
#'     XStARS stability selection. With \code{block_glasso}: the
#'     graphical lasso (\code{glasso::glasso()}) on the correlation
#'     matrix, with the fixed penalties \code{lambda_w} within and
#'     \code{lambda_b} between omic layers, the diagonal unpenalised
#'     and no stability selection; it needs the \pkg{glasso}
#'     package.
#'   \item Returns a \code{CelltypeNetworkResult} object.
#' }
#'
#' @references
#' Albanese A, Kohlen W, Behrouzi P (2026). Multi-omics network
#'   reconstruction with collaborative graphical lasso.
#'   \emph{Bioinformatics} 42(7):btag477.
#'   \doi{10.1093/bioinformatics/btag477}
#'
#' Friedman J, Hastie T, Tibshirani R (2008). Sparse inverse
#'   covariance estimation with the graphical lasso.
#'   \emph{Biostatistics} 9(3):432--441.
#'   \doi{10.1093/biostatistics/kxm045}
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
#' if (requireNamespace("glasso", quietly = TRUE)) {
#'     network <- run_celltype_network(
#'         scee,
#'         metabolites,
#'         celltype = "Monocyte",
#'         method = "block_glasso",
#'         stability = FALSE,
#'         min_cells = 5L,
#'         top_var_genes = 3L,
#'         transform = "log_cpm"
#'     )
#'     network
#' }
run_celltype_network <- function(x, metabolites = NULL, celltype,
                                  celltype_col = "cell_type",
                                  exposure = NULL,
                                  covariates = NULL,
                                  method = c("coglasso",
                                             "block_glasso"),
                                  stability = TRUE,
                                  nlambda_w = 15L,
                                  nlambda_b = 15L,
                                  subsample_ratio = 0.8,
                                  rep_num = 20L,
                                  min_cells = 10L,
                                  top_var_genes = NULL,
                                  lambda_w = 0.3,
                                  lambda_b = 0.5,
                                  precomputed_network = NULL,
                                  transform = c("vst", "log_cpm"),
                                  genes = NULL) {

    stopifnot(is(x, "SingleCellExposomeExperiment"))
    celltype <- .single_name(celltype, "celltype")

    ## --- Handle precomputed network (named list) ---
    if (!is.null(precomputed_network)) {
        return(.wrap_precomputed(precomputed_network, celltype))
    }

    if (!is.null(exposure)) exposure <- .single_name(exposure, "exposure")
    if (!is.null(covariates)) covariates <- as.character(covariates)
    method <- match.arg(method)
    transform <- match.arg(transform)
    if (!is.matrix(metabolites)) {
        stop("'metabolites' must be a matrix with one row for each donor. ",
            "It can be left out only when 'precomputed_network' is given.",
            call. = FALSE)
    }
    if (method == "coglasso") {
        if (!is.numeric(subsample_ratio) || length(subsample_ratio) != 1L ||
            is.na(subsample_ratio) || subsample_ratio <= 0 ||
            subsample_ratio >= 1)
            stop("'subsample_ratio' must be a single number strictly ",
                 "between 0 and 1.")
        if (!is.numeric(rep_num) || length(rep_num) != 1L ||
            is.na(rep_num) || rep_num < 2 || rep_num != round(rep_num))
            stop("'rep_num' must be a single whole number of at least 2.")
        ## asked for before the pseudobulk and the transform are computed
        .require_coglasso("method")
    }

    ## --- Assemble cross-omic data ---
    assembled <- .assemble_crossomic(
        x, metabolites, celltype,
        celltype_col = celltype_col,
        exposure     = exposure,
        covariates   = covariates,
        min_cells    = min_cells,
        top_var_genes = top_var_genes,
        transform    = transform,
        genes        = genes)

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
        ## Build the penalty path once. coglasso::bs() would also run its
        ## own selection (xestars by default), which the XStARS step below
        ## would then repeat and overwrite.
        cg_path <- coglasso::coglasso(
            X,
            p = assembled$n_transcripts,
            nlambda_w = nlambda_w,
            nlambda_b = nlambda_b,
            verbose = FALSE
        )

        ## Select the penalties by XStARS, with the subsampling requested
        cg_stars <- coglasso::xstars(
            cg_path,
            stars_thresh = 0.1,
            stars_subsample_ratio = subsample_ratio,
            rep_num = as.integer(rep_num),
            verbose = FALSE
        )

        ## Extract precision and adjacency
        ## coglasso >= 1.1.0 uses sel_icov / sel_adj
        precision <- as.matrix(cg_stars$sel_icov)
        adjacency <- as.matrix(cg_stars$sel_adj)
        diag(adjacency) <- 0L

        ## Edge-level selection frequencies across the XStARS subsamples
        ## at the selected penalties (a symmetric feature-by-feature
        ## matrix in [0, 1]). sel_variability is a single number describing
        ## the whole selected model, so it goes to the metadata instead.
        stab_mat <- if (stability && !is.null(cg_stars$merge)) {
            as.matrix(cg_stars$merge)
        } else {
            matrix(nrow = 0, ncol = 0)
        }
        selection <- list(
            subsample_ratio   = subsample_ratio,
            rep_num           = as.integer(rep_num),
            selected_lambda_w = cg_stars$sel_lambda_w,
            selected_lambda_b = cg_stars$sel_lambda_b,
            selected_c        = cg_stars$sel_c,
            stars_variability = cg_stars$sel_variability
        )

    } else {
        ## block_glasso
        S <- cov(X)
        bg <- .block_glasso(S, blocks,
                             lambda_w = lambda_w,
                             lambda_b = lambda_b)
        precision <- bg$wi
        adjacency <- bg$adj
        stab_mat <- matrix(nrow = 0, ncol = 0)
        selection <- list(lambda_w = lambda_w, lambda_b = lambda_b)
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
            constant_features = assembled$constant_features,
            nlambda_w      = nlambda_w,
            nlambda_b      = nlambda_b,
            selection      = selection
        )
    )
}


#' Network among exposure-associated features
#'
#' Pre-selects the transcripts and metabolites associated with an exposure,
#' by stability selection or univariate testing, and estimates a cross-omic
#' network among the selected features only. Cheng et al. (2024) analysed
#' transcriptomic and metabolomic signatures of air pollution in two such
#' stages: features associated with the exposures were identified first,
#' and conditional multi-omic networks estimated afterwards.
#'
#' @param x \code{\linkS4class{SingleCellExposomeExperiment}}. Each
#'   cell's donor is read from the \code{colData} column named by
#'   \code{\link{sampleColumn}(x)}.
#' @param metabolites Numeric matrix, or a data frame of numeric columns;
#'   rows = donors, cols = metabolites. Row names must be donor IDs. A
#'   matrix without columns, or a missing or infinite value for an analysed
#'   donor, is an error, and so is a column with the name of another
#'   column or of a candidate gene: two nodes cannot have the same name.
#'   Leave the exposure itself out of the columns: a column that holds the
#'   values of the exposure, or a monotone transform of them, has a rank
#'   correlation of 1 or -1 with it, is always selected and becomes a node.
#' @param celltype Character; cell type to analyse.
#' @param exposure Character; exposure variable name, a column of
#'   \code{exposureData(x)}; another name is an error.
#' @param celltype_col Character; cell type column.
#'   Default \code{"cell_type"}.
#' @param covariates Character vector or NULL; columns of
#'   \code{exposureData(x)}. The selection is made given the covariates
#'   (see Details), and they are regressed out of the transcript and
#'   metabolite values before the network is estimated; another name is an
#'   error, and so is the exposure itself. Donors with a missing covariate
#'   value are left out, with a message. A covariate enters the selection
#'   by its ranks and the adjustment of the values as the number it is, so
#'   either way it is taken as ordered; give a covariate with more than two
#'   unordered categories as one 0/1 column for each category but one.
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
#'   \eqn{(2\pi_{thr} - 1)p \ge 1}). Here \eqn{p} is the number of
#'   candidate features that vary across the donors with a value of the
#'   exposure.
#' @param min_cells Integer; minimum cells per pseudobulk.
#' @param network_method Character; \code{"coglasso"} (default) or
#'   \code{"block_glasso"}, the estimator for the network on the
#'   selected features. \code{"coglasso"} needs at least one selected
#'   transcript and metabolite; when only one omic layer is selected the
#'   block graphical lasso is used instead, with a warning. Otherwise it
#'   requires the \pkg{coglasso} package, version 1.1.0 or later (an error
#'   when it is not installed), and chooses the penalties from a 10 by 10
#'   grid with the defaults of \code{coglasso::xstars()}. The estimator used
#'   is recorded in \code{metadata(result)$network_method}.
#' @param lambda_w Numeric; graphical-lasso penalty for pairs of
#'   features in the same omic layer, on the correlation scale, for the
#'   block graphical lasso (\code{network_method = "block_glasso"}, or the
#'   fallback when only one omic layer is selected). Default 0.3.
#' @param lambda_b Numeric; graphical-lasso penalty for pairs of
#'   features in different omic layers, on the correlation scale, for the
#'   block graphical lasso. Default 0.5.
#' @param transform Character; the transform of the pseudobulk counts,
#'   \code{"vst"} (default) or \code{"log_cpm"}, as in
#'   \code{\link{run_celltype_network}}. The call stops when \code{"vst"}
#'   cannot be computed; it does not switch to the other transform.
#' @param genes Character vector or \code{NULL}; row names of \code{x}, the
#'   genes that are candidates for the selection, taken after the transform
#'   as in \code{\link{run_celltype_network}}. Fewer candidates leave fewer
#'   tests to adjust for and a smaller \code{stability_q}. A name that is
#'   not a row name of \code{x}, or that is the name of more than one row,
#'   is an error. Default \code{NULL} (all genes).
#'
#' @return A \code{\linkS4class{CelltypeNetworkResult}}
#'   restricted to exposure-associated features. Its metadata record
#'   the expression transform, the covariates, the donors of the network
#'   (\code{donors}) and the number of them with a value of the exposure
#'   (\code{n_donors_exposure}), on whom the selection was made, the
#'   selection settings, for stability selection \code{pfer_bound}, in
#'   \code{constant_features} the candidates that are constant across the
#'   donors with a value of the exposure and so could not be selected, and
#'   in \code{selection} the penalties of the network:
#'   \code{lambda_w} and \code{lambda_b} for the block graphical lasso, or
#'   the size of the grid and the \code{selected_lambda_w},
#'   \code{selected_lambda_b} and \code{selected_c} that coglasso
#'   selected, with \code{stars_variability}.
#'
#' @details
#' Stage 1: Feature selection, among the donors with a value of the
#' exposure. Each transcript and metabolite is screened by its Spearman
#' correlation with the exposure. With \code{covariates} it is the partial
#' rank correlation: the ranks of the covariates are regressed out of the
#' ranks of the feature and of the exposure, and the residuals are
#' correlated.
#' A feature without variation across these donors, such as a gene without
#' a count in the cell type, is neither tested nor ranked by either method,
#' and \eqn{p} below counts the features that vary.
#' \code{"univariate"} tests the correlation (asymptotic \emph{t}
#' approximation, with one degree of freedom fewer for each covariate that
#' is not constant or a combination of the others among these donors) with
#' BH correction over these features.
#' The approximation is poor when both the exposure and a feature take one
#' value in all but a few donors, as an exposure that only a few donors
#' have, or that is above its limit of detection in a few, and a gene with
#' a count in one or two donors: the p-value of such a gene is far too
#' small whenever its donors are among those few, and with many such genes
#' among the candidates some pass the BH correction together. For an
#' exposure of this kind, name the candidate genes with \code{genes} and
#' leave out those with a count in few donors.
#' \code{"stability"} applies stability selection (Meinshausen and
#' Buhlmann 2010) to a marginal screen: in each of 100 random
#' subsamples of half the donors with an observed exposure, the
#' \code{stability_q} features with the largest absolute
#' correlation with the exposure are selected, and features selected
#' in a proportion of at least \code{pi_threshold} of subsamples are
#' retained. The screen remains defined when features outnumber
#' donors. Under the exchangeability and better-than-random-guessing
#' assumptions of Meinshausen and Buhlmann, the expected number of
#' falsely selected features is at most
#' \eqn{q^2 / ((2\pi_{thr} - 1)p)}, reported as \code{pfer_bound}.
#' Candidates with a count in very few donors are hardly ever among the
#' \code{stability_q} highest of a subsample, which the exchangeability
#' assumption does not allow for; with many of them among the candidates
#' the bound is too small. \code{genes} can leave them out.
#' At least three features must be selected.
#'
#' Transcript values are log-CPM or VST values and therefore relative
#' to each donor's library size. When the exposure strongly changes
#' genes that carry a large share of the counts, the remaining genes
#' shift in the opposite direction and can pass the screen; with few
#' measured genes this can dominate the univariate screen.
#'
#' Stage 2: Network estimation on the selected features only, with
#' collaborative graphical lasso model selection by XStARS
#' (\code{"coglasso"}) or the block graphical lasso
#' (\code{"block_glasso"}). It is one network for all the donors analysed.
#' The exposure is not regressed out of the features and is not a node, so
#' an edge can reflect the association that two features share with the
#' exposure, and the network does not estimate how the dependence between
#' features changes with the exposure.
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
#' if (requireNamespace("glasso", quietly = TRUE)) {
#'     net <- run_exposure_network(scee, metabolites, celltype = "Mono",
#'         exposure = "PM2.5",
#'         selection_method = "univariate",
#'         network_method = "block_glasso",
#'         transform = "log_cpm")
#'     net
#' }
run_exposure_network <- function(x, metabolites, celltype,
                                  exposure,
                                  celltype_col = "cell_type",
                                  covariates = NULL,
                                  selection_method = c(
                                      "stability",
                                      "univariate"),
                                  fdr_threshold = 0.05,
                                  pi_threshold = 0.6,
                                  stability_q = NULL,
                                  min_cells = 10L,
                                  network_method = c("coglasso",
                                                     "block_glasso"),
                                  lambda_w = 0.3,
                                  lambda_b = 0.5,
                                  transform = c("vst", "log_cpm"),
                                  genes = NULL) {

    stopifnot(is(x, "SingleCellExposomeExperiment"))
    celltype <- .single_name(celltype, "celltype")
    exposure <- .single_name(exposure, "exposure")
    selection_method <- match.arg(selection_method)
    network_method <- match.arg(network_method)
    transform <- match.arg(transform)
    if (!is.null(covariates)) covariates <- as.character(covariates)
    .check_adjustment_names(x, c(exposure, covariates))
    ## nothing of the exposure is left to test once it is adjusted for
    if (exposure %in% covariates) {
        stop("The exposure must not also be listed as a covariate.",
            call. = FALSE)
    }

    ## The covariates are regressed out of the features for the network of
    ## stage 2; the exposure is not (we want its signal)
    assembled <- .assemble_crossomic(
        x, metabolites, celltype,
        celltype_col = celltype_col,
        exposure     = NULL,
        covariates   = covariates,
        min_cells    = min_cells,
        transform    = transform,
        ## constant features stay in the matrix as constant columns, which
        ## the selection stage skips
        drop_constant = FALSE,
        genes        = genes,
        keep_unadjusted = TRUE)

    X <- assembled$data_matrix
    donors <- assembled$donors

    ## --- Stage 1: Feature selection ---
    ## among the donors with a value of the exposure, on the ranks of the
    ## values as they are before the adjustment: the covariates enter as
    ## ranks too (.partial_rank_correlation())
    exp_data <- exposureData(x)
    observed <- which(is.finite(exp_data[donors, exposure]))
    y <- exp_data[donors[observed], exposure]
    values <- assembled$unadjusted[observed, , drop = FALSE]
    adjustment <- if (length(covariates)) {
        as.matrix(exp_data[donors[observed], covariates, drop = FALSE])
    }
    n_covariates <- length(covariates)
    ## Nothing is left to select by when the exposure takes one value in
    ## these donors, or when the ranks of the covariates give its ranks
    if (length(y) > 1L) {
        exposure_ranks <- cbind(rank(y))
        if (n_covariates > 0L) {
            exposure_ranks <- .rank_residuals(exposure_ranks,
                apply(adjustment, 2L, rank))
        }
        if (!isTRUE(stats::sd(exposure_ranks[, 1L]) > 0)) {
            stop("The exposure does not vary across the donors with a ",
                "value of it", if (n_covariates > 0L)
                    " once the covariates are adjusted for",
                ", so nothing can be selected by it.", call. = FALSE)
        }
    }
    p <- ncol(X)
    selected <- rep(FALSE, p)

    pfer_bound <- NA_real_
    ## A feature that is constant across these donors, such as a gene
    ## without a count in the cell type, is neither tested nor ranked. It
    ## cannot be selected, so it does not count among the candidates of the
    ## bound of stability selection or of the default stability_q.
    varies <- apply(values, 2L, function(v) isTRUE(stats::sd(v) > 0))
    p_varying <- sum(varies)
    if (p_varying < 1L) {
        stop("No candidate feature varies across the donors with a value ",
            "of the exposure.", call. = FALSE)
    }

    if (selection_method == "univariate") {
        ## Spearman test per feature, given the covariates; constant
        ## features are not tested and not selected
        if (length(observed) < n_covariates + 4L) {
            stop("The test of the selection needs at least ",
                n_covariates + 4L, " donors with a value of the exposure; ",
                length(observed), " have one.", call. = FALSE)
        }
        pvals <- rep(NA_real_, p)
        pvals[varies] <- .partial_rank_correlation(
            values[, varies, drop = FALSE], y, adjustment)$p_value
        padj <- p.adjust(pvals, method = "BH")
        selected <- !is.na(padj) & padj < fdr_threshold

    } else {
        ## Stability selection (Meinshausen and Buhlmann 2010) over a
        ## marginal screen, which stays defined when p exceeds n
        if (!is.numeric(pi_threshold) || length(pi_threshold) != 1L ||
            is.na(pi_threshold) || pi_threshold <= 0.5 ||
            pi_threshold > 1)
            stop("pi_threshold must be a single number in (0.5, 1].")
        n_half <- floor(length(observed) / 2)
        if (n_half < max(5L, n_covariates + 4L))
            stop("Stability selection needs at least ",
                 2L * max(5L, n_covariates + 4L), " donors with an ",
                 "observed exposure; use selection_method = ",
                 "'univariate'.")
        if (is.null(stability_q))
            stability_q <- max(1,
                floor(sqrt((2 * pi_threshold - 1) * p_varying)))
        if (!is.numeric(stability_q) || length(stability_q) != 1L ||
            is.na(stability_q) || stability_q < 1 ||
            stability_q > p_varying || stability_q != round(stability_q))
            stop("stability_q must be a single integer between 1 and ",
                 "the number of features that vary across the donors (",
                 p_varying, ").")
        stability_q <- as.integer(stability_q)

        n_subsamples <- 100L
        sel_freq <- numeric(p)
        for (b in seq_len(n_subsamples)) {
            idx <- sample.int(length(observed), n_half)
            sel_freq <- sel_freq + .top_marginal_features(
                values[idx, , drop = FALSE], y[idx], stability_q,
                adjustment[idx, , drop = FALSE])
        }

        sel_prob <- sel_freq / n_subsamples
        selected <- sel_prob >= pi_threshold
        pfer_bound <- stability_q^2 / ((2 * pi_threshold - 1) * p_varying)
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
    if (use_coglasso) .require_coglasso("network_method")
    if (network_method == "coglasso" && !use_coglasso) {
        warning(
            "Only one omic layer was selected for cell type '", celltype,
            "', so the block graphical lasso is used instead of coglasso.",
            call. = FALSE
        )
    }
    if (use_coglasso) {
        cg_path <- coglasso::coglasso(
            sel_X,
            p = n_tx_sel,
            nlambda_w = 10L,
            nlambda_b = 10L,
            verbose = FALSE
        )
        cg_stars <- coglasso::xstars(cg_path, verbose = FALSE)
        precision <- as.matrix(cg_stars$sel_icov)
        adjacency <- as.matrix(cg_stars$sel_adj)
        diag(adjacency) <- 0L
        stab_mat <- matrix(nrow = 0, ncol = 0)
        selection <- list(nlambda_w = 10L, nlambda_b = 10L,
            selected_lambda_w = cg_stars$sel_lambda_w,
            selected_lambda_b = cg_stars$sel_lambda_b,
            selected_c = cg_stars$sel_c,
            stars_variability = cg_stars$sel_variability)
    } else {
        bg <- .block_glasso(S, sel_blocks, lambda_w = lambda_w,
            lambda_b = lambda_b)
        precision <- bg$wi
        adjacency <- bg$adj
        stab_mat <- matrix(nrow = 0, ncol = 0)
        selection <- list(lambda_w = lambda_w, lambda_b = lambda_b)
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
        method           = "exposure_restricted",
        metadata         = list(
            n_donors          = length(donors),
            n_donors_exposure = length(observed),
            donors            = donors,
            n_features        = n_selected,
            n_transcripts     = n_tx_sel,
            n_metabolites     = n_met_sel,
            exposure          = exposure,
            covariates        = covariates,
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
            selection         = selection,
            transform         = assembled$transform,
            n_total_features  = p,
            constant_features = assembled$node_info$feature[!varies],
            edge_counts       = edge_counts
        )
    )
}

# coglasso from version 1.1.0: xstars() and the argument 'p' of coglasso()
# are new in that version
.require_coglasso <- function(argument) {
    if (!.is_installed("coglasso") ||
            .installed_version("coglasso") < "1.1.0") {
        stop("Package 'coglasso' (>= 1.1.0) is required for ", argument,
            " = \"coglasso\". Install it with install.packages('coglasso'), ",
            "or use ", argument, " = \"block_glasso\".", call. = FALSE)
    }
    invisible(TRUE)
}

# ============================================================
# Wrap a precomputed network into CelltypeNetworkResult
# ============================================================

#' @keywords internal
.wrap_precomputed <- function(net, celltype) {
    if (!is.list(net))
        stop("precomputed_network must be a list with $adjacency.",
             call. = FALSE)
    if (is.null(net$adjacency))
        stop("precomputed_network list must have '$adjacency'.",
             call. = FALSE)
    adj  <- as.matrix(net$adjacency)
    if (nrow(adj) != ncol(adj)) {
        stop("'adjacency' of precomputed_network must be square.",
            call. = FALSE)
    }
    if (is.logical(adj)) {
        storage.mode(adj) <- "integer"
    }
    ## edges are between distinct features: no self-loops
    diag(adj) <- 0L
    prec <- if (!is.null(net$precision))
        as.matrix(net$precision) else matrix(NA_real_, nrow(adj), ncol(adj))
    if (!identical(dim(prec), dim(adj))) {
        stop("'precision' of precomputed_network must have the dimensions ",
            "of 'adjacency'.", call. = FALSE)
    }
    stab <- if (!is.null(net$stability))
        as.matrix(net$stability) else matrix(nrow = 0, ncol = 0)
    fnames <- if (!is.null(net$feature_names))
        net$feature_names else colnames(adj)
    if (is.null(fnames)) fnames <- sprintf("F%d", seq_len(ncol(adj)))
    fnames <- as.character(fnames)
    ## Block info
    if (!is.null(net$block_sizes)) {
        sizes <- unlist(net$block_sizes)
        if (!is.numeric(sizes) || !all(is.finite(sizes)) ||
                any(sizes < 0) || any(sizes != round(sizes)) ||
                is.null(names(sizes)) || anyNA(names(sizes)) ||
                !all(nzchar(names(sizes)))) {
            stop("'block_sizes' of precomputed_network must be a vector of ",
                "non-negative whole numbers named by omic layer.",
                call. = FALSE)
        }
        blk <- rep(names(sizes), sizes)
    } else if (!is.null(net$feature_blocks)) {
        blk <- net$feature_blocks
    } else {
        blk <- rep("unknown", length(fnames))
    }
    blk <- as.character(blk)
    if (length(fnames) != ncol(adj) || anyNA(fnames) ||
            !all(nzchar(fnames)) || anyDuplicated(fnames)) {
        stop("'feature_names' of precomputed_network, or the column names ",
            "of 'adjacency', must give one unique, non-missing and ",
            "non-empty name per row of 'adjacency'.", call. = FALSE)
    }
    if (length(blk) != length(fnames) || anyNA(blk) || !all(nzchar(blk))) {
        stop("'block_sizes' or 'feature_blocks' of precomputed_network must ",
            "label every feature, with no missing or empty labels.",
            call. = FALSE)
    }

    p <- length(fnames)
    dimnames(adj)  <- list(fnames, fnames)
    dimnames(prec) <- list(fnames, fnames)
    if (nrow(stab) == p && ncol(stab) == p)
        dimnames(stab) <- list(fnames, fnames)

    ## Same node table as for an estimated network: the omic layer of each
    ## feature and an integer block index, in order of first appearance
    node_info <- S4Vectors::DataFrame(
        feature    = fnames,
        omic_layer = as.character(blk),
        block      = match(blk, unique(blk))
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
            source     = "list"
        )
    )
}
