# R/mixture.R
# Mixture exposure modelling at cell-type level
#
# NOTE: This implements a *simplified* quantile-scored linear
# mixture approach. It is NOT a full qgcomp or WQS
# implementation. For formal mixture methods, use the
# dedicated `qgcomp` or `gWQS` packages.

#' @include AllClasses.R
#' @include utils.R
#' @importFrom stats lm coef quantile var as.formula
NULL

#' Cell-type-specific mixture exposure analysis
#'
#' Estimates the relative contribution of multiple correlated
#' exposures to gene expression within a specific cell type
#' using a simplified quantile-scored linear model on
#' pseudobulk data.
#'
#' @param x A \code{\linkS4class{SingleCellExposomeExperiment}}.
#' @param exposures Character vector. Exposure variable names.
#' @param celltype Character. Which cell type to analyse.
#' @param celltype_col Character. Column with cell type labels.
#' @param sample_col Character. Column with donor IDs.
#' @param target_genes Character vector (optional). Genes to
#'   test. Default: top 20 most variable genes.
#' @param covariates Character vector (optional). Additional
#'   covariates from \code{exposureData}; a name that is not a
#'   column of \code{exposureData} is an error.
#' @param q Integer. Number of quantile bins. Default 4. Tied
#'   exposure values (for example at a limit of detection) can
#'   give fewer bins.
#' @param min_cells Integer. Minimum cells per donor. Default 10.
#'
#' @return A \code{list} with components:
#' \describe{
#'   \item{weights}{Named numeric vector. Relative contribution
#'     of each exposure (sums to 1): the mean over target genes
#'     of the absolute coefficients of the quantile-scored
#'     exposures, normalised.}
#'   \item{mixture_coef}{Overall mixture coefficient (sum of
#'     quantile-scored exposure coefficients, averaged across
#'     target genes).}
#'   \item{gene_results}{\code{DataFrame} with per-gene mixture
#'     effect estimates.}
#'   \item{method}{Character. Always \code{"quantile_linear"}.}
#'   \item{celltype}{Character. Cell type analysed.}
#'   \item{n_donors}{Integer. Donors with complete exposure and
#'     covariate data used in the analysis.}
#'   \item{n_genes}{Integer. Genes with a fitted model.}
#' }
#'
#' @details
#' This is a \strong{simplified screening method}, not a formal
#' quantile-based g-computation (Keil et al. 2020) or weighted
#' quantile sum regression (Carrico et al. 2015). The approach:
#'
#' \enumerate{
#'   \item Donors with a missing exposure or covariate are
#'     removed, and exposures are scored into quantile bins
#'     (1 to \code{q}; fewer when quantiles are tied).
#'   \item For each target gene, a linear model is fit on
#'     pseudobulk log-CPM with all quantile-scored exposures
#'     and covariates as predictors.
#'   \item Weights are the absolute exposure coefficients,
#'     averaged across target genes and normalised to sum to 1.
#'     Because every exposure is on the same quantile scale,
#'     the coefficients are directly comparable.
#' }
#'
#' \strong{Limitations:}
#' \itemize{
#'   \item No bootstrap confidence intervals on weights.
#'   \item Does not separate positive/negative weight
#'     directions as in formal qgcomp.
#'   \item Log-CPM may have heteroscedasticity at low counts.
#' }
#'
#' For rigorous mixture analyses suitable for publication,
#' consider \code{qgcomp::qgcomp.noboot()} or
#' \code{gWQS::gwqs()} applied to the pseudobulk data
#' returned by this package's internal aggregation.
#'
#' @references
#' Keil AP et al. (2020). A quantile-based g-computation
#' approach to addressing the effects of exposure mixtures.
#' \emph{Environ Health Perspect} 128:047004.
#' \doi{10.1289/EHP5838}
#'
#' Carrico C et al. (2015). Characterization of weighted
#' quantile sum regression for highly correlated data in a
#' risk analysis setting. \emph{J Agric Biol Environ Stat}
#' 20:100-120.
#'
#' @export
#' @examples
#' library(SingleCellExperiment); library(S4Vectors)
#' set.seed(1)
#' counts <- matrix(rpois(5000, 8), nrow = 50,
#'     dimnames = list(paste0("G", 1:50), paste0("c", 1:100)))
#' sce <- SingleCellExperiment(assays = list(counts = counts),
#'     colData = DataFrame(cell_id = paste0("c", 1:100),
#'         donor_id = rep(paste0("D", 1:10), each = 10),
#'         cell_type = rep(c("Mono", "NK"), 50)))
#' exp_mat <- matrix(rnorm(30), nrow = 10,
#'     dimnames = list(paste0("D", 1:10), c("E1", "E2", "E3")))
#' scee <- build_scee(sce, exp_mat, sample_col = "donor_id")
#' mix <- run_sc_mixture(scee, exposures = c("E1", "E2", "E3"),
#'     celltype = "Mono", celltype_col = "cell_type",
#'     sample_col = "donor_id", min_cells = 3L)
#' mix$weights
run_sc_mixture <- function(x, exposures, celltype,
                            celltype_col = "cell_type",
                            sample_col = "donor_id",
                            target_genes = NULL,
                            covariates = NULL,
                            q = 4L, min_cells = 10L) {
    stopifnot(is(x, "SingleCellExposomeExperiment"))

    exp_data <- slot(x, "exposureData")
    missing_exp <- setdiff(exposures, colnames(exp_data))
    if (length(missing_exp) > 0) {
        stop("Exposures not found: ",
             paste(missing_exp, collapse = ", "))
    }
    covariates <- as.character(covariates)
    missing_cov <- setdiff(covariates, colnames(exp_data))
    if (length(missing_cov) > 0) {
        stop("Covariates not found: ",
             paste(missing_cov, collapse = ", "))
    }

    cd <- SummarizedExperiment::colData(x)
    counts_mat <- SummarizedExperiment::assay(x, "counts")
    samples <- as.character(cd[[sample_col]])
    cell_types <- as.character(cd[[celltype_col]])

    ## Pseudobulk via shared utility
    pb_result <- .pseudobulk_aggregate(
        counts_mat, samples, cell_types, celltype,
        min_cells = min_cells)

    if (is.null(pb_result)) {
        stop("Cell type '", celltype, "' not found or no ",
             "donors with >= ", min_cells, " cells")
    }

    valid <- pb_result$valid_donors
    complete <- stats::complete.cases(
        exp_data[valid, c(exposures, covariates), drop = FALSE])
    valid <- valid[complete]
    pb <- pb_result$pb_mat

    if (length(valid) < 5L) {
        stop("Need >= 5 donors with ", min_cells,
             "+ cells and complete exposure data for '", celltype,
             "'. Found: ", length(valid))
    }

    ## Log-CPM via shared utility
    log_cpm <- .log_cpm(pb)

    if (is.null(target_genes)) {
        gene_vars <- apply(log_cpm[, valid, drop = FALSE], 1, stats::var)
        n_top <- min(20L, nrow(pb))
        target_genes <- names(sort(gene_vars,
                                    decreasing = TRUE))[
                                        seq_len(n_top)]
    }
    ## Keep only genes present in matrix
    target_genes <- intersect(target_genes, rownames(log_cpm))
    if (length(target_genes) == 0) {
        stop("No target genes found in expression matrix")
    }

    ## Quantile-score exposures on one shared scale; tied quantiles
    ## (for example at a detection limit) give fewer bins
    exp_sub <- exp_data[valid, exposures, drop = FALSE]
    exp_q <- vapply(exposures, function(name) {
        values <- exp_sub[, name]
        breaks <- unique(stats::quantile(
            values, probs = seq(0, 1, length.out = q + 1)))
        if (length(breaks) < 2L) {
            stop("Exposure '", name, "' does not vary among the ",
                 "analysed donors.")
        }
        as.integer(cut(values, breaks = breaks, include.lowest = TRUE))
    }, integer(length(valid)))
    exp_q <- matrix(exp_q, nrow = length(valid),
                    dimnames = list(valid, exposures))

    ## Build predictor data frame (quantile-scored exposures and
    ## covariates on their original scale)
    pred_df <- as.data.frame(exp_q)
    for (cov in covariates) {
        pred_df[[cov]] <- exp_data[valid, cov]
    }

    fml <- stats::reformulate(c(exposures, covariates), response = "y")

    ## Per-gene: fit quantile-scored linear model
    gene_fits <- lapply(target_genes, function(gene) {
        df <- data.frame(y = log_cpm[gene, valid], pred_df,
                         check.names = FALSE)

        fit <- tryCatch(lm(fml, data = df),
                         error = function(e) NULL)
        if (is.null(fit)) return(NULL)

        exp_coefs <- stats::coef(fit)[exposures]
        exp_coefs[is.na(exp_coefs)] <- 0
        list(gene = gene, coefficients = exp_coefs)
    })
    gene_fits <- Filter(Negate(is.null), gene_fits)

    if (length(gene_fits)) {
        coefficient_matrix <- do.call(rbind, lapply(gene_fits,
            function(fit) fit$coefficients))
        gene_results <- data.frame(
            gene = vapply(gene_fits, function(fit) fit$gene, character(1)),
            mixture_effect = rowSums(coefficient_matrix),
            stringsAsFactors = FALSE)
        mean_abs <- colMeans(abs(coefficient_matrix))
        weights <- if (sum(mean_abs) > 0) mean_abs / sum(mean_abs) else
            rep(1 / length(exposures), length(exposures))
        mixture_coef <- mean(gene_results$mixture_effect)
    } else {
        gene_results <- NULL
        weights <- rep(NA_real_, length(exposures))
        mixture_coef <- NA_real_
    }
    names(weights) <- exposures

    list(
        weights = weights,
        mixture_coef = mixture_coef,
        method = "quantile_linear",
        celltype = celltype,
        n_donors = length(valid),
        n_genes = length(gene_fits),
        gene_results = if (!is.null(gene_results))
            S4Vectors::DataFrame(gene_results) else
            S4Vectors::DataFrame()
    )
}
