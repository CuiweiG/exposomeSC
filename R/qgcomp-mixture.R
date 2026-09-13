# R/qgcomp-mixture.R
# Formal mixture analysis via quantile g-computation

#' @include AllClasses.R
#' @include utils.R
#' @importFrom stats as.formula
NULL

#' Formal mixture analysis via quantile g-computation
#'
#' Wraps \code{qgcomp::qgcomp.noboot()} or
#' \code{qgcomp::qgcomp.boot()} on pseudobulk data for
#' cell-type-specific exposure mixture analysis. Unlike
#' \code{\link{run_sc_mixture}} (simplified screening), it
#' reports a confidence interval and p-value for the overall
#' mixture effect and, without bootstrapping, the positive and
#' negative exposure weights estimated by \pkg{qgcomp}.
#'
#' @param x A \code{\linkS4class{SingleCellExposomeExperiment}}.
#' @param exposures Character vector. Exposure variable names.
#' @param celltype Character. Target cell type.
#' @param celltype_col Character. Column with cell type labels.
#' @param sample_col Character. Column with donor IDs.
#' @param target_genes Character vector (optional). Genes to
#'   summarise. Default: top 20 most variable.
#' @param covariates Character vector (optional). Columns of
#'   \code{exposureData}; an unknown name is an error.
#' @param q Integer. Number of quantile bins. Default 4.
#' @param bootstrap Logical. If \code{TRUE}, use
#'   \code{qgcomp.boot()} for bootstrap confidence intervals.
#'   Default \code{FALSE} (faster; uses asymptotic CI).
#' @param B Integer. Number of bootstrap replicates (only if
#'   \code{bootstrap = TRUE}). Default 200.
#' @param min_cells Integer. Default 10.
#'
#' @return A list with components:
#' \describe{
#'   \item{positive_weights}{Named numeric. Exposures that
#'     increase the response, with relative weights;
#'     \code{NULL} when \code{bootstrap = TRUE}, because
#'     \code{qgcomp.boot()} does not estimate weights.}
#'   \item{negative_weights}{Named numeric. Exposures that
#'     decrease the response; \code{NULL} when
#'     \code{bootstrap = TRUE}.}
#'   \item{mixture_effect}{Overall mixture effect estimate
#'     (\eqn{\psi}).}
#'   \item{mixture_ci}{95\% confidence interval for \eqn{\psi}.}
#'   \item{mixture_pvalue}{P-value for \eqn{\psi}.}
#'   \item{fit}{The \code{qgcomp} model object (for
#'     plotting and further analysis).}
#'   \item{method}{Character: \code{"qgcomp.noboot"} or
#'     \code{"qgcomp.boot"}.}
#'   \item{celltype, n_donors, n_genes}{Analysis metadata.}
#' }
#'
#' @details
#' Quantile g-computation (Keil et al. 2020) estimates the
#' overall effect of jointly increasing all exposures by one
#' quantile, and without bootstrapping decomposes it into
#' positive and negative contributions. The weights are
#' interpretable only when the exposure effects are of a
#' similar shape (Keil et al. 2020).
#'
#' The pseudobulk response is the mean log-CPM of the target
#' genes over donors with complete exposure and covariate data.
#'
#' \strong{Comparison with \code{run_sc_mixture}:}
#' \tabular{lll}{
#'   Feature \tab run_sc_mixture \tab run_mixture_qgcomp \cr
#'   CI for overall effect \tab No \tab Yes (asymptotic or bootstrap) \cr
#'   Pos/neg weights \tab No \tab Yes (without bootstrap) \cr
#'   P-value for overall effect \tab No \tab Yes \cr
#'   Speed \tab Fast \tab Moderate \cr
#'   Minimum n \tab 5 \tab 10+ recommended \cr
#' }
#'
#' @references
#' Keil AP et al. (2020). A quantile-based g-computation
#' approach to addressing the effects of exposure mixtures.
#' \emph{Environ Health Perspect} 128:047004.
#' \doi{10.1289/EHP5838}
#'
#' @export
#' @examples
#' # Requires the qgcomp package
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
#' \donttest{
#' mix <- run_mixture_qgcomp(scee,
#'     exposures = c("E1", "E2", "E3"),
#'     celltype = "Mono", celltype_col = "cell_type",
#'     sample_col = "donor_id", min_cells = 3L)
#' mix$positive_weights
#' mix$negative_weights
#' }
run_mixture_qgcomp <- function(x, exposures, celltype,
                                celltype_col = "cell_type",
                                sample_col = "donor_id",
                                target_genes = NULL,
                                covariates = NULL,
                                q = 4L, bootstrap = FALSE,
                                B = 200L, min_cells = 10L) {
    stopifnot(is(x, "SingleCellExposomeExperiment"))

    if (!requireNamespace("qgcomp", quietly = TRUE))
        stop("Package 'qgcomp' required. ",
             "install.packages('qgcomp')")

    exp_data <- slot(x, "exposureData")
    missing_exp <- setdiff(exposures, colnames(exp_data))
    if (length(missing_exp) > 0)
        stop("Exposures not found: ",
             paste(missing_exp, collapse = ", "))

    cd <- SummarizedExperiment::colData(x)
    counts_mat <- SummarizedExperiment::assay(x, "counts")
    samples <- as.character(cd[[sample_col]])
    cell_types <- as.character(cd[[celltype_col]])

    pb_result <- .pseudobulk_aggregate(
        counts_mat, samples, cell_types, celltype,
        min_cells = min_cells)
    if (is.null(pb_result))
        stop("No valid donors for '", celltype, "'")

    covariates <- as.character(covariates)
    missing_cov <- setdiff(covariates, colnames(exp_data))
    if (length(missing_cov) > 0)
        stop("Covariates not found: ",
             paste(missing_cov, collapse = ", "))

    valid <- pb_result$valid_donors
    valid <- valid[stats::complete.cases(
        exp_data[valid, c(exposures, covariates), drop = FALSE])]
    if (length(valid) < 5L)
        stop("Need >= 5 donors with complete exposure data. Found: ",
             length(valid))

    log_cpm <- .log_cpm(pb_result$pb_mat)

    if (is.null(target_genes)) {
        gv <- apply(log_cpm, 1, stats::var)
        n_top <- min(20L, nrow(log_cpm))
        target_genes <- names(sort(gv,
            decreasing = TRUE))[seq_len(n_top)]
    }
    target_genes <- intersect(target_genes, rownames(log_cpm))
    if (length(target_genes) == 0)
        stop("No target genes found")

    ## Response: mean log-CPM across target genes
    y <- colMeans(log_cpm[target_genes, valid, drop = FALSE])

    ## Build data frame
    df <- data.frame(y = y)
    for (exp_name in exposures) {
        df[[exp_name]] <- exp_data[valid, exp_name]
    }
    for (cov in covariates) {
        df[[cov]] <- exp_data[valid, cov]
    }

    ## Build formula
    fml <- stats::reformulate(c(exposures, covariates), response = "y")

    ## Run qgcomp
    fit <- tryCatch({
        if (bootstrap) {
            qgcomp::qgcomp.boot(fml, data = df,
                expnms = exposures, q = q, B = B,
                family = stats::gaussian())
        } else {
            qgcomp::qgcomp.noboot(fml, data = df,
                expnms = exposures, q = q,
                family = stats::gaussian())
        }
    }, error = function(e) {
        stop("qgcomp failed: ", conditionMessage(e))
    })

    ## Extract results directly from qgcomp object
    mixture_effect <- fit$psi
    mixture_se <- sqrt(fit$var.psi)
    mixture_ci <- as.numeric(fit$ci)
    ## fit$pval follows fit$coef, whose first element is the intercept
    psi_index <- match("psi1", names(fit$coef))
    mixture_pval <- if (!is.na(psi_index) &&
                        length(fit$pval) >= psi_index) {
        fit$pval[psi_index]
    } else {
        2 * stats::pnorm(abs(mixture_effect / mixture_se),
            lower.tail = FALSE)
    }

    list(
        positive_weights = fit$pos.weights,
        negative_weights = fit$neg.weights,
        mixture_effect = mixture_effect,
        mixture_ci = mixture_ci,
        mixture_pvalue = mixture_pval,
        fit = fit,
        method = if (bootstrap) "qgcomp.boot" else
            "qgcomp.noboot",
        celltype = celltype,
        n_donors = length(valid),
        n_genes = length(target_genes)
    )
}
