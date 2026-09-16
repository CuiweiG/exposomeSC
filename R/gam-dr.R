# R/gam-dr.R
# GAM-based dose-response: nonparametric alternative

#' @include AllClasses.R
#' @include utils.R
#' @importFrom stats AIC BIC p.adjust lm
NULL

#' GAM-based dose-response analysis
#'
#' Fits generalized additive models (GAMs) for cell-type-
#' specific dose-response curves, allowing the data to
#' determine the shape of the exposure-expression
#' relationship without pre-specifying a polynomial degree.
#' Unlike \code{\link{run_dose_response}}, the shape is not
#' restricted to a polynomial.
#'
#' @param x A \code{\linkS4class{SingleCellExposomeExperiment}}.
#' @param exposure Character. Exposure variable name.
#' @param celltype Character. Target cell type.
#' @param celltype_col Character. Column with cell type labels.
#' @param sample_col Character. Column with donor IDs.
#' @param target_genes Character vector (optional). Genes to
#'   test. Default: top 50 most variable.
#' @param k Integer. Maximum basis dimension for the smooth
#'   term. Default \code{min(n_donors - 1, 10)}. Higher k
#'   allows more wiggliness; the REML smoothing penalty
#'   limits overfitting.
#' @param min_cells Integer. Min cells per donor. Default 10.
#' @param loocv Logical. If \code{TRUE}, compute leave-one-out
#'   cross-validated R² for the linear and GAM models.
#'   Default \code{TRUE}.
#'
#' @return A \code{DataFrame} with columns: gene, edf
#'   (effective degrees of freedom for smooth term),
#'   p_smooth (significance of smooth term vs intercept),
#'   p_nonlinear (approximate test that a fully penalised smooth
#'   added to a linear exposure term is zero, that is, of
#'   departure from linearity; see Details),
#'   AIC_linear, AIC_gam, R2_linear, R2_gam,
#'   R2_loocv_linear (if loocv=TRUE), R2_loocv_gam,
#'   deviance_explained, celltype, n_donors.
#'
#' @details
#' For each gene, two models are fit:
#' \enumerate{
#'   \item Linear: \code{y ~ exposure} (OLS)
#'   \item GAM: \code{y ~ s(exposure, k = k)} (penalized
#'     thin-plate regression spline via \code{mgcv::gam})
#' }
#'
#' The GAM smooth is estimated using restricted maximum
#' likelihood (REML), which selects the smoothness penalty.
#' The effective degrees of freedom (edf) indicate the
#' complexity: edf close to 1 means approximately linear, and
#' larger values mean more curvature.
#'
#' \code{p_nonlinear} comes from a third model,
#' \code{y ~ exposure + s(exposure, bs = "tp", m = c(2, 0))},
#' whose smooth has no unpenalised null space, so the linear
#' trend is carried by the parametric term and the smooth can
#' only represent departures from linearity. Its p-value is the
#' approximate Wald-type test for smooth terms of Wood (2013).
#' Comparing the linear and GAM fits by an F-test on their
#' residual deviances is not valid here: the GAM's effective
#' degrees of freedom are estimated and often close to one.
#'
#' Leave-one-out cross-validated R² (loocv) provides an
#' honest estimate of out-of-sample prediction accuracy,
#' unlike apparent R² which inflates with model complexity.
#' This is critical for small-sample dose-response where
#' overfitting is a real risk.
#'
#' \strong{Why GAM over polynomial:}
#' \itemize{
#'   \item Polynomials impose a global parametric form
#'   \item GAMs adapt locally -- can capture threshold effects,
#'     plateaus, U-shapes
#'   \item Penalised smoothness selection limits overfitting
#'   \item edf provides an interpretable complexity measure
#' }
#'
#' @references
#' Wood SN (2017). Generalized Additive Models: An
#' Introduction with R. 2nd ed. Chapman and Hall/CRC.
#'
#' Wood SN (2013). On p-values for smooth components of an
#' extended generalized additive model. \emph{Biometrika}
#' 100:221-228. \doi{10.1093/biomet/ass048}
#'
#' Vandenberg LN et al. (2012). Hormones and endocrine-
#' disrupting chemicals: low-dose effects and nonmonotonic
#' dose responses. \emph{Endocr Rev} 33:378-455.
#' \doi{10.1210/er.2011-1050}
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
#' exp_mat <- matrix(rnorm(20), nrow = 10,
#'     dimnames = list(paste0("D", 1:10), c("E1", "E2")))
#' scee <- build_scee(sce, exp_mat, sample_col = "donor_id")
#' if (requireNamespace("mgcv", quietly = TRUE)) {
#'     dr <- run_dose_response_gam(scee, exposure = "E1",
#'         celltype = "Mono", celltype_col = "cell_type",
#'         sample_col = "donor_id", min_cells = 3L,
#'         target_genes = paste0("G", 1:3))
#'     head(dr)
#' }
run_dose_response_gam <- function(x, exposure, celltype,
                                    celltype_col = "cell_type",
                                    sample_col = "donor_id",
                                    target_genes = NULL,
                                    k = NULL,
                                    min_cells = 10L,
                                    loocv = TRUE) {
    stopifnot(is(x, "SingleCellExposomeExperiment"))

    if (!requireNamespace("mgcv", quietly = TRUE))
        stop("Package 'mgcv' required. install.packages('mgcv')")

    exp_data <- slot(x, "exposureData")
    if (!exposure %in% colnames(exp_data))
        stop("Exposure '", exposure, "' not found")

    cd <- SummarizedExperiment::colData(x)
    counts_mat <- SummarizedExperiment::assay(x, "counts")
    samples <- as.character(cd[[sample_col]])
    cell_types <- as.character(cd[[celltype_col]])

    pb_result <- .pseudobulk_aggregate(
        counts_mat, samples, cell_types, celltype,
        min_cells = min_cells)
    if (is.null(pb_result))
        stop("No valid donors for '", celltype, "'")

    valid <- pb_result$valid_donors
    n <- length(valid)

    if (n < 6L)
        stop("Need >= 6 donors for GAM. Found: ", n)

    if (is.null(k)) k <- min(n - 1L, 10L)
    k <- max(k, 3L)  # mgcv minimum

    log_cpm <- .log_cpm(pb_result$pb_mat)

    if (is.null(target_genes)) {
        gv <- apply(log_cpm, 1, stats::var)
        n_top <- min(50L, nrow(log_cpm))
        target_genes <- names(sort(gv,
            decreasing = TRUE))[seq_len(n_top)]
    }
    target_genes <- intersect(target_genes, rownames(log_cpm))
    if (length(target_genes) == 0)
        stop("No target genes found")

    dose <- exp_data[valid, exposure]

    results <- lapply(target_genes, function(gene) {
        y <- log_cpm[gene, valid]
        df <- data.frame(y = y, dose = dose)

        ## Linear model
        fit_lin <- tryCatch(lm(y ~ dose, data = df),
            error = function(e) NULL)
        if (is.null(fit_lin)) return(NULL)

        ## GAM
        fit_gam <- tryCatch(
            mgcv::gam(y ~ s(dose, k = k, bs = "tp"),
                data = df, method = "REML"),
            error = function(e) NULL)
        if (is.null(fit_gam)) return(NULL)

        ## Extract summary
        s_gam <- summary(fit_gam)
        edf <- s_gam$s.table[1, "edf"]
        p_smooth <- s_gam$s.table[1, "p-value"]

        ## Departure from linearity: linear term plus a smooth whose
        ## null space is penalised away (m = c(2, 0))
        p_nonlinear <- tryCatch({
            fit_departure <- mgcv::gam(
                y ~ dose + s(dose, k = k, bs = "tp", m = c(2, 0)),
                data = df, method = "REML")
            summary(fit_departure)$s.table[1, "p-value"]
        }, error = function(e) NA_real_)

        ## R² values
        r2_lin <- summary(fit_lin)$r.squared
        r2_gam <- s_gam$r.sq

        ## LOOCV R²
        r2_loocv_lin <- NA_real_
        r2_loocv_gam <- NA_real_
        if (loocv && n >= 6L) {
            ## Leave-one-out for linear
            pred_loo_lin <- numeric(n)
            pred_loo_gam <- numeric(n)
            for (i in seq_len(n)) {
                df_train <- df[-i, ]
                df_test <- df[i, , drop = FALSE]

                f_lin <- tryCatch(
                    lm(y ~ dose, data = df_train),
                    error = function(e) NULL)
                if (!is.null(f_lin)) {
                    pred_loo_lin[i] <- predict(f_lin,
                        newdata = df_test)
                } else {
                    pred_loo_lin[i] <- NA_real_
                }

                f_gam <- tryCatch(
                    mgcv::gam(y ~ s(dose, k = min(k, n - 2L),
                        bs = "tp"),
                        data = df_train, method = "REML"),
                    error = function(e) NULL)
                if (!is.null(f_gam)) {
                    pred_loo_gam[i] <- predict(f_gam,
                        newdata = df_test)
                } else {
                    pred_loo_gam[i] <- NA_real_
                }
            }

            ss_total <- sum((y - mean(y))^2)
            if (ss_total > 0 && !any(is.na(pred_loo_lin))) {
                r2_loocv_lin <- 1 - sum(
                    (y - pred_loo_lin)^2) / ss_total
            }
            if (ss_total > 0 && !any(is.na(pred_loo_gam))) {
                r2_loocv_gam <- 1 - sum(
                    (y - pred_loo_gam)^2) / ss_total
            }
        }

        data.frame(
            gene = gene,
            edf = edf,
            p_smooth = p_smooth,
            p_nonlinear = p_nonlinear,
            AIC_linear = AIC(fit_lin),
            AIC_gam = AIC(fit_gam),
            R2_linear = r2_lin,
            R2_gam = r2_gam,
            R2_loocv_linear = r2_loocv_lin,
            R2_loocv_gam = r2_loocv_gam,
            deviance_explained = s_gam$dev.expl,
            celltype = celltype,
            n_donors = n,
            stringsAsFactors = FALSE)
    })

    out <- do.call(rbind, Filter(Negate(is.null), results))
    if (is.null(out) || nrow(out) == 0) {
        return(S4Vectors::DataFrame(
            gene = character(), edf = numeric(),
            p_smooth = numeric(), p_nonlinear = numeric(),
            AIC_linear = numeric(), AIC_gam = numeric(),
            R2_linear = numeric(), R2_gam = numeric(),
            R2_loocv_linear = numeric(),
            R2_loocv_gam = numeric(),
            deviance_explained = numeric(),
            celltype = character(), n_donors = integer()))
    }
    S4Vectors::DataFrame(out)
}
