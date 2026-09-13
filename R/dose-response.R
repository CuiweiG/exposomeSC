# R/dose-response.R
# Nonlinear dose-response modelling at cell-type level

#' @include AllClasses.R
#' @include utils.R
#' @importFrom stats lm poly AIC BIC predict fitted
NULL

#' Cell-type-specific exposure dose-response analysis
#'
#' Fits linear and polynomial dose-response models for each
#' gene within a cell type, testing whether the
#' exposure-expression relationship is nonlinear. Most
#' environmental health studies assume linearity, but
#' endocrine disruptors and air pollutants often show
#' U-shaped or threshold effects (Vandenberg et al. 2012).
#'
#' @param x A \code{\linkS4class{SingleCellExposomeExperiment}}.
#' @param exposure Character. Exposure variable name.
#' @param celltype Character. Target cell type.
#' @param celltype_col Character. Column with cell type labels.
#' @param sample_col Character. Column with donor IDs.
#' @param target_genes Character vector (optional). Genes to
#'   test. Default: top 50 most variable genes.
#' @param max_degree Integer. Maximum polynomial degree to
#'   test. Default 3 (cubic). Higher values risk overfitting
#'   with small sample sizes.
#' @param min_cells Integer. Minimum cells per donor. Default 10.
#'
#' @return A \code{DataFrame} with columns: gene, best_model
#'   ("linear", "quadratic", "cubic"), linear_coef, AIC_linear,
#'   AIC_best, p_nonlinear (F-test vs linear), celltype,
#'   n_donors.
#'
#' @details
#' For each gene, the function fits:
#' \enumerate{
#'   \item \code{y ~ exposure} (linear)
#'   \item \code{y ~ poly(exposure, 2)} (quadratic)
#'   \item \code{y ~ poly(exposure, 3)} (cubic, if max_degree
#'     >= 3)
#' }
#' on pseudobulk log-CPM values. Model selection uses AIC.
#' A likelihood ratio test (via anova) compares the best
#' nonlinear model against the linear baseline.
#'
#' This addresses a critical gap: environmental epidemiology
#' frequently encounters nonlinear dose-response
#' relationships (Vandenberg et al. 2012 \emph{Endocr Rev}),
#' but existing single-cell tools only test linear effects.
#'
#' @references
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
#' dr <- run_dose_response(scee, exposure = "E1",
#'     celltype = "Mono", celltype_col = "cell_type",
#'     sample_col = "donor_id", min_cells = 3L,
#'     max_degree = 2L)
#' head(dr)
run_dose_response <- function(x, exposure, celltype,
                               celltype_col = "cell_type",
                               sample_col = "donor_id",
                               target_genes = NULL,
                               max_degree = 3L,
                               min_cells = 10L) {
    stopifnot(is(x, "SingleCellExposomeExperiment"))
    max_degree <- min(max_degree, 3L)

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
        stop("Cell type '", celltype, "' not found or no ",
             "donors with >= ", min_cells, " cells")

    valid <- pb_result$valid_donors
    pb <- pb_result$pb_mat

    ## Need enough donors for polynomial fitting
    min_needed <- max_degree + 2L
    if (length(valid) < min_needed)
        stop("Need >= ", min_needed, " donors for degree-",
             max_degree, " polynomial. Found: ",
             length(valid))

    log_cpm <- .log_cpm(pb)

    if (is.null(target_genes)) {
        gene_vars <- apply(log_cpm, 1, stats::var)
        n_top <- min(50L, nrow(pb))
        target_genes <- names(sort(gene_vars,
            decreasing = TRUE))[seq_len(n_top)]
    }
    target_genes <- intersect(target_genes, rownames(log_cpm))
    if (length(target_genes) == 0)
        stop("No target genes found in expression matrix")

    dose <- exp_data[valid, exposure]
    model_names <- c("linear", "quadratic", "cubic")[
        seq_len(max_degree)]

    results <- lapply(target_genes, function(gene) {
        y <- log_cpm[gene, valid]
        df <- data.frame(y = y, dose = dose)

        fits <- list()
        aics <- numeric()

        ## Linear
        fits[["linear"]] <- tryCatch(
            lm(y ~ dose, data = df),
            error = function(e) NULL)
        if (is.null(fits[["linear"]])) return(NULL)
        aics["linear"] <- AIC(fits[["linear"]])

        ## Quadratic
        if (max_degree >= 2L) {
            fits[["quadratic"]] <- tryCatch(
                lm(y ~ poly(dose, 2, raw = TRUE), data = df),
                error = function(e) NULL)
            if (!is.null(fits[["quadratic"]]))
                aics["quadratic"] <- AIC(fits[["quadratic"]])
        }

        ## Cubic
        if (max_degree >= 3L && length(valid) >= 6L) {
            fits[["cubic"]] <- tryCatch(
                lm(y ~ poly(dose, 3, raw = TRUE), data = df),
                error = function(e) NULL)
            if (!is.null(fits[["cubic"]]))
                aics["cubic"] <- AIC(fits[["cubic"]])
        }

        best <- names(which.min(aics))
        linear_coef <- coef(fits[["linear"]])["dose"]

        ## F-test: best nonlinear vs linear
        p_nonlinear <- NA_real_
        if (best != "linear" && !is.null(fits[[best]])) {
            f_test <- tryCatch(
                anova(fits[["linear"]], fits[[best]]),
                error = function(e) NULL)
            if (!is.null(f_test) && nrow(f_test) == 2)
                p_nonlinear <- f_test[2, "Pr(>F)"]
        }

        data.frame(
            gene = gene,
            best_model = best,
            linear_coef = as.numeric(linear_coef),
            AIC_linear = aics["linear"],
            AIC_best = min(aics),
            p_nonlinear = p_nonlinear,
            celltype = celltype,
            n_donors = length(valid),
            stringsAsFactors = FALSE)
    })

    out <- do.call(rbind, Filter(Negate(is.null), results))
    if (is.null(out) || nrow(out) == 0) {
        return(S4Vectors::DataFrame(
            gene = character(), best_model = character(),
            linear_coef = numeric(), AIC_linear = numeric(),
            AIC_best = numeric(), p_nonlinear = numeric(),
            celltype = character(), n_donors = integer()))
    }
    S4Vectors::DataFrame(out)
}
