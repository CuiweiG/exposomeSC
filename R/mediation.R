# R/mediation.R
# Causal mediation: Exposure -> Cell composition -> Gene expression

#' @include AllClasses.R
#' @include utils.R
#' @importFrom stats lm as.formula p.adjust
NULL

#' Mediation analysis: exposure through cell composition
#'
#' Tests whether an exposure affects gene expression
#' directly or indirectly through changes in cell type
#' composition. This causal decomposition addresses a
#' fundamental question in single-cell exposomics: does
#' the exposure alter gene expression within existing cell
#' types (direct effect), or does it shift cell type
#' proportions (mediated/indirect effect)?
#'
#' @param x A \code{\linkS4class{SingleCellExposomeExperiment}}.
#' @param exposure Character. Exposure variable name.
#' @param celltype_col Character. Column with cell type labels.
#' @param sample_col Character. Column with donor IDs.
#' @param mediator_celltype Character. Cell type whose
#'   proportion serves as the mediator.
#' @param outcome_celltype Character. Cell type for outcome
#'   gene expression.
#' @param target_genes Character vector (optional). Genes to
#'   test. Default: top 30 most variable.
#' @param covariates Character vector (optional).
#' @param n_sims Integer. Number of simulations for confidence
#'   intervals (Imai et al. method). Default 1000.
#' @param min_cells Integer. Default 10.
#'
#' @return A \code{DataFrame} with columns: gene,
#'   ACME (average causal mediation effect),
#'   ACME_ci_lo, ACME_ci_hi, ACME_p,
#'   ADE (average direct effect),
#'   ADE_ci_lo, ADE_ci_hi, ADE_p,
#'   total_effect, prop_mediated,
#'   mediator_celltype, outcome_celltype, n_donors.
#'
#' @details
#' The causal model is:
#' \preformatted{
#' Exposure -> Cell proportion (mediator) -> Gene expression
#'     |                                        ^
#'     +---------- direct effect ---------------+
#' }
#'
#' This uses the counterfactual framework of Imai, Keele &
#' Tingley (2010), implemented via the \code{mediation}
#' package:
#'
#' \enumerate{
#'   \item \strong{Mediator model}: cell type proportion ~
#'     exposure + covariates (linear regression)
#'   \item \strong{Outcome model}: pseudobulk log-CPM ~
#'     exposure + cell proportion + covariates (linear)
#'   \item \strong{Mediation}: estimate ACME (indirect),
#'     ADE (direct), and proportion mediated via
#'     quasi-Bayesian simulation
#' }
#'
#' \strong{Causal assumptions (sequential ignorability):}
#' \itemize{
#'   \item No unmeasured exposure-outcome confounders
#'   \item No unmeasured mediator-outcome confounders
#'   \item Exposure does not affect confounders of
#'     mediator-outcome
#' }
#'
#' These are strong assumptions. Report sensitivity analyses
#' (e.g., E-values) in publications.
#'
#' @references
#' Imai K, Keele L, Tingley D (2010). A general approach to
#' causal mediation analysis. \emph{Psychol Methods}
#' 15:309-334. \doi{10.1037/a0020761}
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
#'         cell_type = sample(c("Mono", "NK", "T"), 100, TRUE)))
#' exp_mat <- matrix(rnorm(20), nrow = 10,
#'     dimnames = list(paste0("D", 1:10), c("E1", "E2")))
#' scee <- build_scee(sce, exp_mat, sample_col = "donor_id")
#' \donttest{
#' med <- run_mediation(scee, exposure = "E1",
#'     celltype_col = "cell_type", sample_col = "donor_id",
#'     mediator_celltype = "Mono", outcome_celltype = "NK",
#'     min_cells = 2L)
#' head(med)
#' }
run_mediation <- function(x, exposure,
                           celltype_col = "cell_type",
                           sample_col = "donor_id",
                           mediator_celltype,
                           outcome_celltype,
                           target_genes = NULL,
                           covariates = NULL,
                           n_sims = 1000L,
                           min_cells = 10L) {
    stopifnot(is(x, "SingleCellExposomeExperiment"))

    if (!requireNamespace("mediation", quietly = TRUE))
        stop("Package 'mediation' required. ",
             "install.packages('mediation')")

    exp_data <- slot(x, "exposureData")
    if (!exposure %in% colnames(exp_data))
        stop("Exposure '", exposure, "' not found")

    cd <- SummarizedExperiment::colData(x)
    counts_mat <- SummarizedExperiment::assay(x, "counts")
    samples <- as.character(cd[[sample_col]])
    cell_types <- as.character(cd[[celltype_col]])
    donors <- unique(samples)

    ## Compute mediator: proportion of mediator_celltype
    ct_levels <- sort(unique(cell_types))
    if (!mediator_celltype %in% ct_levels)
        stop("mediator_celltype '", mediator_celltype,
             "' not found")
    if (!outcome_celltype %in% ct_levels)
        stop("outcome_celltype '", outcome_celltype,
             "' not found")

    ## Per-donor cell proportions
    prop_mediator <- vapply(donors, function(d) {
        mask <- samples == d
        if (sum(mask) < min_cells) return(NA_real_)
        sum(cell_types[mask] == mediator_celltype) / sum(mask)
    }, numeric(1))

    ## Pseudobulk for outcome cell type
    pb_result <- .pseudobulk_aggregate(
        counts_mat, samples, cell_types, outcome_celltype,
        min_cells = min_cells)
    if (is.null(pb_result))
        stop("No valid donors for '", outcome_celltype, "'")

    valid_pb <- pb_result$valid_donors
    valid <- intersect(donors[!is.na(prop_mediator)], valid_pb)
    if (length(valid) < 8L)
        stop("Need >= 8 donors for mediation. Found: ",
             length(valid))

    log_cpm <- .log_cpm(pb_result$pb_mat)

    if (is.null(target_genes)) {
        gv <- apply(log_cpm[, valid, drop = FALSE], 1,
            stats::var)
        n_top <- min(30L, nrow(log_cpm))
        target_genes <- names(sort(gv,
            decreasing = TRUE))[seq_len(n_top)]
    }
    target_genes <- intersect(target_genes, rownames(log_cpm))
    if (length(target_genes) == 0)
        stop("No target genes found")

    ## Build analysis data frame
    df <- data.frame(
        donor = valid,
        exposure = exp_data[valid, exposure],
        mediator = prop_mediator[valid],
        stringsAsFactors = FALSE)
    if (!is.null(covariates)) {
        for (cov in covariates) {
            if (cov %in% colnames(exp_data))
                df[[cov]] <- exp_data[valid, cov]
        }
    }

    cov_str <- if (!is.null(covariates) &&
                   any(covariates %in% colnames(df))) {
        paste("+", paste(intersect(covariates, colnames(df)),
            collapse = " + "))
    } else ""

    results <- lapply(target_genes, function(gene) {
        df$y <- log_cpm[gene, valid]

        ## Mediator model
        fml_med <- as.formula(paste("mediator ~ exposure",
            cov_str))
        ## Outcome model
        fml_out <- as.formula(paste("y ~ exposure + mediator",
            cov_str))

        fit_med <- tryCatch(lm(fml_med, data = df),
            error = function(e) NULL)
        fit_out <- tryCatch(lm(fml_out, data = df),
            error = function(e) NULL)

        if (is.null(fit_med) || is.null(fit_out)) return(NULL)

        ## Run mediation
        med <- tryCatch(
            mediation::mediate(fit_med, fit_out,
                treat = "exposure", mediator = "mediator",
                sims = n_sims),
            error = function(e) NULL)

        if (is.null(med)) return(NULL)

        data.frame(
            gene = gene,
            ACME = med$d0,         # avg causal mediation effect
            ACME_ci_lo = med$d0.ci[1],
            ACME_ci_hi = med$d0.ci[2],
            ACME_p = med$d0.p,
            ADE = med$z0,          # avg direct effect
            ADE_ci_lo = med$z0.ci[1],
            ADE_ci_hi = med$z0.ci[2],
            ADE_p = med$z0.p,
            total_effect = med$tau.coef,
            prop_mediated = med$n0,
            mediator_celltype = mediator_celltype,
            outcome_celltype = outcome_celltype,
            n_donors = length(valid),
            stringsAsFactors = FALSE)
    })

    out <- do.call(rbind, Filter(Negate(is.null), results))
    if (is.null(out) || nrow(out) == 0) {
        return(S4Vectors::DataFrame(
            gene = character(), ACME = numeric(),
            ACME_ci_lo = numeric(), ACME_ci_hi = numeric(),
            ACME_p = numeric(), ADE = numeric(),
            ADE_ci_lo = numeric(), ADE_ci_hi = numeric(),
            ADE_p = numeric(), total_effect = numeric(),
            prop_mediated = numeric(),
            mediator_celltype = character(),
            outcome_celltype = character(),
            n_donors = integer()))
    }

    out$ACME_padj <- p.adjust(out$ACME_p, "BH")
    out$ADE_padj <- p.adjust(out$ADE_p, "BH")
    S4Vectors::DataFrame(out)
}
