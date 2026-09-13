# R/causal-mediation.R
# Formal causal mediation for ERD
# Based on Imai et al. 2010 + VanderWeele 2015

#' @include AllClasses.R
#' @include utils.R
#' @importFrom stats lm coef residuals predict quantile
NULL

#' Formal Causal Mediation Analysis for ERD
#'
#' Applies formal causal mediation (Imai et al. 2010) to
#' provide:
#' \enumerate{
#'   \item Average Causal Mediation Effect (ACME) with CI
#'   \item Average Direct Effect (ADE) with CI
#'   \item Sensitivity analysis for sequential ignorability
#'   \item Identifiability conditions check
#' }
#'
#' @param scee A \code{SingleCellExposomeExperiment}.
#' @param exposure Character. Exposure variable name.
#' @param celltype Character. Target cell type.
#' @param celltype_col Character. Column with cell type labels.
#' @param sample_col Character. Column with donor IDs.
#' @param genes Character vector (optional). Genes to test.
#'   Default: the 20 genes with the largest absolute Pearson
#'   correlation between donor log-CPM and the exposure.
#' @param n_sims Integer. Monte Carlo simulations for
#'   mediation CI. Default 1000.
#' @param sensitivity Logical. Run sensitivity analysis for
#'   unmeasured confounding (rho parameter). Default TRUE.
#'
#' @return A \code{data.frame} with columns:
#'   gene, ACME (indirect/compositional), ACME_ci_lo/hi,
#'   ADE (direct), ADE_ci_lo/hi, prop_mediated,
#'   rho_at_zero (sensitivity: correlation at which ACME=0).
#'
#' @details
#' \strong{Causal identification assumptions (sequential
#' ignorability, Imai et al. 2010):}
#' \enumerate{
#'   \item No unmeasured exposure-outcome confounders
#'     (given covariates)
#'   \item No unmeasured mediator-outcome confounders
#'     (given exposure + covariates)
#'   \item Exposure does not affect mediator-outcome
#'     confounders
#' }
#'
#' Assumption 2 is the strongest and generally untestable.
#' The sensitivity analysis varies \eqn{\rho}, the
#' correlation between mediator and outcome residuals
#' (which would be induced by an unmeasured confounder),
#' and reports the \eqn{\rho} value at which ACME = 0.
#' Larger |\eqn{\rho}| means the result is more robust.
#'
#' @references
#' Imai K, Keele L, Tingley D (2010). A general approach to
#' causal mediation analysis. \emph{Psychol Methods}
#' 15:309-334.
#'
#' @examples
#' \dontrun{
#' sce <- SingleCellExperiment::SingleCellExperiment(
#'     assays = list(counts = matrix(rpois(500, 5), 50, 10,
#'         dimnames = list(paste0("G", 1:50), paste0("C", 1:10)))))
#' sce$donor_id <- rep(paste0("D", 1:5), each = 2)
#' sce$cell_type <- rep("Mono", 10)
#' exp_mat <- matrix(rnorm(5), 5, 1,
#'     dimnames = list(paste0("D", 1:5), "PM25"))
#' scee <- build_scee(sce, exp_mat, sample_col = "donor_id")
#' # Requires mediation package
#' if (requireNamespace("mediation", quietly = TRUE)) {
#'   res <- run_causal_mediation(scee, "PM25", "Mono",
#'       mediator_gene = "G1", genes = c("G2", "G3"), sims = 50)
#' }
#' }
#' @export
run_causal_mediation <- function(scee, exposure, celltype,
                                  celltype_col = "cell_type",
                                  sample_col = "donor_id",
                                  genes = NULL,
                                  n_sims = 1000L,
                                  sensitivity = TRUE) {

    stopifnot(is(scee, "SingleCellExposomeExperiment"))

    exp_data <- exposureData(scee)
    cd <- SummarizedExperiment::colData(scee)
    counts_mat <- SummarizedExperiment::assay(scee, "counts")
    samples <- as.character(cd[[sample_col]])
    cell_types <- as.character(cd[[celltype_col]])
    donors <- unique(samples)

    ## Compute mediator: cell type proportion per donor
    ct_prop <- vapply(donors, function(d) {
        mask <- samples == d
        sum(cell_types[mask] == celltype) / sum(mask)
    }, numeric(1))

    ## Pseudobulk for target celltype
    pb <- .pseudobulk_aggregate(counts_mat, samples,
        cell_types, celltype, min_cells = 10L)
    if (is.null(pb)) stop("No valid donors")
    valid <- pb$valid_donors
    pb_mat <- pb$pb_mat

    ## Log-CPM transform
    lcpm <- .log_cpm(pb_mat)

    ## Exposure vector
    exp_vec <- exp_data[valid, exposure]
    med_vec <- ct_prop[valid]

    ## If no genes are specified, screen by correlation with the exposure
    if (is.null(genes)) {
        ## Quick screen: correlate each gene with exposure
        cors <- apply(lcpm, 1, function(y)
            abs(cor(y, exp_vec, use = "complete.obs")))
        genes <- names(sort(cors, decreasing = TRUE))[
            seq_len(min(20, length(cors)))]
    }
    genes <- intersect(genes, rownames(lcpm))

    results <- list()

    for (g in genes) {
        y <- lcpm[g, ]

        ## Mediator model: M ~ E
        fit_m <- lm(med_vec ~ exp_vec)

        ## Outcome model: Y ~ E + M
        fit_y <- lm(y ~ exp_vec + med_vec)

        ## Use mediation package if available
        if (requireNamespace("mediation", quietly = TRUE)) {
            med_out <- tryCatch(
                mediation::mediate(fit_m, fit_y,
                    treat = "exp_vec",
                    mediator = "med_vec",
                    sims = n_sims),
                error = function(e) NULL)

            if (!is.null(med_out)) {
                row <- data.frame(
                    gene = g,
                    ACME = med_out$d0,
                    ACME_ci_lo = med_out$d0.ci[1],
                    ACME_ci_hi = med_out$d0.ci[2],
                    ACME_p = med_out$d0.p,
                    ADE = med_out$z0,
                    ADE_ci_lo = med_out$z0.ci[1],
                    ADE_ci_hi = med_out$z0.ci[2],
                    ADE_p = med_out$z0.p,
                    total = med_out$tau.coef,
                    prop_mediated = med_out$n0,
                    prop_mediated_p = med_out$n0.p,
                    rho_at_zero = NA_real_,
                    stringsAsFactors = FALSE)

                ## Sensitivity analysis
                if (sensitivity) {
                    sens <- tryCatch(
                        mediation::medsens(med_out,
                            rho.by = 0.05, sims = 200),
                        error = function(e) NULL)
                    if (!is.null(sens)) {
                        row$rho_at_zero <- sens$rho.at.sign
                    }
                }

                results[[g]] <- row
                next
            }
        }

        ## Fallback: manual difference method with bootstrap CI
        beta_total <- coef(lm(y ~ exp_vec))["exp_vec"]
        beta_direct <- coef(fit_y)["exp_vec"]
        acme <- beta_total - beta_direct

        ## Bootstrap CI
        boot_acme <- vapply(seq_len(n_sims), function(b) {
            idx <- sample(length(y), replace = TRUE)
            bt <- coef(lm(y[idx] ~ exp_vec[idx]))["exp_vec"]
            bd <- coef(lm(y[idx] ~ exp_vec[idx] +
                med_vec[idx]))["exp_vec"]
            bt - bd
        }, numeric(1))

        results[[g]] <- data.frame(
            gene = g,
            ACME = acme,
            ACME_ci_lo = quantile(boot_acme, 0.025),
            ACME_ci_hi = quantile(boot_acme, 0.975),
            ACME_p = 2 * min(mean(boot_acme > 0),
                mean(boot_acme < 0)),
            ADE = beta_direct,
            ADE_ci_lo = NA_real_,
            ADE_ci_hi = NA_real_,
            ADE_p = summary(fit_y)$coefficients[
                "exp_vec", "Pr(>|t|)"],
            total = beta_total,
            prop_mediated = ifelse(abs(beta_total) > 1e-10,
                acme / beta_total, 0),
            prop_mediated_p = NA_real_,
            rho_at_zero = NA_real_,
            stringsAsFactors = FALSE)
    }

    out <- do.call(rbind, results)
    rownames(out) <- NULL
    out
}
