# R/state-coupling.R
# Cell-State Coupling
# Exposure-driven correlation-trajectory analysis

#' @importFrom stats cor quantile lm pchisq
NULL

#' Cell-State Coupling: Correlation-Trajectory Analysis
#'
#' Tests whether gene-protein coupling varies along a
#' continuous cell state trajectory AND whether this
#' variation depends on a donor-level exposure.
#'
#' @param scee A \code{SingleCellExposomeExperiment} with
#'   an altExp containing protein data.
#' @param gene Character. Gene name.
#' @param protein Character. Protein name.
#' @param exposure Character. Exposure variable.
#' @param celltype Character. Cell type.
#' @param state_col Character. Column in \code{colData}
#'   containing continuous cell state (pseudotime, PC1,
#'   or any continuous trajectory coordinate).
#' @param celltype_col,sample_col Character. Column names.
#' @param altexp_name Character. altExp name. Default "CITE".
#' @param n_bins Integer. Number of quantile bins of the state
#'   variable over all cells of the cell type. Default 5. Tied
#'   quantiles give fewer bins.
#' @param min_cells_per_bin Integer. Minimum cells per bin.
#'   Default 20.
#' @param min_donors Integer. Minimum donors. Default 10.
#'
#' @return A \code{data.frame} with one row per state bin:
#'   bin, mean_state, beta0 (baseline coupling at this state),
#'   beta1 (exposure effect at this state), se_beta1, p_beta1,
#'   n_donors, mean_r. With at least two bins, the attribute
#'   \code{slope_heterogeneity} holds the test of whether
#'   \code{beta1} differs between bins (statistic, df1, df2,
#'   pvalue, method). \code{method} names the test that was
#'   used: the F-test, the likelihood-ratio test, or
#'   \code{"not estimable"} with \code{NA} values when neither
#'   model could be fitted.
#'
#' @details
#' \strong{Mathematical model:}
#'
#' For each state bin \eqn{t}:
#' \enumerate{
#'   \item Within each donor \eqn{d}, compute Spearman
#'     correlation \eqn{r_d(t)} between gene and protein
#'     for cells in bin \eqn{t}
#'   \item Fisher z-transform: \eqn{z_d(t) = \text{arctanh}(r_d(t))}
#'   \item Meta-regression: \eqn{z_d(t) = \beta_0(t) +
#'     \beta_1(t) \cdot E_d + \epsilon_d}
#' }
#'
#' \eqn{\beta_1(t)} is the exposure effect on coupling
#' at state \eqn{t}.
#'
#' \strong{Slope heterogeneity.} The same donors contribute to
#' every bin, so the per-bin estimates are not independent. All
#' donor-by-bin correlations therefore enter one multilevel
#' meta-regression (\code{metafor::rma.mv}) with bin-specific
#' intercepts and exposure slopes, a random donor effect shared
#' across bins and a residual heterogeneity term. The F-test of
#' the bin-by-exposure terms asks whether the exposure effect on
#' coupling differs between bins. A non-significant result does
#' not show that the effect is the same in every state.
#'
#' When the correlations carry little signal the two variance
#' components are estimated at zero. The fit then tries
#' \pkg{metafor}'s optimisers in turn and keeps the first that
#' converges with finite coefficients. If the resulting
#' coefficient covariance is still unusable for the F-test, the
#' same question is put as a likelihood-ratio test between the
#' full and the additive model, both fitted by maximum
#' likelihood, and \code{method} says so.
#'
#' @examples
#' set.seed(1)
#' donors <- sprintf("D%02d", 1:12)
#' donor <- rep(donors, each = 60)
#' gene <- matrix(stats::rpois(2 * length(donor), 10), nrow = 2,
#'     dimnames = list(c("Gene1", "Gene2"), paste0("c", seq_along(donor))))
#' protein <- matrix(stats::rpois(2 * length(donor), 5), nrow = 2,
#'     dimnames = list(c("Prot1", "Prot2"), colnames(gene)))
#' sce <- SingleCellExperiment::SingleCellExperiment(
#'     assays = list(counts = gene),
#'     colData = S4Vectors::DataFrame(donor = donor, celltype = "T",
#'         pseudotime = stats::runif(length(donor))))
#' SingleCellExperiment::altExp(sce, "CITE") <-
#'     SummarizedExperiment::SummarizedExperiment(
#'         assays = list(counts = protein))
#' exp_mat <- matrix(seq(0, 2, length.out = 12), ncol = 1,
#'     dimnames = list(donors, "exposure"))
#' scee <- build_scee(sce, exp_mat, sample_col = "donor")
#' if (requireNamespace("metafor", quietly = TRUE)) {
#'     run_state_coupling(scee, gene = "Gene1", protein = "Prot1",
#'         exposure = "exposure", celltype = "T",
#'         state_col = "pseudotime", n_bins = 3L,
#'         min_cells_per_bin = 10L, min_donors = 10L)
#' }
#' @export
run_state_coupling <- function(scee, gene, protein, exposure,
                                celltype,
                                state_col,
                                celltype_col = "celltype",
                                sample_col = "donor",
                                altexp_name = "CITE",
                                n_bins = 5L,
                                min_cells_per_bin = 20L,
                                min_donors = 10L) {

    stopifnot(is(scee, "SingleCellExposomeExperiment"))

    if (!requireNamespace("metafor", quietly = TRUE))
        stop("Package 'metafor' required")

    exp_data <- exposureData(scee)
    if (!exposure %in% colnames(exp_data))
        stop("Exposure '", exposure, "' not in exposureData")
    exp_vec <- setNames(exp_data[, exposure],
        rownames(exp_data))

    cd <- SummarizedExperiment::colData(scee)
    for (column in c(celltype_col, sample_col, state_col)) {
        if (!column %in% colnames(cd))
            stop("Column '", column, "' not in colData")
    }
    if (!altexp_name %in% SingleCellExperiment::altExpNames(scee))
        stop("altExp '", altexp_name, "' not found")
    prot_se <- SingleCellExperiment::altExp(scee, altexp_name)
    if (!gene %in% rownames(scee))
        stop("Gene '", gene, "' not found")
    if (!protein %in% rownames(prot_se))
        stop("Protein '", protein, "' not found in altExp")

    ct_idx <- which(cd[[celltype_col]] == celltype)
    if (length(ct_idx) == 0)
        stop("No cells for ", celltype)

    ## Get state variable
    state <- as.numeric(cd[[state_col]][ct_idx])
    donors <- as.character(cd[[sample_col]][ct_idx])

    ## Get gene + protein data
    gene_vals <- as.numeric(
        SummarizedExperiment::assay(scee, "counts")[gene, ct_idx])
    prot_vals <- as.numeric(
        SummarizedExperiment::assay(prot_se, "counts")[protein, ct_idx])

    ## Quantile bins of the state over all cells of the cell type; tied
    ## quantiles give fewer bins
    bin_breaks <- unique(quantile(state, probs = seq(0, 1,
        length.out = n_bins + 1), na.rm = TRUE))
    if (length(bin_breaks) < 2L)
        stop("The state variable does not vary within '", celltype, "'")
    bin_labels <- seq_len(length(bin_breaks) - 1L)
    bins <- cut(state, breaks = bin_breaks, labels = bin_labels,
        include.lowest = TRUE)

    ## For each bin: per-donor correlation + meta-regression
    results <- list()
    long_rows <- list()

    for (b in bin_labels) {
        b_idx <- which(bins == b)
        if (length(b_idx) < min_cells_per_bin * 2) next

        b_donors <- unique(donors[b_idx])
        donor_z <- numeric()
        donor_v <- numeric()
        donor_exp <- numeric()
        donor_id <- character()

        for (d in b_donors) {
            d_idx <- b_idx[donors[b_idx] == d]
            n_d <- length(d_idx)
            if (n_d < min_cells_per_bin) next

            g_d <- gene_vals[d_idx]
            p_d <- prot_vals[d_idx]
            if (sd(g_d) < 1e-10 || sd(p_d) < 1e-10) next

            r <- cor(g_d, p_d, method = "spearman")
            r <- max(min(r, 0.999), -0.999)
            z <- atanh(r)
            v <- 1.06 / (n_d - 3)

            e_d <- exp_vec[d]
            if (is.na(e_d)) next

            donor_z <- c(donor_z, z)
            donor_v <- c(donor_v, v)
            donor_exp <- c(donor_exp, e_d)
            donor_id <- c(donor_id, d)
        }

        if (length(donor_z) < min_donors) next

        ## Meta-regression
        fit <- tryCatch(
            metafor::rma(yi = donor_z, vi = donor_v,
                mods = ~ donor_exp,
                method = "REML", test = "knha"),
            error = function(e) NULL)

        if (is.null(fit)) next

        results[[as.character(b)]] <- data.frame(
            bin = as.integer(b),
            mean_state = mean(state[b_idx], na.rm = TRUE),
            beta0 = fit$beta[1, 1],
            beta1 = fit$beta[2, 1],
            se_beta1 = fit$se[2],
            p_beta1 = fit$pval[2],
            n_donors = length(donor_z),
            mean_r = mean(tanh(donor_z)),
            stringsAsFactors = FALSE)
        long_rows[[as.character(b)]] <- data.frame(
            bin = as.integer(b), donor = donor_id, yi = donor_z,
            vi = donor_v, exposure = unname(donor_exp),
            stringsAsFactors = FALSE)
    }

    if (length(results) == 0) return(data.frame())

    out <- do.call(rbind, results)
    rownames(out) <- NULL

    ## Do the exposure slopes differ between bins? The same donors
    ## contribute to every bin, so one multilevel model is fitted
    if (nrow(out) >= 2L) {
        long <- do.call(rbind, long_rows)
        rownames(long) <- NULL
        long$bin <- factor(long$bin)
        long$observation <- seq_len(nrow(long))
        attr(out, "slope_heterogeneity") <- .slope_heterogeneity_test(long)
    }

    out
}

## Fit the multilevel meta-regression behind the slope-heterogeneity test,
## trying metafor's optimisers in turn. When the donor-by-bin correlations
## carry little signal the donor and residual variance components are
## estimated at zero, and whether an optimiser stops there cleanly depends
## on the platform: the same data converged on Linux and Windows and failed
## on macOS. metafor signals an optimiser that did not converge as an
## error, so the search moves to the next optimiser; a fit is kept only
## when every coefficient is finite.
.fit_rma_mv <- function(mods, data, method) {
    controls <- list(
        list(optimizer = "nlminb"),
        list(optimizer = "optim", optmethod = "BFGS"),
        list(optimizer = "optim", optmethod = "Nelder-Mead"))
    for (ctrl in controls) {
        fit <- tryCatch(
            metafor::rma.mv(yi = yi, V = vi, mods = mods,
                random = list(~ 1 | donor, ~ 1 | observation),
                data = data, method = method, test = "t",
                control = ctrl),
            error = function(e) NULL)
        if (!is.null(fit) && all(is.finite(fit$beta)))
            return(fit)
    }
    NULL
}

## The test of whether the exposure slope differs between bins. The Wald
## F-test of the bin-by-exposure terms is used when the fitted model gives
## it a usable coefficient covariance. A variance component estimated at
## zero can leave that covariance singular, in which case the same
## question is put as a likelihood-ratio test between the full and the
## additive model, both fitted by maximum likelihood, which does not need
## the covariance. The method field of the result says which was used.
.slope_heterogeneity_test <- function(long) {
    not_estimable <- data.frame(statistic = NA_real_, df1 = NA_real_,
                                df2 = NA_real_, pvalue = NA_real_,
                                method = "not estimable",
                                stringsAsFactors = FALSE)

    joint <- .fit_rma_mv(~ bin * exposure, long, "REML")
    if (!is.null(joint)) {
        interaction_terms <- grep(":exposure$", rownames(joint$beta))
        wald <- tryCatch(
            stats::anova(joint, btt = interaction_terms),
            error = function(e) NULL)
        if (!is.null(wald) && is.finite(wald$QMp)) {
            return(data.frame(statistic = wald$QM,
                              df1 = wald$QMdf[1],
                              df2 = wald$QMdf[2],
                              pvalue = wald$QMp,
                              method = paste("rma.mv bin-by-exposure",
                                             "F-test with a random",
                                             "donor effect"),
                              stringsAsFactors = FALSE))
        }
    }

    full <- .fit_rma_mv(~ bin * exposure, long, "ML")
    reduced <- .fit_rma_mv(~ bin + exposure, long, "ML")
    if (is.null(full) || is.null(reduced)) return(not_estimable)
    lrt <- tryCatch(stats::anova(full, reduced), error = function(e) NULL)
    if (is.null(lrt) || !is.finite(lrt$pval)) return(not_estimable)
    data.frame(statistic = lrt$LRT,
               df1 = lrt$parms.f - lrt$parms.r,
               df2 = NA_real_,
               pvalue = lrt$pval,
               method = paste("rma.mv bin-by-exposure likelihood-ratio",
                              "test (ML) with a random donor effect"),
               stringsAsFactors = FALSE)
}
