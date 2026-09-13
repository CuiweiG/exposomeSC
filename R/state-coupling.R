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
#' @param n_bins Integer. Number of state bins per donor.
#'   Default 5.
#' @param min_cells_per_bin Integer. Minimum cells per bin.
#'   Default 20.
#' @param min_donors Integer. Minimum donors. Default 10.
#'
#' @return A \code{data.frame} with one row per state bin:
#'   bin, mean_state, beta0 (baseline coupling at this state),
#'   beta1 (exposure effect at this state), p_beta1,
#'   n_donors, mean_r.
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
#' AT state \eqn{t}. A significant interaction between
#' state and exposure (\eqn{\beta_1} varies across bins)
#' indicates state-dependent rewiring.
#'
#' Global test: ANOVA on \eqn{\beta_1(t)} across bins.
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
    exp_vec <- setNames(exp_data[, exposure],
        rownames(exp_data))

    cd <- SummarizedExperiment::colData(scee)
    ct_idx <- which(cd[[celltype_col]] == celltype)
    if (length(ct_idx) == 0)
        stop("No cells for ", celltype)

    ## Get state variable
    state <- as.numeric(cd[[state_col]][ct_idx])
    donors <- as.character(cd[[sample_col]][ct_idx])

    ## Get gene + protein data
    gene_vals <- as.numeric(
        SummarizedExperiment::assay(scee, "counts")[gene, ct_idx])
    prot_se <- SingleCellExperiment::altExp(scee, altexp_name)
    prot_vals <- as.numeric(
        SummarizedExperiment::assay(prot_se, "counts")[protein, ct_idx])

    ## Create state bins (quantile-based)
    bin_breaks <- quantile(state, probs = seq(0, 1,
        length.out = n_bins + 1), na.rm = TRUE)
    bin_labels <- seq_len(n_bins)
    bins <- cut(state, breaks = bin_breaks, labels = bin_labels,
        include.lowest = TRUE)

    ## For each bin: per-donor correlation + meta-regression
    results <- list()

    for (b in bin_labels) {
        b_idx <- which(bins == b)
        if (length(b_idx) < min_cells_per_bin * 2) next

        b_donors <- unique(donors[b_idx])
        donor_z <- numeric()
        donor_v <- numeric()
        donor_exp <- numeric()

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
    }

    if (length(results) == 0) return(data.frame())

    out <- do.call(rbind, results)
    rownames(out) <- NULL

    ## Global test: is beta1 varying across bins?
    if (nrow(out) >= 3) {
        ## Cochran Q-like test on beta1 heterogeneity
        beta1s <- out$beta1
        se1s <- out$se_beta1
        w <- 1 / se1s^2
        beta1_avg <- sum(w * beta1s) / sum(w)
        Q <- sum(w * (beta1s - beta1_avg)^2)
        Q_df <- nrow(out) - 1
        Q_p <- 1 - pchisq(Q, Q_df)

        attr(out, "heterogeneity_Q") <- Q
        attr(out, "heterogeneity_p") <- Q_p
        attr(out, "interpretation") <-
            if (Q_p < 0.05)
                "Significant state-dependent rewiring"
            else
                "Coupling change is consistent across states"
    }

    out
}
