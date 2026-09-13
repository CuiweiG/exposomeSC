# R/power.R
# Power analysis for sc-ExWAS study design

#' @include AllClasses.R
#' @importFrom stats rnorm rpois var power.t.test sd pnorm qnorm
NULL

#' Power analysis for sc-ExWAS study design
#'
#' Estimates statistical power for detecting an exposure
#' effect at the pseudobulk level, or computes the minimum
#' number of donors needed to achieve a target power. This
#' is essential for study design -- single-cell studies are
#' expensive, and researchers need to know whether their
#' planned sample size is adequate \emph{before} collecting
#' data.
#'
#' @param n_donors Integer vector. Donor sample sizes to
#'   evaluate. Default: \code{seq(10, 100, by = 10)}.
#' @param effect_size Numeric. Expected log2 fold change per
#'   unit exposure. Default 0.5.
#' @param n_genes Integer. Number of genes tested (for
#'   multiple testing correction). Default 5000.
#' @param n_celltypes Integer. Number of cell types tested.
#'   Default 5.
#' @param alpha Numeric. Significance level after BH
#'   correction. Default 0.05.
#' @param dispersion Numeric. DESeq2-like dispersion
#'   parameter. Default 0.1 (typical for pseudobulk with
#'   >= 20 cells per donor). Higher dispersion = more noise
#'   = less power.
#' @param exposure_sd Numeric. Standard deviation of the
#'   exposure variable across donors. Default 1
#'   (standardised).
#'
#' @return A \code{DataFrame} with columns: n_donors, power,
#'   effect_size, n_tests, alpha_adjusted, dispersion.
#'
#' @details
#' Power is computed analytically using a normal
#' approximation to the Wald test used by DESeq2. The
#' effective per-test alpha is Bonferroni-adjusted for the
#' total number of tests (\code{n_genes * n_celltypes}).
#' This is conservative relative to BH, so the true power
#' under BH will be somewhat higher.
#'
#' The key insight: pseudobulk sc-ExWAS has the same
#' statistical power as a bulk RNA-seq study with n=donors
#' -- the thousands of cells increase precision of per-donor
#' estimates but do \emph{not} increase the effective sample
#' size for donor-level inference. This function makes that
#' reality explicit for study planners.
#'
#' @references
#' Squair JW et al. (2021). Confronting false discoveries in
#' single-cell differential expression. \emph{Nat Commun}
#' 12:5692.
#'
#' @export
#' @examples
#' pw <- estimate_power(n_donors = seq(10, 50, by = 10),
#'     effect_size = 0.5, n_genes = 2000, n_celltypes = 3)
#' pw
estimate_power <- function(n_donors = seq(10, 100, by = 10),
                            effect_size = 0.5,
                            n_genes = 5000L,
                            n_celltypes = 5L,
                            alpha = 0.05,
                            dispersion = 0.1,
                            exposure_sd = 1.0) {
    stopifnot(all(n_donors >= 3))
    stopifnot(effect_size > 0)
    stopifnot(dispersion > 0)

    n_tests <- as.numeric(n_genes) * as.numeric(n_celltypes)
    ## Bonferroni-adjusted alpha (conservative bound)
    alpha_adj <- alpha / n_tests

    results <- lapply(n_donors, function(n) {
        ## Standard error of the exposure effect estimate
        ## Under the DESeq2 GLM, SE ~ sqrt(dispersion) /
        ## (exposure_sd * sqrt(n))
        se <- sqrt(dispersion) / (exposure_sd * sqrt(n))

        ## Power = P(|Z| > z_crit) where Z ~ N(effect/se, 1)
        z_crit <- qnorm(1 - alpha_adj / 2)
        ncp <- abs(effect_size) / se
        power <- pnorm(ncp - z_crit) +
                 pnorm(-ncp - z_crit)

        data.frame(
            n_donors = n,
            power = power,
            effect_size = effect_size,
            n_tests = n_tests,
            alpha_adjusted = alpha_adj,
            dispersion = dispersion,
            stringsAsFactors = FALSE)
    })

    S4Vectors::DataFrame(do.call(rbind, results))
}
