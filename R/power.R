# R/power.R
# Power analysis for sc-ExWAS study design

#' @include AllClasses.R
#' @importFrom stats rnorm rpois var power.t.test sd pnorm qnorm
NULL

#' Power analysis for sc-ExWAS study design
#'
#' Approximates the power to detect a donor-level exposure effect on one gene
#' in one cell type at the pseudobulk level, for a range of donor numbers.
#'
#' @param n_donors Integer vector. Donor sample sizes to
#'   evaluate. Default: \code{seq(10, 100, by = 10)}.
#' @param effect_size Numeric. Expected log2 fold change per
#'   unit exposure. Default 0.5.
#' @param n_genes Integer. Number of genes tested (for
#'   multiple testing correction). Default 5000.
#' @param n_celltypes Integer. Number of cell types tested.
#'   Default 5.
#' @param alpha Numeric. Family-wise significance level; the per-test level
#'   is \code{alpha / (n_genes * n_celltypes)}. Default 0.05.
#' @param dispersion Numeric. Negative-binomial dispersion \eqn{\phi} of
#'   donor pseudobulk counts, whose variance is \eqn{\mu + \phi\mu^2}.
#'   Default 0.1.
#' @param base_mean Numeric. Mean donor pseudobulk count \eqn{\mu} of the
#'   gene. Default 100; \code{Inf} gives the large-count limit.
#' @param exposure_sd Numeric. Standard deviation of the
#'   exposure variable across donors. Default 1
#'   (standardised).
#'
#' @return A \code{DataFrame} with columns: n_donors, power,
#'   effect_size, n_tests, alpha_adjusted, dispersion, base_mean.
#'
#' @details
#' For a negative-binomial GLM with a log link, the large-sample standard
#' error of the exposure coefficient on the \eqn{\log_2} scale is
#' \deqn{SE = \frac{\sqrt{1/\mu + \phi}}{\log(2)\, s_x \sqrt{n}},}
#' where \eqn{s_x} is \code{exposure_sd} and \eqn{n} the number of donors.
#' Power is that of a two-sided Wald test at the Bonferroni per-test level.
#' The approximation treats the dispersion as known and uses the normal
#' distribution, so it is optimistic for small numbers of donors; simulate
#' from the planned design when an accurate figure matters. Bonferroni is
#' conservative relative to Benjamini-Hochberg control of the false
#' discovery rate.
#'
#' Pseudobulk inference has donors, not cells, as its independent units.
#' More cells per donor raise \eqn{\mu} and so shrink the \eqn{1/\mu} term,
#' but they do not increase \eqn{n}.
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
                            base_mean = 100,
                            exposure_sd = 1.0) {
    stopifnot(all(n_donors >= 3))
    stopifnot(effect_size > 0)
    stopifnot(dispersion > 0)
    stopifnot(is.numeric(base_mean), length(base_mean) == 1L,
              base_mean > 0)
    stopifnot(exposure_sd > 0)

    n_tests <- as.numeric(n_genes) * as.numeric(n_celltypes)
    ## Bonferroni-adjusted alpha (conservative bound)
    alpha_adj <- alpha / n_tests
    z_crit <- qnorm(1 - alpha_adj / 2)

    results <- lapply(n_donors, function(n) {
        ## Large-sample SE of the log2 exposure coefficient in a
        ## negative-binomial GLM with a log link
        se <- sqrt(1 / base_mean + dispersion) /
            (log(2) * exposure_sd * sqrt(n))
        ncp <- abs(effect_size) / se
        power <- pnorm(ncp - z_crit) + pnorm(-ncp - z_crit)

        data.frame(
            n_donors = n,
            power = power,
            effect_size = effect_size,
            n_tests = n_tests,
            alpha_adjusted = alpha_adj,
            dispersion = dispersion,
            base_mean = base_mean,
            stringsAsFactors = FALSE)
    })

    S4Vectors::DataFrame(do.call(rbind, results))
}
