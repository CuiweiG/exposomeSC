# R/below-lod.R
# Below-LOD imputation for environmental exposure data

#' @include AllClasses.R
#' @importFrom stats rnorm runif
NULL

#' Impute exposure values below the limit of detection
#'
#' Environmental exposure measurements frequently contain
#' values below the limit of detection (LOD). Naive
#' approaches (substitution with LOD/2, LOD/sqrt(2), or 0)
#' introduce bias. This function provides multiple
#' principled imputation methods.
#'
#' @param exposure_matrix Numeric matrix. Rows = samples,
#'   columns = exposure variables.
#' @param lod Named numeric vector. LOD for each exposure.
#'   Names must match column names of \code{exposure_matrix}.
#'   Exposures not in \code{lod} are left unchanged.
#' @param method Character. Imputation method:
#'   \describe{
#'     \item{\code{"lod_sqrt2"}}{LOD / sqrt(2). Simple,
#'       commonly used. Unbiased under lognormal assumption
#'       when censoring < 30\%.}
#'     \item{\code{"lod_2"}}{LOD / 2. Common but
#'       slightly biased.}
#'     \item{\code{"kaplan_meier"}}{Reverse Kaplan-Meier
#'       (ROS method). Uses observed distribution to impute.
#'       Best for moderate censoring (10-50\%).}
#'     \item{\code{"multiple"}}{Multiple imputation via
#'       truncated lognormal. Returns \code{M} imputed
#'       datasets for subsequent pooling. Gold standard
#'       for heavy censoring.}
#'   }
#' @param M Integer. Number of imputations (only for
#'   \code{method = "multiple"}). Default 5.
#' @param verbose Logical. Print summary of below-LOD
#'   proportions. Default TRUE.
#'
#' @return For \code{method != "multiple"}: imputed numeric
#'   matrix (same dimensions as input).
#'   For \code{method = "multiple"}: list of \code{M} imputed
#'   matrices.
#'
#' @details
#' Below-LOD values are identified as observations <= LOD for
#' each exposure. The percentage of below-LOD values per
#' exposure is reported if \code{verbose = TRUE}.
#'
#' \strong{Choosing a method:}
#' \tabular{lll}{
#'   Censoring \tab Recommended \tab Rationale \cr
#'   < 10\% \tab lod_sqrt2 \tab Minimal bias \cr
#'   10-30\% \tab kaplan_meier \tab Distribution-based \cr
#'   30-50\% \tab multiple \tab Proper uncertainty \cr
#'   > 50\% \tab Exclude variable \tab Too much missing \cr
#' }
#'
#' The \code{kaplan_meier} method implements Regression on
#' Order Statistics (ROS): observed values are probability-
#' plotted, a lognormal line is fit, and below-LOD values
#' are imputed from the fitted distribution.
#'
#' The \code{multiple} method draws from a truncated
#' lognormal distribution (truncated at LOD) fitted to
#' observed data, generating \code{M} complete datasets.
#' Analysis should be run on each dataset and results
#' pooled using Rubin's rules.
#'
#' @references
#' Lubin JH et al. (2004). Epidemiologic evaluation of measurement
#' data in the presence of detection limits.
#' \emph{Environ Health Perspect} 112:1691-1696.
#'
#' Helsel DR (2012). Statistics for Censored Environmental
#' Data Using Minitab and R. 2nd ed. Wiley.
#'
#' @export
#' @examples
#' set.seed(1)
#' mat <- matrix(rlnorm(200, 0, 1), nrow = 40, ncol = 5,
#'     dimnames = list(paste0("S", 1:40),
#'         c("Pb", "Cd", "Hg", "As", "BPA")))
#' lod <- c(Pb = 0.5, Cd = 0.3, Hg = 0.2, As = 0.1, BPA = 0.4)
#' ## Set some values below LOD
#' for (j in seq_len(5)) {
#'     below <- mat[, j] < lod[j]
#'     mat[below, j] <- NA
#' }
#' imputed <- exposure_impute_lod(mat, lod, method = "lod_sqrt2")
#' dim(imputed)
exposure_impute_lod <- function(exposure_matrix, lod,
                                 method = c("lod_sqrt2",
                                            "lod_2",
                                            "kaplan_meier",
                                            "multiple"),
                                 M = 5L, verbose = TRUE) {
    method <- match.arg(method)
    stopifnot(is.matrix(exposure_matrix))
    stopifnot(is.numeric(lod))

    mat <- exposure_matrix
    exp_names <- intersect(names(lod), colnames(mat))

    if (length(exp_names) == 0) {
        warning("No matching exposure names between lod and ",
                "exposure_matrix", call. = FALSE)
        return(mat)
    }

    ## Detect below-LOD (NA or <= LOD)
    for (exp_name in exp_names) {
        col <- mat[, exp_name]
        below <- is.na(col) | col <= lod[exp_name]
        pct <- round(100 * sum(below) / length(col), 1)

        if (verbose) {
            message(exp_name, ": ", sum(below), "/",
                length(col), " below LOD (",
                pct, "%)")
        }

        if (pct > 80) {
            warning(exp_name, ": >80% below LOD. Consider ",
                    "excluding this variable.", call. = FALSE)
        }
    }

    if (method == "multiple") {
        return(.impute_multiple(mat, lod, exp_names, M))
    }

    ## Single imputation methods
    for (exp_name in exp_names) {
        col <- mat[, exp_name]
        below <- is.na(col) | col <= lod[exp_name]
        if (!any(below)) next

        imputed_vals <- switch(method,
            lod_sqrt2 = rep(lod[exp_name] / sqrt(2),
                sum(below)),
            lod_2 = rep(lod[exp_name] / 2, sum(below)),
            kaplan_meier = .impute_ros(col, lod[exp_name],
                sum(below))
        )

        mat[below, exp_name] <- imputed_vals
    }

    mat
}

## ---- Internal: ROS imputation (Helsel 2012) ----
## Implements Regression on Order Statistics with Blom
## plotting positions, as described in Helsel (2012)
## "Statistics for Censored Environmental Data", Ch. 6.
.impute_ros <- function(x, lod, n_below) {
    ## Regression on Order Statistics
    observed <- x[!is.na(x) & x > lod]
    n_obs <- length(observed)
    n_total <- n_obs + n_below

    if (n_obs < 3) {
        ## Fall back to LOD/sqrt(2) if too few observed
        return(rep(lod / sqrt(2), n_below))
    }

    ## Step 1: Assign ranks to ALL observations
    ## Censored values get ranks 1..n_below
    ## Observed values get ranks (n_below+1)..n_total
    ## (sorted by value within observed)
    obs_sorted <- sort(observed)

    ## Step 2: Blom plotting positions (Helsel 2012, p.56)
    ## pp_i = (rank_i - 3/8) / (n + 1/4)
    ## This is preferred over (i-0.5)/n for normal scores
    ranks_censored <- seq_len(n_below)
    ranks_observed <- n_below + seq_len(n_obs)
    pp_obs <- (ranks_observed - 3 / 8) / (n_total + 1 / 4)

    ## Step 3: Normal scores regression on observed data
    ## Regress log(observed) on qnorm(plotting_position)
    z_obs <- stats::qnorm(pp_obs)
    log_obs <- log(obs_sorted)

    fit <- tryCatch(
        stats::lm(log_obs ~ z_obs),
        error = function(e) NULL)

    if (is.null(fit)) {
        ## Fallback to MLE
        mu <- mean(log(observed))
        sigma <- stats::sd(log(observed))
        if (sigma <= 0) sigma <- 0.1
    } else {
        ## Intercept = mu, slope = sigma of lognormal
        mu <- stats::coef(fit)[1]
        sigma <- abs(stats::coef(fit)[2])
        if (sigma <= 0) sigma <- 0.1
    }

    ## Step 4: Impute censored values using the fitted line
    ## Censored values get Blom positions for ranks 1..n_below
    pp_cens <- (ranks_censored - 3 / 8) / (n_total + 1 / 4)
    z_cens <- stats::qnorm(pp_cens)
    imputed <- exp(mu + sigma * z_cens)

    ## Ensure all imputed values are below LOD
    imputed <- pmin(imputed, lod * 0.999)
    ## Ensure positive
    imputed <- pmax(imputed, .Machine$double.eps)
    imputed
}

## ---- Internal: Multiple imputation ----
.impute_multiple <- function(mat, lod, exp_names, M) {
    result <- vector("list", M)

    for (m in seq_len(M)) {
        mat_m <- mat
        for (exp_name in exp_names) {
            col <- mat[, exp_name]
            below <- is.na(col) | col <= lod[exp_name]
            if (!any(below)) next

            observed <- col[!is.na(col) & col > lod[exp_name]]
            n_obs <- length(observed)

            if (n_obs < 3) {
                ## Fallback: uniform(0, LOD)
                mat_m[below, exp_name] <- stats::runif(
                    sum(below), 0, lod[exp_name])
                next
            }

            log_obs <- log(observed)
            mu <- mean(log_obs)
            sigma <- stats::sd(log_obs)
            if (sigma <= 0) sigma <- 0.1

            ## Draw from truncated lognormal
            ## Rejection sampling
            draws <- numeric(0)
            max_iter <- sum(below) * 100
            iter <- 0
            while (length(draws) < sum(below) &&
                   iter < max_iter) {
                candidates <- exp(stats::rnorm(
                    sum(below) * 2, mu, sigma))
                candidates <- candidates[
                    candidates <= lod[exp_name] &
                    candidates > 0]
                draws <- c(draws, candidates)
                iter <- iter + 1
            }

            if (length(draws) >= sum(below)) {
                mat_m[below, exp_name] <- draws[
                    seq_len(sum(below))]
            } else {
                ## Fallback
                mat_m[below, exp_name] <- lod[exp_name] /
                    sqrt(2)
            }
        }
        result[[m]] <- mat_m
    }

    result
}
