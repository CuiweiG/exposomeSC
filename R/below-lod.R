# R/below-lod.R
# Below-LOD imputation for environmental exposure data

#' @include AllClasses.R
#' @importFrom stats rnorm runif
NULL

#' Impute exposure values below the limit of detection
#'
#' Environmental exposure measurements often contain values reported only as
#' below the limit of detection (LOD), that is, left-censored values. This
#' function replaces them by substitution, by regression on order statistics
#' or by multiple imputation from a censored log-normal model. Missing
#' values are kept apart from censored ones and are never imputed.
#'
#' @param exposure_matrix Numeric matrix. Rows = samples,
#'   columns = exposure variables.
#' @param lod Named numeric vector of positive, finite detection limits, one
#'   per exposure. Names must match column names of
#'   \code{exposure_matrix}. Exposures not in \code{lod} are left unchanged.
#' @param method Character. Imputation method:
#'   \describe{
#'     \item{\code{"lod_sqrt2"}}{Substitution with LOD / sqrt(2).}
#'     \item{\code{"lod_2"}}{Substitution with LOD / 2.}
#'     \item{\code{"ros"}}{Regression on order statistics (ROS) for one
#'       detection limit per exposure; see Details.}
#'     \item{\code{"multiple"}}{Multiple imputation from a log-normal model
#'       fitted by censored maximum likelihood to a bootstrap sample of the
#'       data, as proposed by Lubin et al. (2004); returns \code{M}
#'       completed matrices.}
#'   }
#' @param censored \code{NULL} (default) or a logical matrix with the same
#'   dimensions as \code{exposure_matrix}, \code{TRUE} where a value was
#'   reported as below the LOD. The numeric value at a censored position is
#'   ignored and may be \code{NA}. With \code{NULL}, a value is censored when
#'   it is not \code{NA} and is at or below the LOD. Either way, \code{NA}
#'   values that are not marked as censored are missing and stay \code{NA}.
#' @param M Integer. Number of imputations (only for
#'   \code{method = "multiple"}). Default 5.
#' @param verbose Logical. Report the numbers of censored and missing
#'   values per exposure. Default TRUE.
#'
#' @return For \code{method != "multiple"}: the numeric matrix with censored
#'   values replaced, of the same dimensions as the input. For
#'   \code{method = "multiple"}: an \code{S4Vectors::SimpleList} of \code{M}
#'   such matrices. In both cases a \code{data.frame} with one row per
#'   exposure in \code{lod} (\code{exposure}, \code{n_detected},
#'   \code{n_censored}, \code{n_missing}, \code{method}) records what was
#'   done; it is the \code{"lod_summary"} attribute of the matrix, or
#'   \code{metadata(result)$lod_summary} for the list.
#'
#' @details
#' \strong{Substitution} is simple, but it distorts variances and
#' associations, increasingly so as the censored fraction grows
#' (Lubin et al. 2004).
#'
#' \strong{ROS.} With a single detection limit every censored value lies
#' below every detected value, so the censored values take the lowest ranks
#' among the non-missing values of the exposure. Blom plotting positions
#' \eqn{(i - 3/8) / (n + 1/4)} are computed for all \eqn{n} non-missing
#' values, the logarithms of the detected values are regressed on the normal
#' quantiles of their plotting positions, and each censored value is
#' replaced by the value of the fitted line at one of the lowest plotting
#' positions, capped just below the LOD (Helsel 2012). ROS reproduces the
#' distribution of the censored values, not their individual values: the
#' imputed values are assigned to the censored samples in row order, so they
#' carry no information on which censored sample had the larger true value.
#' Use the ROS values for summaries of the distribution; for association
#' analyses prefer \code{method = "multiple"}. ROS is not the reverse
#' Kaplan-Meier estimator, which is not implemented here.
#'
#' \strong{Multiple imputation.} For each of the \code{M} imputations the
#' non-missing values of an exposure are resampled with replacement, the
#' mean and standard deviation of the log values are estimated by maximum
#' likelihood with the censored values contributing
#' \eqn{P(X \le \mathrm{LOD})}, and every censored value is drawn from the
#' fitted log-normal distribution truncated to \eqn{(0, \mathrm{LOD}]}. The
#' bootstrap step carries the uncertainty of the fitted distribution into
#' the imputations (Lubin et al. 2004). Analyse each completed matrix and
#' combine the estimates with Rubin's rules (Rubin 1987). The draws use the
#' random-number stream of the session, so call \code{set.seed()} first for
#' reproducible imputations.
#'
#' ROS and multiple imputation assume a log-normal distribution and need at
#' least three detected values; an exposure with fewer is imputed with
#' LOD / sqrt(2), with a warning, and this is recorded in
#' \code{lod_summary}. An exposure with more than 80\% of its non-missing
#' values censored gives a warning.
#'
#' Other implementations: \code{rexposome::ilod()} and
#' \code{rexposome::imputeLOD()} impute values below the LOD with the QRILC
#' or MinProb methods of the \pkg{imputeLCMD} package, which are not the
#' methods offered here.
#'
#' @references
#' Lubin JH, Colt JS, Camann D, et al. (2004). Epidemiologic evaluation of
#' measurement data in the presence of detection limits.
#' \emph{Environ Health Perspect} 112(17):1691--1696.
#' \doi{10.1289/ehp.7199}
#'
#' Helsel DR (2012). Statistics for Censored Environmental Data Using Minitab
#' and R. 2nd ed. Wiley. \doi{10.1002/9781118162729}
#'
#' Rubin DB (1987). Multiple Imputation for Nonresponse in Surveys. Wiley.
#' \doi{10.1002/9780470316696}
#'
#' @export
#' @examples
#' set.seed(1)
#' mat <- matrix(rlnorm(200, 0, 1), nrow = 40, ncol = 5,
#'     dimnames = list(paste0("S", 1:40),
#'         c("Pb", "Cd", "Hg", "As", "BPA")))
#' lod <- c(Pb = 0.5, Cd = 0.3, Hg = 0.2, As = 0.1, BPA = 0.4)
#' ## Values reported as below the LOD, and one genuinely missing value
#' censored <- sweep(mat, 2, lod[colnames(mat)], "<")
#' mat[censored] <- NA
#' mat["S1", "Pb"] <- NA
#' censored["S1", "Pb"] <- FALSE
#' imputed <- exposure_impute_lod(mat, lod, method = "ros",
#'     censored = censored)
#' attr(imputed, "lod_summary")
#' imputations <- exposure_impute_lod(mat, lod, method = "multiple",
#'     censored = censored, M = 3L, verbose = FALSE)
#' length(imputations)
exposure_impute_lod <- function(exposure_matrix, lod,
                                 method = c("lod_sqrt2",
                                            "lod_2",
                                            "ros",
                                            "multiple"),
                                 censored = NULL,
                                 M = 5L, verbose = TRUE) {
    method <- match.arg(method)
    if (!is.matrix(exposure_matrix) || !is.numeric(exposure_matrix)) {
        stop("exposure_matrix must be a numeric matrix.")
    }
    if (!is.numeric(lod) || is.null(names(lod)) || anyNA(names(lod)) ||
            any(!nzchar(names(lod))) || anyDuplicated(names(lod)) ||
            any(!is.finite(lod)) || any(lod <= 0)) {
        stop("lod must be a named numeric vector of positive, finite ",
             "detection limits with unique names.")
    }
    if (method == "multiple" &&
            (!is.numeric(M) || length(M) != 1L || is.na(M) || M < 1 ||
                M != round(M))) {
        stop("M must be one positive integer.")
    }
    censored <- .lod_censoring(exposure_matrix, lod, censored)

    mat <- exposure_matrix
    exp_names <- intersect(names(lod), colnames(mat))
    if (length(exp_names) == 0) {
        warning("No matching exposure names between lod and ",
                "exposure_matrix", call. = FALSE)
        return(mat)
    }

    summary <- .lod_summary(mat, censored, exp_names, method)
    for (i in seq_len(nrow(summary))) {
        row <- summary[i, ]
        if (verbose) {
            message(row$exposure, ": ", row$n_censored, " below LOD, ",
                    row$n_missing, " missing, ", row$n_detected, " detected")
        }
        n_observed <- row$n_censored + row$n_detected
        if (n_observed > 0 && row$n_censored / n_observed > 0.8) {
            warning(row$exposure, ": >80% below LOD. Consider ",
                    "excluding this variable.", call. = FALSE)
        }
        if (row$method == "lod_sqrt2" && method %in% c("ros", "multiple")) {
            warning(row$exposure, ": fewer than 3 detected values; ",
                    "imputed with LOD/sqrt(2) instead of ", method, ".",
                    call. = FALSE)
        }
    }

    if (method == "multiple") {
        imputations <- lapply(seq_len(M), function(m) {
            .impute_lod_once(mat, lod, censored, summary, "multiple")
        })
        names(imputations) <- paste0("imputation", seq_len(M))
        result <- S4Vectors::SimpleList(imputations)
        S4Vectors::metadata(result)$lod_summary <- summary
        return(result)
    }
    result <- .impute_lod_once(mat, lod, censored, summary, method)
    attr(result, "lod_summary") <- summary
    result
}

## ---- Internal: censoring indicator ----
## A logical matrix marking censored values. Values not in an exposure with
## an LOD are never censored.
.lod_censoring <- function(mat, lod, censored) {
    in_lod <- colnames(mat) %in% names(lod)
    if (is.null(censored)) {
        limits <- matrix(lod[colnames(mat)], nrow = nrow(mat),
                         ncol = ncol(mat), byrow = TRUE)
        censored <- !is.na(mat) & !is.na(limits) & mat <= limits
        censored[, !in_lod] <- FALSE
        return(censored)
    }
    if (!is.matrix(censored) || !is.logical(censored) ||
            !identical(dim(censored), dim(mat))) {
        stop("censored must be NULL or a logical matrix with the ",
             "dimensions of exposure_matrix.")
    }
    if (!is.null(dimnames(censored)) && !is.null(dimnames(mat)) &&
            !identical(unname(dimnames(censored)), unname(dimnames(mat)))) {
        stop("censored must have the row and column names of ",
             "exposure_matrix.")
    }
    if (anyNA(censored)) {
        stop("censored must not contain NA.")
    }
    if (any(censored[, !in_lod])) {
        stop("censored marks values of exposures that have no LOD in lod.")
    }
    for (name in colnames(mat)[in_lod]) {
        detected <- !censored[, name] & !is.na(mat[, name])
        if (any(mat[detected, name] <= lod[[name]])) {
            stop("Exposure '", name, "' has values at or below its LOD ",
                 "that are not marked as censored.")
        }
    }
    dimnames(censored) <- dimnames(mat)
    censored
}

## ---- Internal: per-exposure counts and the method actually used ----
.lod_summary <- function(mat, censored, exp_names, method) {
    n_censored <- colSums(censored[, exp_names, drop = FALSE])
    n_missing <- colSums(is.na(mat[, exp_names, drop = FALSE]) &
                         !censored[, exp_names, drop = FALSE])
    n_detected <- nrow(mat) - n_censored - n_missing
    used <- rep(method, length(exp_names))
    used[n_censored == 0] <- "none"
    used[n_censored > 0 & n_detected < 3 &
             method %in% c("ros", "multiple")] <- "lod_sqrt2"
    data.frame(
        exposure = exp_names,
        n_detected = as.integer(n_detected),
        n_censored = as.integer(n_censored),
        n_missing = as.integer(n_missing),
        method = used,
        row.names = NULL,
        stringsAsFactors = FALSE
    )
}

## ---- Internal: one completed matrix ----
.impute_lod_once <- function(mat, lod, censored, summary, method) {
    for (i in seq_len(nrow(summary))) {
        name <- summary$exposure[[i]]
        below <- censored[, name]
        limit <- lod[[name]]
        mat[below, name] <- switch(summary$method[[i]],
            none = mat[below, name],
            lod_sqrt2 = rep(limit / sqrt(2), sum(below)),
            lod_2 = rep(limit / 2, sum(below)),
            ros = .impute_ros(mat[!below & !is.na(mat[, name]), name],
                              limit, sum(below)),
            multiple = .impute_censored_lognormal(
                mat[!below & !is.na(mat[, name]), name], limit, sum(below)
            )
        )
    }
    mat
}

## ---- Internal: ROS imputation (Helsel 2012) ----
## Regression on order statistics with Blom plotting positions for a single
## detection limit: censored values hold ranks 1..n_below of the
## n_obs + n_below non-missing values.
.impute_ros <- function(observed, lod, n_below) {
    n_obs <- length(observed)
    n_total <- n_obs + n_below

    ## Blom plotting positions (Helsel 2012): (rank - 3/8) / (n + 1/4)
    pp_obs <- (n_below + seq_len(n_obs) - 3 / 8) / (n_total + 1 / 4)
    pp_cens <- (seq_len(n_below) - 3 / 8) / (n_total + 1 / 4)

    ## Regress log(detected) on the normal quantiles of their positions;
    ## the intercept and slope estimate the log-scale mean and SD.
    fit <- stats::lm(log(sort(observed)) ~ stats::qnorm(pp_obs))
    mu <- stats::coef(fit)[[1]]
    sigma <- stats::coef(fit)[[2]]
    imputed <- exp(mu + sigma * stats::qnorm(pp_cens))

    ## Keep imputed values positive and below the LOD
    pmax(pmin(imputed, lod * 0.999), .Machine$double.eps)
}

## ---- Internal: one multiple-imputation draw (Lubin et al. 2004) ----
## Bootstrap the non-missing values, fit a left-censored log-normal by
## maximum likelihood, then draw the censored values from the fitted
## distribution truncated to (0, lod].
.impute_censored_lognormal <- function(observed, lod, n_below) {
    n_total <- length(observed) + n_below
    is_detect <- c(rep(TRUE, length(observed)), rep(FALSE, n_below))
    values <- c(observed, rep(NA_real_, n_below))
    params <- NULL
    for (attempt in seq_len(100L)) {
        pick <- sample.int(n_total, n_total, replace = TRUE)
        boot_detected <- values[pick][is_detect[pick]]
        if (length(unique(boot_detected)) >= 2L) {
            params <- .censored_lognormal_mle(
                log(boot_detected), log(lod), sum(!is_detect[pick])
            )
            break
        }
    }
    if (is.null(params)) {
        ## no usable bootstrap sample: use the fit to the observed data
        params <- .censored_lognormal_mle(log(observed), log(lod), n_below)
    }
    upper <- stats::pnorm(log(lod), params[["mu"]], params[["sigma"]])
    exp(stats::qnorm(stats::runif(n_below, 0, upper),
                     params[["mu"]], params[["sigma"]]))
}

## Maximum likelihood for a normal sample (log values) with n_censored
## observations known only to lie below log_lod.
.censored_lognormal_mle <- function(y, log_lod, n_censored) {
    negative_loglik <- function(theta) {
        sigma <- exp(theta[[2]])
        -(sum(stats::dnorm(y, theta[[1]], sigma, log = TRUE)) +
              n_censored * stats::pnorm(log_lod, theta[[1]], sigma,
                                        log.p = TRUE))
    }
    start <- c(mean(y), log(max(stats::sd(y), 1e-3)))
    fit <- stats::optim(start, negative_loglik, method = "BFGS")
    if (fit$convergence != 0L) {
        warning("Censored log-normal fit did not converge; ",
                "imputations may be unreliable.", call. = FALSE)
    }
    c(mu = fit$par[[1]], sigma = exp(fit$par[[2]]))
}
