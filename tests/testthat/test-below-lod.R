# tests/testthat/test-below-lod.R

.lod_matrix <- function() {
    set.seed(42)
    mat <- matrix(rlnorm(200), nrow = 40, ncol = 5,
        dimnames = list(paste0("S", 1:40),
            c("Pb", "Cd", "Hg", "As", "BPA")))
    lod <- c(Pb = 0.5, Cd = 0.3)
    censored <- matrix(FALSE, nrow(mat), ncol(mat), dimnames = dimnames(mat))
    censored[, "Pb"] <- mat[, "Pb"] < 0.5
    censored[, "Cd"] <- mat[, "Cd"] < 0.3
    mat[censored] <- NA
    list(mat = mat, lod = lod, censored = censored)
}

test_that("substitution replaces censored values only", {
    d <- .lod_matrix()
    d$mat["S1", "Pb"] <- NA
    d$censored["S1", "Pb"] <- FALSE

    result <- exposure_impute_lod(d$mat, d$lod, method = "lod_sqrt2",
        censored = d$censored, verbose = FALSE)

    expect_equal(dim(result), dim(d$mat))
    expect_true(is.na(result["S1", "Pb"]))
    expect_equal(unique(result[d$censored[, "Pb"], "Pb"]), 0.5 / sqrt(2))
    expect_equal(unique(result[d$censored[, "Cd"], "Cd"]), 0.3 / sqrt(2))
    detected <- !d$censored & !is.na(d$mat)
    expect_identical(result[detected], d$mat[detected])

    half <- exposure_impute_lod(d$mat, d$lod, method = "lod_2",
        censored = d$censored, verbose = FALSE)
    expect_equal(unique(half[d$censored[, "Pb"], "Pb"]), 0.25)

    summary <- attr(result, "lod_summary")
    expect_identical(summary$exposure, c("Pb", "Cd"))
    expect_identical(summary$n_censored,
        as.integer(colSums(d$censored[, c("Pb", "Cd")])))
    expect_identical(summary$n_missing, c(1L, 0L))
})

test_that("without an indicator, NA is missing and values at the LOD are censored", {
    mat <- matrix(c(0.2, 0.5, NA, 0.9, 1.4, 2.0), ncol = 1,
        dimnames = list(paste0("S", 1:6), "E1"))
    result <- exposure_impute_lod(mat, c(E1 = 0.5), method = "lod_2",
        verbose = FALSE)
    expect_equal(result[, "E1"], c(0.25, 0.25, NA, 0.9, 1.4, 2.0),
        ignore_attr = TRUE)
    expect_identical(attr(result, "lod_summary")$n_missing, 1L)
    expect_identical(attr(result, "lod_summary")$n_censored, 2L)
})

test_that("ROS follows the Blom-position regression", {
    set.seed(123)
    n <- 100
    x <- rlnorm(n, 1, 0.5)
    lod_val <- stats::qlnorm(0.3, 1, 0.5)
    mat <- matrix(x, ncol = 1, dimnames = list(paste0("S", seq_len(n)), "E1"))
    below <- mat[, 1] <= lod_val
    mat[below, 1] <- NA

    result <- exposure_impute_lod(mat, c(E1 = lod_val), method = "ros",
        censored = matrix(below, ncol = 1, dimnames = dimnames(mat)),
        verbose = FALSE)
    imputed <- result[below, "E1"]

    n_below <- sum(below)
    observed <- sort(x[!below])
    pp <- (seq_len(n) - 3 / 8) / (n + 1 / 4)
    fit <- lm(log(observed) ~ qnorm(pp[-seq_len(n_below)]))
    expected <- exp(coef(fit)[[1]] + coef(fit)[[2]] *
        qnorm(pp[seq_len(n_below)]))
    expect_equal(unname(imputed), pmin(expected, lod_val * 0.999))
    expect_true(all(imputed < lod_val) && all(imputed > 0))
    expect_true(all(diff(imputed) >= 0))
    expect_identical(attr(result, "lod_summary")$method, "ros")
})

test_that("ROS plotting positions ignore missing values", {
    d <- .lod_matrix()
    with_missing <- d$mat
    with_missing[c("S2", "S3"), "Pb"] <- NA
    censored <- d$censored
    censored[c("S2", "S3"), "Pb"] <- FALSE
    censored[, "Cd"] <- FALSE
    keep <- setdiff(rownames(d$mat), c("S2", "S3"))

    full <- exposure_impute_lod(with_missing, d$lod["Pb"], method = "ros",
        censored = censored, verbose = FALSE)
    dropped <- exposure_impute_lod(with_missing[keep, ], d$lod["Pb"],
        method = "ros", censored = censored[keep, ], verbose = FALSE)
    expect_equal(full[keep, "Pb"], dropped[, "Pb"])
    expect_true(all(is.na(full[c("S2", "S3"), "Pb"])))
})

test_that("too few detected values fall back to LOD/sqrt(2) with a warning", {
    mat <- matrix(c(NA, NA, NA, NA, 1.2, 1.5), ncol = 1,
        dimnames = list(paste0("S", 1:6), "E1"))
    censored <- matrix(c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE), ncol = 1,
        dimnames = dimnames(mat))
    expect_warning(
        result <- exposure_impute_lod(mat, c(E1 = 1), method = "ros",
            censored = censored, verbose = FALSE),
        "fewer than 3 detected values"
    )
    expect_equal(unique(result[1:4, "E1"]), 1 / sqrt(2))
    expect_identical(attr(result, "lod_summary")$method, "lod_sqrt2")
})

test_that("multiple imputation returns M completed matrices", {
    d <- .lod_matrix()
    d$mat["S1", "Pb"] <- NA
    d$censored["S1", "Pb"] <- FALSE

    set.seed(1)
    result <- exposure_impute_lod(d$mat, d$lod, method = "multiple",
        censored = d$censored, M = 3L, verbose = FALSE)
    set.seed(1)
    again <- exposure_impute_lod(d$mat, d$lod, method = "multiple",
        censored = d$censored, M = 3L, verbose = FALSE)

    expect_s4_class(result, "SimpleList")
    expect_length(result, 3L)
    expect_identical(as.list(result), as.list(again))
    expect_false(identical(result[[1]], result[[2]]))
    for (m in seq_along(result)) {
        completed <- result[[m]]
        expect_equal(dim(completed), dim(d$mat))
        expect_true(is.na(completed["S1", "Pb"]))
        pb <- completed[d$censored[, "Pb"], "Pb"]
        expect_true(all(pb > 0 & pb <= 0.5))
        detected <- !d$censored & !is.na(d$mat)
        expect_identical(completed[detected], d$mat[detected])
    }
    expect_identical(S4Vectors::metadata(result)$lod_summary$method,
        c("multiple", "multiple"))
})

test_that("the censored log-normal fit recovers known parameters", {
    set.seed(2026)
    y <- rnorm(5000, mean = 1, sd = 0.5)
    log_lod <- qnorm(0.3, 1, 0.5)
    fit <- exposomeSC:::.censored_lognormal_mle(
        y[y > log_lod], log_lod, sum(y <= log_lod))
    expect_equal(unname(fit[["mu"]]), 1, tolerance = 0.03)
    expect_equal(unname(fit[["sigma"]]), 0.5, tolerance = 0.03)
})

test_that("multiple imputation reproduces the censored distribution", {
    set.seed(7)
    n <- 4000
    x <- rlnorm(n, 0, 1)
    lod_val <- qlnorm(0.4, 0, 1)
    mat <- matrix(x, ncol = 1, dimnames = list(NULL, "E1"))
    below <- matrix(x <= lod_val, ncol = 1, dimnames = list(NULL, "E1"))
    mat[below] <- NA
    result <- exposure_impute_lod(mat, c(E1 = lod_val), method = "multiple",
        censored = below, M = 2L, verbose = FALSE)
    imputed <- result[[1]][below[, 1], "E1"]
    ## mean of a standard log-normal truncated at the LOD
    z <- log(lod_val)
    expected_mean <- exp(0.5) * pnorm(z - 1) / pnorm(z)
    expect_equal(mean(imputed), expected_mean, tolerance = 0.05)
})

test_that("invalid inputs are rejected", {
    d <- .lod_matrix()
    expect_error(
        exposure_impute_lod(d$mat, d$lod, method = "kaplan_meier"),
        "should be one of"
    )
    expect_error(
        exposure_impute_lod(d$mat, c(Pb = 0), verbose = FALSE),
        "positive, finite"
    )
    expect_error(
        exposure_impute_lod(d$mat, d$lod, censored = d$censored[, 1:2],
            verbose = FALSE),
        "dimensions"
    )
    bad <- d$censored
    bad[1, 1] <- NA
    expect_error(
        exposure_impute_lod(d$mat, d$lod, censored = bad, verbose = FALSE),
        "must not contain NA"
    )
    bad <- d$censored
    bad[1, "Hg"] <- TRUE
    expect_error(
        exposure_impute_lod(d$mat, d$lod, censored = bad, verbose = FALSE),
        "no LOD"
    )
    low <- d$mat
    low[which(!d$censored[, "Pb"])[1], "Pb"] <- 0.1
    expect_error(
        exposure_impute_lod(low, d$lod, censored = d$censored,
            verbose = FALSE),
        "not marked as censored"
    )
})

test_that("the censored fit copes with detected values that are all equal", {
    fit <- exposomeSC:::.censored_lognormal_mle(log(c(1, 1, 1)), log(0.9),
        7L)
    expect_true(all(is.finite(fit)))
    expect_lt(abs(fit[["mu"]]), 5)
    expect_gt(fit[["sigma"]], 0.01)
    expect_lt(fit[["sigma"]], 10)
    set.seed(3)
    mat <- cbind(A = c(1, 1, 1, 0.5, 0.4, 0.3, 0.6, 0.2, 0.7, 0.8))
    rownames(mat) <- paste0("D", seq_len(10))
    result <- exposure_impute_lod(mat, c(A = 0.9), method = "multiple",
        M = 3L, verbose = FALSE)
    imputed <- unlist(lapply(result, function(x) x[4:10, "A"]))
    expect_true(all(imputed > 0 & imputed <= 0.9))
    expect_gt(min(imputed), 1e-3)
})

test_that("a censored fit that fails is skipped or reported", {
    ## a detected value of zero has no finite logarithm, so no fit exists
    expect_null(exposomeSC:::.censored_lognormal_mle(c(-Inf, 0, 1), log(0.9),
        2L))
    ## resamples that hold the zero fail and are skipped, the others fit
    set.seed(1)
    drawn <- exposomeSC:::.impute_censored_lognormal(c(0, 1, 2), 0.9, 4L)
    expect_length(drawn, 4L)
    expect_true(all(drawn > 0 & drawn <= 0.9))
    ## when every fit fails the exposure is named in the error
    expect_null(exposomeSC:::.impute_censored_lognormal(c(0, 0, 0), 0.9, 7L))
    mat <- cbind(A = c(0, 0, 0, 0.5, 0.4))
    censored <- cbind(A = c(FALSE, FALSE, FALSE, TRUE, TRUE))
    summary <- data.frame(exposure = "A", method = "multiple")
    expect_error(
        exposomeSC:::.impute_lod_once(mat, c(A = 0.9), censored, summary,
            "multiple"),
        "Exposure 'A'")
})

test_that("infinite values are refused and draws stay inside (0, LOD]", {
    mat <- cbind(A = c(Inf, 1.5, 2, 0.5, 0.4))
    expect_error(
        exposure_impute_lod(mat, c(A = 0.9), method = "multiple", M = 2L,
            verbose = FALSE),
        "infinite")
    ## a truncation point 100 standard deviations below the mean
    set.seed(1)
    drawn <- exposomeSC:::.draw_truncated_lognormal(5L, 1, 10, 0.1)
    expect_true(all(drawn > 0 & drawn <= 1))
})

test_that("an infinite value at a censored position is ignored", {
    mat <- cbind(A = c(-Inf, 1.5, 2, 3, NA, 1.2))
    censored <- cbind(A = c(TRUE, FALSE, FALSE, FALSE, TRUE, FALSE))
    imputed <- exposure_impute_lod(mat, c(A = 0.9), method = "lod_2",
        censored = censored, verbose = FALSE)
    expect_equal(imputed[, "A"], c(0.45, 1.5, 2, 3, 0.45, 1.2))
    ## without an indicator, -Inf is at or below the LOD and so censored
    expect_equal(
        exposure_impute_lod(cbind(A = c(-Inf, 1.5, 2)), c(A = 0.9),
            method = "lod_2", verbose = FALSE)[, "A"],
        c(0.45, 1.5, 2))
    ## a detected infinite value is an error, whatever the method
    expect_error(
        exposure_impute_lod(cbind(A = c(Inf, 1.5, 2)), c(A = 0.9),
            method = "lod_sqrt2", verbose = FALSE),
        "infinite values that are not censored")
})

test_that("a censored fit that reached the maximum is kept", {
    ## both starts reach the maximum here, and the fit from the first start
    ## is kept: the second start replaces it only when it is better by more
    ## than 1e-6 relative
    set.seed(3)
    x <- stats::rlnorm(1000, 0, 2)
    lod <- stats::quantile(x, 0.5, names = FALSE)
    y <- log(x[x > lod])
    k <- sum(x <= lod)
    nll <- function(theta) {
        -(sum(stats::dnorm(y, theta[1], exp(theta[2]), log = TRUE)) +
            k * stats::pnorm(log(lod), theta[1], exp(theta[2]), log.p = TRUE))
    }
    spread <- max(stats::sd(y), abs(mean(y) - log(lod)) / 2, 0.05)
    first <- stats::optim(c(mean(y), log(spread)), nll, method = "BFGS")
    fit <- exposomeSC:::.censored_lognormal_mle(y, log(lod), k)
    expect_equal(unname(fit), c(first$par[1], exp(first$par[2])))
})

test_that("the censored fit is not misled by many censored values", {
    ## from the mean of the detected values the optimiser reported
    ## convergence far below the data, at means of about -380 and -850;
    ## the likelihood is highest at these values
    fit <- exposomeSC:::.censored_lognormal_mle(rep(0.5, 6), 0, 1000L)
    expect_equal(unname(fit), c(-3.591, 1.430), tolerance = 1e-3)
    fit <- exposomeSC:::.censored_lognormal_mle(
        log(c(1.284, 1.275, 1.312, 1.170)), 0, 1000L)
    expect_equal(unname(fit), c(-1.830, 0.690), tolerance = 1e-3)
    mat <- cbind(A = c(rep(exp(0.5), 6), rep(NA, 1000)))
    censored <- cbind(A = rep(c(FALSE, TRUE), c(6, 1000)))
    set.seed(1)
    expect_warning(
        result <- exposure_impute_lod(mat, c(A = 1), method = "multiple",
            censored = censored, M = 2L, verbose = FALSE),
        ">80% below LOD")
    imputed <- unlist(lapply(result, function(x) x[-(1:6), "A"]))
    expect_true(all(imputed > 0 & imputed <= 1))
    expect_gt(stats::median(imputed), 1e-3)
    ## four equal values far above the LOD: BFGS stops at its iteration
    ## limit from both starts, and only the Nelder-Mead retry from the
    ## spread rule reaches the maximum
    fit <- exposomeSC:::.censored_lognormal_mle(rep(10, 4), 0, 100L)
    expect_equal(unname(fit), c(-38.886, 22.110), tolerance = 1e-3)
})

test_that("lod names must be columns of the exposure matrix", {
    mat <- cbind(PM2.5 = c(0.1, 0.5, 1, 2, 3), NO2 = c(0.2, 0.4, 3, 5, 7))
    rownames(mat) <- paste0("D", 1:5)
    expect_error(
        exposure_impute_lod(mat, lod = c(PM25 = 0.3), method = "multiple",
            verbose = FALSE),
        "No name of lod is a column of exposure_matrix")
    expect_warning(
        imputed <- exposure_impute_lod(mat, lod = c(PM2.5 = 0.3, PM25 = 0.3),
            method = "lod_2", verbose = FALSE),
        "ignored: PM25")
    expect_identical(attr(imputed, "lod_summary")$exposure, "PM2.5")
    expect_equal(imputed[1, "PM2.5"], 0.15)
})
