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
