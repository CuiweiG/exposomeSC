# tests/testthat/test-below-lod.R

test_that("exposure_impute_lod works with lod_sqrt2", {
    for (f in list.files("../../R", full.names = TRUE))
        source(f, local = TRUE)

    set.seed(42)
    mat <- matrix(rlnorm(200), nrow = 40, ncol = 5,
        dimnames = list(paste0("S", 1:40),
            c("Pb", "Cd", "Hg", "As", "BPA")))
    lod <- c(Pb = 0.5, Cd = 0.3)
    mat[mat[, "Pb"] < 0.5, "Pb"] <- NA
    mat[mat[, "Cd"] < 0.3, "Cd"] <- NA

    result <- exposure_impute_lod(mat, lod,
        method = "lod_sqrt2", verbose = FALSE)

    expect_equal(dim(result), dim(mat))
    expect_false(any(is.na(result[, c("Pb", "Cd")])))
    expect_true(all(result[, "Pb"] > 0))
    ## Imputed values should be LOD/sqrt(2) = 0.3536
    imputed_pb <- result[is.na(mat[, "Pb"]), "Pb"]
    expect_equal(unique(imputed_pb), 0.5 / sqrt(2),
        tolerance = 1e-6)
})

test_that("multiple imputation returns M datasets", {
    for (f in list.files("../../R", full.names = TRUE))
        source(f, local = TRUE)

    set.seed(42)
    mat <- matrix(rlnorm(200), nrow = 40, ncol = 5,
        dimnames = list(paste0("S", 1:40),
            c("Pb", "Cd", "Hg", "As", "BPA")))
    lod <- c(Pb = 0.5, Cd = 0.3)
    mat[mat[, "Pb"] < 0.5, "Pb"] <- NA
    mat[mat[, "Cd"] < 0.3, "Cd"] <- NA

    result <- exposure_impute_lod(mat, lod,
        method = "multiple", M = 3, verbose = FALSE)

    expect_type(result, "list")
    expect_length(result, 3)
    expect_equal(dim(result[[1]]), dim(mat))
    ## Each imputation should differ
    expect_false(identical(result[[1]], result[[2]]))
})

test_that("ROS imputation gives values below LOD", {
    for (f in list.files("../../R", full.names = TRUE))
        source(f, local = TRUE)

    set.seed(42)
    mat <- matrix(rlnorm(200), nrow = 40, ncol = 5,
        dimnames = list(paste0("S", 1:40),
            c("Pb", "Cd", "Hg", "As", "BPA")))
    lod <- c(Pb = 0.5)
    mat[mat[, "Pb"] < 0.5, "Pb"] <- NA

    result <- exposure_impute_lod(mat, lod,
        method = "kaplan_meier", verbose = FALSE)

    imputed <- result[is.na(mat[, "Pb"]), "Pb"]
    expect_true(all(imputed < 0.5))
    expect_true(all(imputed > 0))
})

test_that("ROS uses Blom positions and lognormal fit", {
    for (f in list.files("../../R", full.names = TRUE))
        source(f, local = TRUE)

    set.seed(123)
    ## Generate lognormal data with known parameters
    n <- 100
    true_mu <- 1.0; true_sigma <- 0.5
    x <- rlnorm(n, true_mu, true_sigma)
    lod_val <- stats::qlnorm(0.3, true_mu, true_sigma)
    ## About 30% censored

    mat <- matrix(x, ncol = 1,
        dimnames = list(paste0("S", seq_len(n)), "E1"))
    mat[mat[, 1] <= lod_val, 1] <- NA
    lod <- c(E1 = lod_val)

    result <- exposure_impute_lod(mat, lod,
        method = "kaplan_meier", verbose = FALSE)

    imputed <- result[is.na(mat[, "E1"]), "E1"]

    ## All imputed values below LOD
    expect_true(all(imputed < lod_val))
    expect_true(all(imputed > 0))

    ## Imputed values should be monotonically increasing
    ## (ordered by Blom plotting positions)
    expect_true(all(diff(imputed) >= 0))

    ## Imputed values should roughly match lognormal
    ## distribution: mean of imputed should be reasonable
    ## (not too far from LOD/sqrt(2) = ~0.71 * LOD)
    expect_true(mean(imputed) < lod_val)
    expect_true(mean(imputed) > lod_val * 0.1)
})
