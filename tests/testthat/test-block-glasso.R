# tests/testthat/test-block-glasso.R

chain_data <- function(n = 2000L, p = 6L, seed = 1L) {
    omega <- diag(p)
    for (i in seq_len(p - 1L)) {
        omega[i, i + 1L] <- omega[i + 1L, i] <- 0.4
    }
    set.seed(seed)
    list(
        omega = omega,
        x = MASS::mvrnorm(n, rep(0, p), solve(omega)),
        blocks = rep(c("transcript", "metabolite"), each = p / 2L)
    )
}

test_that(".block_glasso is the graphical lasso on the correlation scale", {
    skip_if_not_installed("glasso")
    d <- chain_data()
    S <- stats::cov(d$x)
    fit <- exposomeSC:::.block_glasso(S, d$blocks,
        lambda_w = 0.2, lambda_b = 0.2)

    rho <- matrix(0.2, 6L, 6L)
    diag(rho) <- 0
    direct <- glasso::glasso(stats::cov2cor(S), rho = rho,
        penalize.diagonal = FALSE)$wi
    expect_equal(fit$wi, (direct + t(direct)) / 2)
    expect_true(isSymmetric(fit$wi))
    expect_identical(fit$rho, rho)
})

test_that(".block_glasso recovers a known sparse precision matrix", {
    skip_if_not_installed("glasso")
    d <- chain_data()
    fit <- exposomeSC:::.block_glasso(stats::cov(d$x), d$blocks,
        lambda_w = 0.2, lambda_b = 0.2)
    truth <- (d$omega != 0) * 1L
    diag(truth) <- 0L
    expect_identical(unname(fit$adj), unname(truth))

    off <- row(d$omega) != col(d$omega)
    edge <- off & truth == 1L
    expect_gt(min(abs(fit$wi[edge])), 0.2)
    expect_lt(max(abs(fit$wi[off & truth == 0L])), 1e-10)
})

test_that("the between-block penalty acts only on cross-block pairs", {
    skip_if_not_installed("glasso")
    d <- chain_data()
    fit <- exposomeSC:::.block_glasso(stats::cov(d$x), d$blocks,
        lambda_w = 0.2, lambda_b = 10)
    cross <- outer(d$blocks, d$blocks, "!=")
    expect_true(all(fit$adj[cross] == 0L))
    expect_true(all(fit$rho[cross] == 10))
    expect_true(all(fit$rho[!cross & row(cross) != col(cross)] == 0.2))
    ## the within-block chain edges survive
    expect_identical(fit$adj[1, 2], 1L)
    expect_identical(fit$adj[5, 6], 1L)
})

test_that(".block_glasso refuses bad penalties and constant features", {
    skip_if_not_installed("glasso")
    expect_error(
        exposomeSC:::.block_glasso(diag(3), 1:3, lambda_w = -1),
        "non-negative"
    )
    expect_error(
        exposomeSC:::.block_glasso(diag(3), 1:3, lambda_b = c(0.1, 0.2)),
        "single"
    )
    expect_error(exposomeSC:::.block_glasso(diag(3), 1:2), "one entry")
    S <- diag(3)
    S[2, 2] <- 0
    expect_error(exposomeSC:::.block_glasso(S, 1:3), "positive, finite")
})
