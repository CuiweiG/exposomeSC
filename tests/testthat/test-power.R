test_that("estimate_power returns valid DataFrame", {
    pw <- estimate_power(
        n_donors = c(10, 20, 50),
        effect_size = 0.5,
        n_genes = 1000,
        n_celltypes = 3)
    expect_s4_class(pw, "DataFrame")
    expect_equal(nrow(pw), 3L)
    expect_true(all(pw$power >= 0 & pw$power <= 1))
    ## More donors = more power
    expect_true(pw$power[3] >= pw$power[1])
})

test_that("estimate_power monotonically increases with n", {
    pw <- estimate_power(
        n_donors = seq(10, 100, by = 10),
        effect_size = 0.3)
    diffs <- diff(pw$power)
    expect_true(all(diffs >= 0))
})

test_that("estimate_power with large effect has high power", {
    pw <- estimate_power(
        n_donors = 50,
        effect_size = 2.0,
        n_genes = 100,
        n_celltypes = 1)
    expect_true(pw$power[1] > 0.8)
})

test_that("estimate_power errors on bad input", {
    expect_error(estimate_power(n_donors = 2))
    expect_error(estimate_power(effect_size = 0))
    expect_error(estimate_power(dispersion = 0))
})

test_that("estimate_power uses the log2-scale negative-binomial SE", {
    pw <- estimate_power(n_donors = 20, effect_size = 0.5, n_genes = 2000,
                         n_celltypes = 3, dispersion = 0.1, base_mean = 100)
    se <- sqrt(1 / 100 + 0.1) / (log(2) * sqrt(20))
    z <- stats::qnorm(1 - 0.05 / 6000 / 2)
    expected <- stats::pnorm(0.5 / se - z) + stats::pnorm(-0.5 / se - z)
    expect_equal(pw$power, expected)
    ## A simulated negative-binomial GLM gives about 0.56 for this design
    expect_lt(pw$power, 0.7)
    low <- estimate_power(n_donors = 20, base_mean = 20)$power
    high <- estimate_power(n_donors = 20, base_mean = Inf)$power
    expect_lt(low, high)
})
