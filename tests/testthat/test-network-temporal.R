# tests/testthat/test-network-temporal.R

test_that("TemporalNetwork S4 class works", {
    tn <- new("TemporalNetwork")
    expect_s4_class(tn, "TemporalNetwork")
    expect_output(show(tn), "TemporalNetwork")
})

test_that("run_temporal_network rejects invalid input", {
    # Mismatched lengths
    expect_error(
        run_temporal_network(
            list(a = "x", b = "y"),
            list(a = "x"),
            celltype = "Mono"),
        "same length")

    # Too few time points
    expect_error(
        run_temporal_network(
            list(a = "x"),
            list(a = "x"),
            celltype = "Mono"))
})
