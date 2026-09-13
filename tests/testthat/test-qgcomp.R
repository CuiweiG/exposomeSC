# tests/testthat/test-qgcomp.R

test_that("run_mixture_qgcomp works", {
    skip_if_not_installed("qgcomp")

    set.seed(42)
    counts <- matrix(rpois(10000, 10), nrow = 100,
        dimnames = list(paste0("G", 1:100), paste0("c", 1:100)))
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = counts),
        colData = S4Vectors::DataFrame(
            cell_id = paste0("c", 1:100),
            donor_id = rep(paste0("D", 1:10), each = 10),
            cell_type = rep(c("Mono", "NK"), 50)))
    exp_mat <- matrix(rnorm(30), nrow = 10,
        dimnames = list(paste0("D", 1:10),
            c("E1", "E2", "E3")))

    scee <- build_scee(sce, exp_mat, sample_col = "donor_id")

    mix <- run_mixture_qgcomp(scee,
        exposures = c("E1", "E2", "E3"),
        celltype = "Mono", celltype_col = "cell_type",
        sample_col = "donor_id", q = 2L, min_cells = 3L)

    expect_type(mix, "list")
    expect_true("positive_weights" %in% names(mix))
    expect_true("negative_weights" %in% names(mix))
    expect_true("mixture_pvalue" %in% names(mix))
    expect_true("mixture_ci" %in% names(mix))
    expect_equal(length(mix$mixture_ci), 2)

    ## The mixture p-value is that of psi, not of the intercept
    psi_index <- match("psi1", names(mix$fit$coef))
    expect_equal(mix$mixture_pvalue, mix$fit$pval[psi_index])
    expect_false(isTRUE(all.equal(mix$mixture_pvalue, mix$fit$pval[1])))

    ## Weights should exist
    all_weights <- c(mix$positive_weights,
        mix$negative_weights)
    expect_true(length(all_weights) > 0)
})
