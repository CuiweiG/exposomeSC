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
        q = 2L, min_cells = 3L)

    expect_s4_class(mix, "DataFrame")
    expect_identical(nrow(mix), 1L)
    expect_identical(
        colnames(mix),
        c("celltype", "mixture_effect", "mixture_se", "ci_lower",
          "ci_upper", "mixture_pvalue", "n_donors", "n_genes", "method")
    )
    expect_identical(mix$method, "qgcomp.noboot")
    expect_lt(mix$ci_lower, mix$ci_upper)
    meta <- S4Vectors::metadata(mix)
    expect_true(all(c("positive_weights", "negative_weights", "fit",
                      "target_genes") %in% names(meta)))
    expect_equal(mix$mixture_effect, unname(meta$fit$psi))
    ## a t interval on the residual degrees of freedom, wider than a
    ## normal interval with the same standard error
    df_residual <- meta$fit$fit$df.residual
    expect_equal(c(mix$ci_lower, mix$ci_upper), mix$mixture_effect +
        c(-1, 1) * stats::qt(0.975, df_residual) * mix$mixture_se)
    expect_gt(mix$ci_upper - mix$ci_lower,
        2 * stats::qnorm(0.975) * mix$mixture_se)

    ## The mixture p-value is that of psi, not of the intercept
    psi_index <- match("psi1", names(meta$fit$coef))
    expect_equal(mix$mixture_pvalue, unname(meta$fit$pval[psi_index]))
    expect_false(isTRUE(all.equal(mix$mixture_pvalue,
                                  unname(meta$fit$pval[1]))))

    ## Weights should exist
    all_weights <- c(meta$positive_weights, meta$negative_weights)
    expect_true(length(all_weights) > 0)
})

test_that("the qgcomp mixture p-value uses the residual degrees of freedom", {
    skip_if_not_installed("qgcomp")
    set.seed(42)
    counts <- matrix(rpois(10000, 10), nrow = 100,
        dimnames = list(paste0("G", 1:100), paste0("c", 1:100)))
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = counts),
        colData = S4Vectors::DataFrame(
            donor_id = rep(paste0("D", 1:10), each = 10),
            cell_type = rep(c("Mono", "NK"), 50)))
    exp_mat <- matrix(rnorm(50), nrow = 10,
        dimnames = list(paste0("D", 1:10), c("E1", "E2", "E3", "C1", "C2")))
    scee <- build_scee(sce, exp_mat, sample_col = "donor_id")
    mix <- run_mixture_qgcomp(scee, exposures = c("E1", "E2", "E3"),
        covariates = c("C1", "C2"), celltype = "Mono", q = 2L,
        min_cells = 3L)
    fit <- S4Vectors::metadata(mix)$fit
    ## ten donors, intercept, three exposures and two covariates
    expect_equal(fit$fit$df.residual, 4)
    t_value <- mix$mixture_effect / mix$mixture_se
    expect_equal(mix$mixture_pvalue,
        2 * stats::pt(-abs(t_value), df = 4))
    expect_equal(c(mix$ci_lower, mix$ci_upper), mix$mixture_effect +
        c(-1, 1) * stats::qt(0.975, 4) * mix$mixture_se)
    ## five donors and five coefficients leave no residual degrees of
    ## freedom: the interval and the p-value are NaN, without a warning
    keep <- sce$donor_id %in% paste0("D", 1:5)
    scee5 <- build_scee(sce[, keep], exp_mat[1:5, ], sample_col = "donor_id")
    expect_no_warning(
        mix5 <- run_mixture_qgcomp(scee5, exposures = c("E1", "E2", "E3"),
            covariates = "C1", celltype = "Mono", q = 2L, min_cells = 3L))
    expect_true(all(is.nan(c(mix5$ci_lower, mix5$ci_upper,
        mix5$mixture_pvalue))))
})
