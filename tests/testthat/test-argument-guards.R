# tests/testthat/test-argument-guards.R

.guard_fixture <- function() {
    set.seed(7)
    donors <- sprintf("D%02d", 1:12)
    donor <- rep(donors, each = 20)
    counts <- matrix(stats::rpois(200 * length(donor), 5), nrow = 200,
        dimnames = list(paste0("G", 1:200), paste0("c", seq_along(donor))))
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = counts),
        colData = S4Vectors::DataFrame(donor_id = donor, cell_type = "Mono"))
    exp_mat <- cbind(E1 = stats::rnorm(12), age = stats::rnorm(12, 50, 10))
    rownames(exp_mat) <- donors
    build_scee(sce, exp_mat, sample_col = "donor_id")
}

test_that("run_sc_exwas rejects arguments it does not use", {
    scee <- .guard_fixture()
    expect_error(
        run_sc_exwas(scee, exposure = "E1", celltype_col = "cell_type",
                     covariats = "age"),
        "Unused argument\\(s\\) in run_sc_exwas\\(\\): covariats")
})

test_that("run_multi_exwas rejects arguments run_sc_exwas does not accept", {
    scee <- .guard_fixture()
    expect_error(
        run_multi_exwas(scee, exposures = "E1", covariats = "age"),
        "Unused argument\\(s\\) in run_multi_exwas\\(\\): covariats")
})

test_that("the run_multi_exwas pass-through list matches run_sc_exwas", {
    method <- methods::getMethod("run_sc_exwas",
                                 "SingleCellExposomeExperiment")
    body_expr <- body(method@.Data)
    local_fn <- eval(body_expr[[2]][[3]])
    method_args <- setdiff(names(formals(local_fn)), "...")
    set_by_wrapper <- c("x", "exposure", "celltype_col", "sample_col",
                        "covariates", "min_cells", "min_donors",
                        "filter_genes")
    expect_setequal(exposomeSC:::.run_sc_exwas_passthrough,
                    setdiff(method_args, set_by_wrapper))
})
