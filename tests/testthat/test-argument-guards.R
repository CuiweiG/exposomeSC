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
    body_expr <- body(method)
    local_fn <- eval(body_expr[[2]][[3]])
    method_args <- setdiff(names(formals(local_fn)), "...")
    set_by_wrapper <- c("x", "exposure", "celltype_col", "covariates",
                        "min_cells", "min_donors", "filter_genes")
    expect_setequal(exposomeSC:::.run_sc_exwas_passthrough,
                    setdiff(method_args, set_by_wrapper))
})

.formula_name_fixture <- function() {
    set.seed(11)
    donors <- sprintf("D%02d", 1:12)
    donor <- rep(donors, each = 20)
    cells <- paste0("c", seq_along(donor))
    counts <- matrix(stats::rpois(50 * length(donor), 5), nrow = 50,
        dimnames = list(paste0("G", 1:50), cells))
    proteins <- matrix(stats::rpois(2 * length(donor), 5), nrow = 2,
        dimnames = list(c("P1", "P2"), cells))
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = counts),
        colData = S4Vectors::DataFrame(donor_id = donor, cell_type = "Mono",
            row.names = cells))
    SingleCellExperiment::altExp(sce, "CITE") <-
        SummarizedExperiment::SummarizedExperiment(
            assays = list(counts = proteins))
    exp_mat <- cbind(E1 = stats::rnorm(12), "PCB-153" = stats::rnorm(12),
        age = stats::rnorm(12, 50, 10), "age-years" = stats::rnorm(12, 50, 10))
    rownames(exp_mat) <- donors
    build_scee(sce, exp_mat, sample_col = "donor_id")
}

test_that("formula names must be syntactic R names", {
    expect_invisible(
        exposomeSC:::.check_formula_names(c("age", "PM2.5"), "covariate"))
    expect_invisible(exposomeSC:::.check_formula_names(NULL, "covariate"))
    expect_error(
        exposomeSC:::.check_formula_names(
            c("age", "1-OHP", "BPA (ng/mL)"), "exposure"),
        paste("Non-syntactic exposure name(s) cannot be used in a model",
            "formula: '1-OHP', 'BPA (ng/mL)'. Rename the variable(s) in",
            "exposureData(x), for example to 'X1.OHP', 'BPA..ng.mL.'."),
        fixed = TRUE)
})

test_that("covariate and exposure names a formula would misread are refused", {
    scee <- .formula_name_fixture()
    covariate_error <- "covariate name\\(s\\) .*'age-years'.*'age.years'"
    expect_error(
        run_sc_exwas(scee, exposure = "E1", covariates = "age-years"),
        covariate_error)
    expect_error(
        run_sc_exwas_pb_offset(scee, exposure = "E1",
            covariates = "age-years"),
        covariate_error)
    expect_error(
        run_sc_mixture(scee, exposures = c("E1", "PCB-153"),
            celltype = "Mono", min_cells = 3L),
        "exposure name\\(s\\) .*'PCB-153'.*'PCB.153'")
    expect_error(
        run_sc_mixture(scee, exposures = c("E1", "age"), celltype = "Mono",
            covariates = "age-years", min_cells = 3L),
        covariate_error)
    ## the exposure of run_sc_exwas() is renamed before the model is built
    result <- suppressWarnings(
        run_sc_exwas(scee, exposure = "PCB-153", covariates = "age"))
    expect_gt(nrow(result), 0L)
    expect_true(all(result$exposure == "PCB-153"))
})

test_that("run_mixture_qgcomp refuses names a formula would misread", {
    skip_if_not_installed("qgcomp")
    scee <- .formula_name_fixture()
    expect_error(
        run_mixture_qgcomp(scee, exposures = c("E1", "PCB-153"),
            celltype = "Mono", min_cells = 3L),
        "exposure name\\(s\\) .*'PCB-153'")
    expect_error(
        run_mixture_qgcomp(scee, exposures = c("E1", "age"),
            celltype = "Mono", covariates = "age-years", min_cells = 3L),
        "covariate name\\(s\\) .*'age-years'")
})

test_that("run_sc_exwas_glmm refuses covariate names a formula would misread", {
    skip_if_not_installed("lme4")
    skip_if_not_installed("lmerTest")
    scee <- .formula_name_fixture()
    expect_error(
        run_sc_exwas_glmm(scee, exposure = "E1", celltype_col = "cell_type",
            covariates = "age-years", target_genes = "G1"),
        "covariate name\\(s\\) .*'age-years'")
})

test_that("run_cell_coupling refuses moderator names a formula would misread", {
    skip_if_not_installed("metafor")
    skip_if_not_installed("ppcor")
    scee <- .formula_name_fixture()
    expect_error(
        run_cell_coupling(scee, genes = "G1", proteins = "P1",
            exposure = "E1", celltype = "Mono", altexp_name = "CITE",
            min_cells = 10L, min_donors = 5L, meta_moderators = "age-years"),
        "moderator name\\(s\\) .*'age-years'")
})

test_that("run_sc_exwas checks its thresholds for every backend", {
    scee <- .guard_fixture()
    for (method in c("edgeR", "DESeq2", "voom_dream")) {
        for (min_cells in list(-1, 2.5, NA_real_)) {
            expect_error(
                run_sc_exwas(scee, exposure = "E1", method = method,
                    min_cells = min_cells),
                "min_cells must be one positive integer")
        }
        expect_error(
            run_sc_exwas(scee, exposure = "E1", method = method,
                min_donors = 3.5),
            "min_donors must be one integer >= 3")
        expect_error(
            run_sc_exwas(scee, exposure = "E1", method = method,
                min_donors = 2L),
            "min_donors must be one integer >= 3")
        expect_error(
            run_sc_exwas(scee, exposure = "E1", method = method,
                filter_genes = NA),
            "filter_genes must be TRUE or FALSE")
    }
    expect_error(
        run_sc_exwas(scee, exposure = "E1", min_group_donors = 2.5),
        "min_group_donors must be one integer >= 2")
    ## a whole number that as.integer() cannot keep
    expect_error(run_sc_exwas(scee, exposure = "E1", min_cells = 3e9),
        "min_cells must be one positive integer")
})

test_that("celltype and exposure must each be a single name", {
    scee <- .formula_name_fixture()
    single <- "'celltype' must be a single name"
    expect_error(run_dose_response(scee, exposure = "E1",
        celltype = c("Mono", "NK"), min_cells = 3L), single)
    expect_error(run_dose_response(scee, exposure = c("E1", "age"),
        celltype = "Mono", min_cells = 3L), "'exposure' must be a single name")
    expect_error(run_sc_mixture(scee, exposures = c("E1", "age"),
        celltype = c("Mono", "NK"), min_cells = 3L), single)
    ## a factor label is taken as its value
    by_factor <- run_dose_response(scee, exposure = "E1",
        celltype = factor("Mono"), min_cells = 3L, max_degree = 2L,
        target_genes = "G1")
    expect_identical(by_factor$celltype, "Mono")
    ## so is a number, for a cluster column of integers
    scee$cluster <- 1L
    by_number <- run_dose_response(scee, exposure = "E1", celltype = 1,
        celltype_col = "cluster", min_cells = 3L, max_degree = 2L,
        target_genes = "G1")
    expect_equal(by_number$linear_coef, by_factor$linear_coef)
    ## a factor exposure is its label, not the column of its integer code
    by_name <- run_sc_exwas(scee, exposure = "age")
    expect_equal(run_sc_exwas(scee, exposure = factor("age"))$pvalue,
        by_name$pvalue)
    expect_equal(
        run_dose_response(scee, exposure = factor("age"), celltype = "Mono",
            min_cells = 3L, max_degree = 2L, target_genes = "G1"),
        run_dose_response(scee, exposure = "age", celltype = "Mono",
            min_cells = 3L, max_degree = 2L, target_genes = "G1"))
    skip_if_not_installed("mgcv")
    expect_error(run_dose_response_gam(scee, exposure = "E1",
        celltype = c("Mono", "NK"), min_cells = 3L), single)
    skip_if_not_installed("qgcomp")
    expect_error(run_mixture_qgcomp(scee, exposures = c("E1", "age"),
        celltype = c("Mono", "NK"), min_cells = 3L), single)
    skip_if_not_installed("metafor")
    skip_if_not_installed("ppcor")
    expect_error(run_cell_coupling(scee, genes = "G1", proteins = "P1",
        min_cells = 10L, min_donors = 5L,
        exposure = "E1", celltype = c("Mono", "NK"), altexp_name = "CITE"),
        single)
})
