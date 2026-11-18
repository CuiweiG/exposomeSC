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

.formula_name_fixture <- function(extra = NULL) {
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
    if (!is.null(extra)) {
        ## one more variable, under a name the test chooses
        exp_mat <- cbind(exp_mat, stats::rnorm(12, 50, 10))
        colnames(exp_mat)[ncol(exp_mat)] <- extra
    }
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

test_that("names that R reads specially in a formula are refused", {
    ## make.names() leaves these names as they are
    special <- c(".", "...", "..1", "..2", "..10", "..01")
    expect_identical(make.names(special), special)
    for (name in special) {
        expect_error(
            exposomeSC:::.check_formula_names(c("age", name), "covariate"),
            paste0("The covariate name(s) '", name, "' cannot be used in a ",
                "model formula, where R reads '.', '...' and names such as ",
                "'..1' specially. Rename the variable(s) in exposureData(x)."),
            fixed = TRUE)
    }
    expect_error(
        exposomeSC:::.check_formula_names(c(".", "age", "..1"), "moderator"),
        "The moderator name(s) '.', '..1' cannot be used", fixed = TRUE)
    ## the names next to them are ordinary names
    expect_invisible(exposomeSC:::.check_formula_names(
        c("..", "....", "...1", ".x", "..a", "..1a", "x.", "PM2.5"),
        "covariate"))
})

test_that("the analyses refuse the names R reads specially in a formula", {
    with_name <- .formula_name_fixture
    refused <- function(role, name) {
        paste0("The ", role, " name(s) '", name, "' cannot be used in a ",
            "model formula")
    }
    for (name in c(".", "...", "..1")) {
        renamed <- with_name(name)
        expect_error(run_sc_exwas(renamed, exposure = "E1", covariates = name),
            refused("covariate", name), fixed = TRUE)
        expect_error(
            run_sc_exwas_pb_offset(renamed, exposure = "E1",
                covariates = name),
            refused("covariate", name), fixed = TRUE)
        expect_error(
            run_sc_mixture(renamed, exposures = c("E1", name),
                celltype = "Mono", min_cells = 3L),
            refused("exposure", name), fixed = TRUE)
        expect_error(
            run_sc_mixture(renamed, exposures = "E1", covariates = name,
                celltype = "Mono", min_cells = 3L),
            refused("covariate", name), fixed = TRUE)
    }
    ## such a name can still be the exposure of run_sc_exwas(), which is
    ## renamed before the model is built
    result <- suppressWarnings(run_sc_exwas(with_name("..1"),
        exposure = "..1"))
    expect_gt(nrow(result), 0L)
    expect_true(all(result$exposure == "..1"))

    skip_if_not_installed("lme4")
    skip_if_not_installed("lmerTest")
    expect_error(
        run_sc_exwas_glmm(with_name("."), exposure = "E1", covariates = ".",
            target_genes = "G1"),
        refused("covariate", "."), fixed = TRUE)
    skip_if_not_installed("qgcomp")
    for (bootstrap in c(FALSE, TRUE)) {
        expect_error(
            run_mixture_qgcomp(with_name("."), exposures = "E1",
                covariates = ".", celltype = "Mono", min_cells = 3L,
                bootstrap = bootstrap),
            refused("covariate", "."), fixed = TRUE)
    }
    expect_error(
        run_mixture_qgcomp(with_name("..."), exposures = c("E1", "..."),
            celltype = "Mono", min_cells = 3L),
        refused("exposure", "..."), fixed = TRUE)
    skip_if_not_installed("metafor")
    skip_if_not_installed("ppcor")
    for (name in c(".", "...", "..1")) {
        expect_error(
            run_cell_coupling(with_name(name), genes = "G1", proteins = "P1",
                exposure = "E1", celltype = "Mono", altexp_name = "CITE",
                min_cells = 10L, min_donors = 5L, meta_moderators = name),
            refused("moderator", name), fixed = TRUE)
    }
    ## a formula can name them too; ~ . is refused as a formula
    expect_error(
        run_cell_coupling(with_name("..1"), genes = "G1", proteins = "P1",
            exposure = "E1", celltype = "Mono", altexp_name = "CITE",
            min_cells = 10L, min_donors = 5L, meta_moderators = ~ ..1),
        refused("moderator", "..1"), fixed = TRUE)
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

test_that("an exposure cannot also be a covariate", {
    scee <- .formula_name_fixture()
    overlap <- "must not also be listed as a covariate"
    expect_error(run_sc_mixture(scee, exposures = c("E1", "age"),
        celltype = "Mono", covariates = "E1", min_cells = 3L),
        paste0(overlap, ": E1"))
    expect_error(run_sc_mixture(scee, exposures = c("E1", "E1"),
        celltype = "Mono", min_cells = 3L), "duplicated")
    ## two cell types, so that without the check the F tests would run
    donor <- rep(sprintf("D%02d", 1:10), each = 6)
    counts <- matrix(stats::rpois(20 * length(donor), 8), nrow = 20,
        dimnames = list(paste0("G", 1:20), paste0("c", seq_along(donor))))
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = counts),
        colData = S4Vectors::DataFrame(donor_id = donor,
            cell_type = rep(c("Mono", "NK"), length.out = length(donor))))
    exposures <- cbind(E1 = stats::rnorm(10), age = stats::rnorm(10, 50, 10))
    rownames(exposures) <- sprintf("D%02d", 1:10)
    two_types <- build_scee(sce, exposures, sample_col = "donor_id")
    expect_error(run_interaction_test(two_types, exposure = "E1",
        covariates = "E1", min_cells = 3L, min_donors = 3L), overlap)
    skip_if_not_installed("qgcomp")
    expect_error(run_mixture_qgcomp(scee, exposures = c("E1", "age"),
        celltype = "Mono", covariates = "E1", min_cells = 3L), overlap)
})

test_that("target_genes missing from the data are named in a warning", {
    ## beyond ten names the warning counts the others
    expect_warning(
        exposomeSC:::.match_target_genes(c("G1", paste0("X", 1:12)), "G1"),
        "X10 and 2 more.", fixed = TRUE)
    scee <- .formula_name_fixture()
    absent <- "target_genes not found in the data are not tested: GX"
    expect_warning(
        dose <- run_dose_response(scee, exposure = "E1", celltype = "Mono",
            min_cells = 3L, max_degree = 2L, target_genes = c("G1", "GX")),
        absent)
    expect_identical(dose$gene, "G1")
    expect_warning(run_sc_mixture(scee, exposures = c("E1", "age"),
        celltype = "Mono", min_cells = 3L,
        target_genes = c("G1", "G2", "GX")), absent)
    expect_warning(run_sc_exwas_pb_offset(scee, exposure = "E1",
        target_genes = c("G1", "GX")), absent)
    expect_error(suppressWarnings(run_sc_exwas_pb_offset(scee,
        exposure = "E1", target_genes = "GX")), "No target genes present")
    set.seed(3)
    donor <- rep(sprintf("D%02d", 1:10), each = 6)
    counts <- matrix(stats::rpois(20 * length(donor), 8), nrow = 20,
        dimnames = list(paste0("G", 1:20), paste0("c", seq_along(donor))))
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = counts),
        colData = S4Vectors::DataFrame(donor_id = donor,
            cell_type = rep(c("Mono", "NK"), length.out = length(donor))))
    two_types <- build_scee(sce, matrix(stats::rnorm(10), ncol = 1,
        dimnames = list(sprintf("D%02d", 1:10), "E1")),
        sample_col = "donor_id")
    expect_warning(run_interaction_test(two_types, exposure = "E1",
        min_cells = 3L, min_donors = 3L, target_genes = c("G1", "GX")),
        absent)
    skip_if_not_installed("mgcv")
    expect_warning(run_dose_response_gam(scee, exposure = "E1",
        celltype = "Mono", min_cells = 3L, target_genes = c("G1", "GX")),
        absent)
    skip_if_not_installed("qgcomp")
    expect_warning(run_mixture_qgcomp(scee, exposures = c("E1", "age"),
        celltype = "Mono", min_cells = 3L, target_genes = c("G1", "GX")),
        absent)
    skip_if_not_installed("lme4")
    skip_if_not_installed("lmerTest")
    expect_warning(suppressMessages(run_sc_exwas_glmm(scee, exposure = "E1",
        target_genes = c("G1", "GX"))), absent)
})

test_that("run_sc_exwas refuses the exposure as its own covariate", {
    scee <- .guard_fixture()
    expect_error(run_sc_exwas(scee, exposure = "E1", covariates = "E1"),
        "The exposure must not also be listed as a covariate")
    expect_error(run_sc_exwas(scee, exposure = "E1",
        covariates = c("age", "E1")),
        "The exposure must not also be listed as a covariate")
    expect_error(run_sc_exwas_pb_offset(scee, exposure = "E1",
        covariates = "E1"),
        "The exposure must not also be listed as a covariate")
    skip_if_not_installed("lme4")
    skip_if_not_installed("lmerTest")
    expect_error(run_sc_exwas_glmm(scee, exposure = "E1", covariates = "E1",
        family = "nbinom", target_genes = "G1"),
        "The exposure must not also be listed as a covariate")
})
