.make_edger_scee <- function() {
    set.seed(20260712)
    donors <- paste0("D", seq_len(20))
    exposure <- rep(c(0, 1), each = 10)
    donor <- rep(donors, each = 40)
    celltype <- rep(rep(c("Mono", "NK"), each = 20), 20)
    genes <- paste0("G", seq_len(80))
    counts <- matrix(
        stats::rnbinom(length(genes) * length(donor), mu = 8, size = 5),
        nrow = length(genes),
        dimnames = list(genes, paste0("cell", seq_along(donor)))
    )
    signal_cells <- celltype == "Mono" & exposure[match(donor, donors)] == 1
    counts["G1", signal_cells] <- counts["G1", signal_cells] * 3L
    counts["G80", ] <- 0L
    sparse_cells <- donor %in% donors[11:14]
    counts["G80", sparse_cells] <- stats::rpois(sum(sparse_cells), 20)
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = counts),
        colData = S4Vectors::DataFrame(
            donor_id = donor,
            cell_type = celltype
        )
    )
    exp_data <- cbind(
        exposure = exposure,
        age_scaled = as.numeric(scale(seq_along(donors))),
        study_b = rep(c(0, 1), 10)
    )
    rownames(exp_data) <- donors
    build_scee(sce, exp_data, sample_col = "donor_id")
}

test_that("edgeR sc-ExWAS returns globally adjusted, auditable results", {
    skip_if_not_installed("edgeR")
    scee <- .make_edger_scee()
    result <- run_sc_exwas(
        scee,
        exposure = "exposure",
        celltype_col = "cell_type",
        covariates = c("age_scaled", "study_b"),
        min_cells = 10L,
        min_donors = 15L,
        min_group_donors = 5L,
        robust = FALSE
    )

    expect_s4_class(result, "DataFrame")
    expect_gt(nrow(result), 0L)
    expect_setequal(unique(result$celltype), c("Mono", "NK"))
    expect_true(all(result$method == "edgeR_QL"))
    expect_true(all(result$n_donors == 20L))
    expect_true(all(result$n_unexposed == 10L))
    expect_true(all(result$n_exposed == 10L))
    expect_type(result$se, "double")
    expect_true(all(is.na(result$se)))
    expect_false(any(is.nan(result$se)))
    expect_identical(
        S4Vectors::metadata(result)$parameters$schema_version,
        "exposomeSC_sc_exwas_v3"
    )
    expect_identical(
        S4Vectors::metadata(result)$parameters$se_method,
        "not_available_edgeR_QL"
    )
    expect_identical(
        S4Vectors::metadata(result)$parameters$statistic_type,
        "signed_sqrt_qlf"
    )
    nonzero <- result$statistic != 0
    expect_equal(
        sign(result$statistic[nonzero]),
        sign(result$log2FC[nonzero])
    )
    expect_equal(
        result$padj_global,
        stats::p.adjust(result$pvalue, method = "BH"),
        tolerance = 1e-12
    )
    expect_false(any(result$pvalue == 0, na.rm = TRUE))

    parameters <- S4Vectors::metadata(result)$parameters
    expect_equal(parameters$exposure, "exposure")
    expect_equal(parameters$covariates, c("age_scaled", "study_b"))
    expect_equal(parameters$filter_method, "fixed_support")
    expect_equal(
        parameters$multiplicity,
        "BH across all tested gene-by-cell-type hypotheses"
    )
})

test_that("edgeR sc-ExWAS retains full-transcriptome library sizes", {
    skip_if_not_installed("edgeR")
    scee <- .make_edger_scee()
    result <- run_sc_exwas(
        scee,
        exposure = "exposure",
        celltype_col = "cell_type",
        celltypes = "Mono",
        min_cells = 10L,
        min_donors = 15L,
        min_group_donors = 5L,
        robust = FALSE
    )

    diagnostics <- S4Vectors::metadata(result)$diagnostics$Mono
    counts <- SummarizedExperiment::assay(scee, "counts")
    cd <- SummarizedExperiment::colData(scee)
    expected <- vapply(diagnostics$donor_ids, function(donor) {
        index <- cd$donor_id == donor & cd$cell_type == "Mono"
        sum(counts[, index, drop = FALSE])
    }, numeric(1))
    expect_equal(unname(diagnostics$full_library_size), unname(expected))
    expect_equal(diagnostics$n_genes_input, nrow(scee))
    expect_equal(diagnostics$filter_method, "fixed_support")
    expect_equal(diagnostics$filter_min_cpm, 1)
    expect_equal(diagnostics$filter_min_donors, 5L)
    expect_equal(diagnostics$n_genes_fit, diagnostics$n_genes_tested)
    expect_false("G80" %in% diagnostics$test_feature_ids)
})

test_that("test_features changes only the reported hypothesis family", {
    skip_if_not_installed("edgeR")
    scee <- .make_edger_scee()
    test_family <- c("G1", "G2", "G3")
    arguments <- list(
        x = scee,
        exposure = "exposure",
        celltype_col = "cell_type",
        celltypes = c("Mono", "NK"),
        min_cells = 10L,
        min_donors = 15L,
        min_group_donors = 5L,
        robust = FALSE
    )
    full <- do.call(run_sc_exwas, arguments)
    focused <- do.call(
        run_sc_exwas,
        c(arguments, list(test_features = test_family))
    )

    expect_setequal(unique(focused$gene), test_family)
    expect_equal(nrow(focused), length(test_family) * 2L)

    full_data <- as.data.frame(full)
    focused_data <- as.data.frame(focused)
    matched <- merge(
        focused_data,
        full_data,
        by = c("gene", "celltype"),
        suffixes = c("_focused", "_full")
    )
    expect_equal(
        matched$log2FC_focused,
        matched$log2FC_full,
        tolerance = 1e-12
    )
    expect_equal(
        matched$statistic_focused,
        matched$statistic_full,
        tolerance = 1e-12
    )
    expect_equal(
        matched$pvalue_focused,
        matched$pvalue_full,
        tolerance = 1e-12
    )
    expect_equal(
        matched$logCPM_focused,
        matched$logCPM_full,
        tolerance = 1e-12
    )

    full_diagnostics <- S4Vectors::metadata(full)$diagnostics
    focused_diagnostics <- S4Vectors::metadata(focused)$diagnostics
    for (celltype in c("Mono", "NK")) {
        expect_equal(
            focused_diagnostics[[celltype]]$full_library_size,
            full_diagnostics[[celltype]]$full_library_size
        )
        expect_equal(
            focused_diagnostics[[celltype]]$normalisation_factors,
            full_diagnostics[[celltype]]$normalisation_factors
        )
        expect_equal(
            focused_diagnostics[[celltype]]$n_genes_fit,
            full_diagnostics[[celltype]]$n_genes_fit
        )
        expect_equal(
            focused_diagnostics[[celltype]]$n_genes_tested,
            length(test_family)
        )
        expect_setequal(
            focused_diagnostics[[celltype]]$test_feature_ids,
            test_family
        )

        celltype_rows <- focused$celltype == celltype
        expect_equal(
            focused$padj[celltype_rows],
            stats::p.adjust(
                focused$pvalue[celltype_rows],
                method = "BH"
            ),
            tolerance = 1e-12
        )
    }
    expect_equal(
        focused$padj_global,
        stats::p.adjust(focused$pvalue, method = "BH"),
        tolerance = 1e-12
    )
    parameters <- S4Vectors::metadata(focused)$parameters
    expect_equal(parameters$test_features, test_family)
    expect_equal(parameters$n_test_features_requested, length(test_family))
})

test_that("test_features never silently drops a requested feature", {
    skip_if_not_installed("edgeR")
    expect_error(
        run_sc_exwas(
            .make_edger_scee(),
            exposure = "exposure",
            celltype_col = "cell_type",
            celltypes = "Mono",
            min_cells = 10L,
            min_donors = 15L,
            min_group_donors = 5L,
            test_features = c("G1", "G80"),
            robust = FALSE
        ),
        "did not pass expression filtering.*G80"
    )
})

test_that("fixed-support filtering is invariant to exposure labels", {
    skip_if_not_installed("edgeR")
    scee <- .make_edger_scee()
    permuted <- scee
    permuted_exposure <- exposureData(permuted)
    permuted_exposure[, "exposure"] <- rev(
        permuted_exposure[, "exposure"]
    )
    exposureData(permuted) <- permuted_exposure

    run_fixed_support <- function(object) {
        run_sc_exwas(
            object,
            exposure = "exposure",
            celltype_col = "cell_type",
            celltypes = "Mono",
            min_cells = 10L,
            min_donors = 15L,
            min_group_donors = 5L,
            filter_method = "fixed_support",
            filter_min_cpm = 1,
            filter_min_donors = 5L,
            robust = FALSE
        )
    }
    observed <- run_fixed_support(scee)
    relabelled <- run_fixed_support(permuted)
    observed_diagnostics <- S4Vectors::metadata(observed)$diagnostics$Mono
    relabelled_diagnostics <-
        S4Vectors::metadata(relabelled)$diagnostics$Mono

    expect_equal(
        observed_diagnostics$test_feature_ids,
        relabelled_diagnostics$test_feature_ids
    )
    expect_equal(
        observed_diagnostics$n_genes_fit,
        relabelled_diagnostics$n_genes_fit
    )
    expect_equal(
        observed_diagnostics$full_library_size,
        relabelled_diagnostics$full_library_size
    )
    expect_false("G80" %in% observed_diagnostics$test_feature_ids)
})

test_that("filterByExpr remains available as a compatibility option", {
    skip_if_not_installed("edgeR")
    result <- run_sc_exwas(
        .make_edger_scee(),
        exposure = "exposure",
        celltype_col = "cell_type",
        celltypes = "Mono",
        min_cells = 10L,
        min_donors = 15L,
        min_group_donors = 5L,
        filter_method = "filterByExpr",
        robust = FALSE
    )

    expect_gt(nrow(result), 0L)
    expect_equal(
        S4Vectors::metadata(result)$diagnostics$Mono$filter_method,
        "filterByExpr"
    )
})

test_that("edgeR sc-ExWAS validates design inputs", {
    skip_if_not_installed("edgeR")
    scee <- .make_edger_scee()
    expect_error(
        run_sc_exwas(
            scee,
            exposure = "missing",
            celltype_col = "cell_type"
        ),
        "not found"
    )
    expect_error(
        run_sc_exwas(
            scee,
            exposure = "exposure",
            celltype_col = "cell_type",
            min_group_donors = 1L
        ),
        "min_group_donors"
    )
    expect_error(
        run_sc_exwas(
            scee,
            exposure = "exposure",
            celltype_col = "cell_type",
            test_features = "absent"
        ),
        "not found"
    )
    expect_error(
        run_sc_exwas(
            scee,
            exposure = "exposure",
            celltype_col = "cell_type",
            test_features = c("G1", "G1")
        ),
        "unique"
    )
    expect_error(
        run_sc_exwas(
            scee,
            exposure = "exposure",
            celltype_col = "cell_type",
            filter_min_cpm = 0
        ),
        "filter_min_cpm"
    )
    expect_error(
        run_sc_exwas(
            scee,
            exposure = "exposure",
            celltype_col = "cell_type",
            filter_min_donors = 1.5
        ),
        "filter_min_donors"
    )
})

test_that("empty edgeR results retain the uncertainty schema", {
    skip_if_not_installed("edgeR")
    result <- suppressWarnings(run_sc_exwas(
        .make_edger_scee(),
        exposure = "exposure",
        celltype_col = "cell_type",
        min_cells = 10L,
        min_donors = 100L,
        min_group_donors = 5L,
        robust = FALSE
    ))

    expect_equal(nrow(result), 0L)
    parameters <- S4Vectors::metadata(result)$parameters
    expect_identical(
        parameters$schema_version,
        "exposomeSC_sc_exwas_v3"
    )
    expect_identical(parameters$se_method, "not_available_edgeR_QL")
    expect_identical(parameters$statistic_type, "signed_sqrt_qlf")
    expect_type(result$se, "double")
    expect_identical(colnames(result), exposomeSC:::.sc_exwas_columns)
})

test_that("a covariate named 'exposure' is refused", {
    skip_if_not_installed("edgeR")
    scee <- .make_edger_scee()
    expect_error(
        run_sc_exwas(scee, exposure = "age_scaled", celltype_col = "cell_type",
            covariates = "exposure", min_donors = 15L),
        "cannot be called 'exposure'"
    )
})

test_that("a collinear covariate is dropped with a warning and recorded", {
    skip_if_not_installed("edgeR")
    scee <- .make_edger_scee()
    exp_data <- exposureData(scee)
    exp_data <- cbind(exp_data, exposure_copy = exp_data[, "exposure"])
    scee2 <- build_scee(methods::as(scee, "SingleCellExperiment"), exp_data,
                        sample_col = "donor_id")
    warnings <- character()
    result <- withCallingHandlers(
        run_sc_exwas(scee2, exposure = "exposure",
            celltype_col = "cell_type",
            covariates = "exposure_copy",
            min_donors = 15L, robust = FALSE),
        warning = function(w) {
            warnings <<- c(warnings, conditionMessage(w))
            invokeRestart("muffleWarning")
        }
    )
    ## one warning per cell type, naming the dropped column
    expect_length(warnings, 2L)
    expect_true(all(grepl("exposure_copy", warnings)))
    expect_true(any(startsWith(warnings, "Mono:")))
    expect_true(any(startsWith(warnings, "NK:")))
    diag <- S4Vectors::metadata(result)$diagnostics
    expect_identical(diag$Mono$dropped_design_columns, "exposure_copy")
    expect_identical(diag$Mono$design_columns, c("(Intercept)", "exposure"))
    expect_identical(diag$NK$dropped_design_columns, "exposure_copy")
})

test_that("a cell type that is skipped gets no collinear-covariate warning", {
    skip_if_not_installed("edgeR")
    scee <- .make_edger_scee()
    exp_data <- exposureData(scee)
    exp_data <- cbind(exp_data, exposure_copy = exp_data[, "exposure"])
    sce <- methods::as(scee, "SingleCellExperiment")
    cells <- SummarizedExperiment::colData(sce)
    keep <- cells$cell_type == "Mono" |
        cells$donor_id %in% c("D9", "D10", "D11", "D12")
    scee2 <- build_scee(sce[, keep], exp_data, sample_col = "donor_id")
    warnings <- character()
    withCallingHandlers(
        run_sc_exwas(scee2, exposure = "exposure",
            celltype_col = "cell_type", covariates = "exposure_copy",
            min_donors = 3L, min_group_donors = 2L, robust = FALSE),
        warning = function(w) {
            warnings <<- c(warnings, conditionMessage(w))
            invokeRestart("muffleWarning")
        }
    )
    ## NK keeps four donors, too few residual degrees of freedom to fit
    expect_length(warnings, 2L)
    expect_true(any(startsWith(warnings, "Mono: covariate column(s)")))
    expect_true(any(startsWith(warnings, "Skipping NK: insufficient")))
})

test_that("the residual df diagnostic uses edgeR's adjusted residual df", {
    skip_if_not_installed("edgeR")
    scee <- .make_edger_scee()
    ## G80 is zero in 16 of the 20 donors; keep it in the fit
    result <- run_sc_exwas(scee, exposure = "exposure",
        celltype_col = "cell_type",
        covariates = c("age_scaled", "study_b"),
        min_donors = 15L, filter_min_donors = 2L,
        robust = FALSE)
    diag <- S4Vectors::metadata(result)$diagnostics$Mono
    nominal <- diag$n_donors - length(diag$design_columns)
    expect_true(is.finite(diag$residual_df))
    expect_gt(diag$residual_df, 0)
    expect_lt(diag$residual_df, nominal - 0.5)
})

## ---- the steps of the edgeR backend ------------------------------------

.edger_step_fixture <- function() {
    scee <- .make_edger_scee()
    list(
        scee = scee,
        inputs = exposomeSC:::.sc_exwas_inputs(scee, "exposure", NULL,
            "cell_type", NULL),
        formula = stats::reformulate("exposure"),
        filter = list(genes = TRUE, method = "fixed_support", min_cpm = 1,
            min_donors = 5L, min_count = 10L, min_total_count = 15L)
    )
}

.edger_step_fit <- function(fx, min_cells = 10L) {
    exposomeSC:::.edger_fit_celltype(fx$inputs, "Mono", "exposure", NULL,
        fx$formula, min_cells, 5L, 5L, fx$filter, TRUE)
}

test_that("the pseudobulk step returns the complete donors or warns", {
    fx <- .edger_step_fixture()
    step <- exposomeSC:::.edger_pseudobulk
    stratum <- step(fx$inputs, "Mono", "exposure", NULL, 10L, 5L)
    expect_named(stratum,
        c("donors", "design_info", "pseudobulk", "donor_cells"))
    expect_setequal(stratum$donors, paste0("D", seq_len(20)))
    expect_identical(colnames(stratum$pseudobulk), stratum$donors)
    expect_identical(dim(stratum$pseudobulk), c(80L, 20L))
    expect_true(all(stratum$donor_cells == 20L))
    expect_warning(
        expect_null(step(fx$inputs, "Mono", "exposure", NULL, 50L, 5L)),
        "Skipping Mono: no donors with >= 50 cells.", fixed = TRUE)
    expect_warning(
        expect_null(step(fx$inputs, "Mono", "exposure", NULL, 10L, 21L)),
        "Skipping Mono: only 20 complete donors with >= 10 cells.",
        fixed = TRUE)
    ## a donor without an exposure value is left out
    inputs <- fx$inputs
    inputs$exp_data["D3", "exposure"] <- NA
    fewer <- step(inputs, "Mono", "exposure", NULL, 10L, 5L)
    expect_identical(fewer$donors, setdiff(stratum$donors, "D3"))
    expect_identical(colnames(fewer$pseudobulk), fewer$donors)
    expect_identical(rownames(fewer$design_info$data), fewer$donors)
    expect_length(fewer$donor_cells, 19L)
    expect_warning(
        expect_null(step(inputs, "Mono", "exposure", NULL, 10L, 20L)),
        "Skipping Mono: only 19 complete donors with >= 10 cells.",
        fixed = TRUE)
})

test_that("the design step drops collinear covariates or warns", {
    fx <- .edger_step_fixture()
    step <- exposomeSC:::.edger_design
    design_data <- exposomeSC:::.sc_exwas_design_data
    donors <- paste0("D", seq_len(20))
    exp_data <- cbind(fx$inputs$exp_data,
        exposure_copy = fx$inputs$exp_data[, "exposure"])
    info <- design_data(exp_data, donors, "exposure", NULL)
    design <- step(fx$formula, info, "Mono", 5L)
    expect_identical(colnames(design$design), c("(Intercept)", "exposure"))
    expect_identical(design$dropped_design, character())
    ## a covariate that copies the exposure is left out and named
    covariates <- c("age_scaled", "exposure_copy")
    adjusted <- step(stats::reformulate(c("exposure", covariates)),
        design_data(exp_data, donors, "exposure", covariates), "Mono", 5L)
    expect_identical(colnames(adjusted$design),
        c("(Intercept)", "exposure", "age_scaled"))
    expect_identical(adjusted$dropped_design, "exposure_copy")
    expect_warning(expect_null(step(fx$formula, info, "Mono", 11L)),
        paste("Skipping Mono: binary exposure group has fewer than 11",
            "complete donors."), fixed = TRUE)
    unexposed <- design_data(exp_data, donors[1:10], "exposure", NULL)
    expect_warning(expect_null(step(fx$formula, unexposed, "NK", 5L)),
        "Skipping NK: the exposure is not identifiable.", fixed = TRUE)
    four <- design_data(exp_data, donors[9:12], "exposure", NULL)
    expect_warning(expect_null(step(fx$formula, four, "NK", 2L)),
        "Skipping NK: insufficient residual degrees of freedom.",
        fixed = TRUE)
})

test_that("the fit step filters, normalises and fits, or warns", {
    fx <- .edger_step_fixture()
    stratum <- exposomeSC:::.edger_pseudobulk(fx$inputs, "Mono", "exposure",
        NULL, 10L, 5L)
    design <- exposomeSC:::.edger_design(fx$formula, stratum$design_info,
        "Mono", 5L)
    filter_step <- exposomeSC:::.edger_filter
    filtered <- filter_step(stratum$pseudobulk, design$design, fx$filter,
        "Mono")
    ## library sizes come from the complete transcriptome, before the filter
    expect_identical(filtered$full_library_size,
        colSums(stratum$pseudobulk))
    expect_identical(unname(filtered$dge$samples$lib.size),
        unname(colSums(stratum$pseudobulk)))
    expect_identical(setdiff(rownames(stratum$pseudobulk),
        rownames(filtered$dge)), "G80")
    strict <- fx$filter
    strict$min_total_count <- 1e9
    expect_warning(
        expect_null(filter_step(stratum$pseudobulk, design$design, strict,
            "Mono")),
        "Skipping Mono: no gene passed expression filtering.", fixed = TRUE)

    fit_step <- exposomeSC:::.edger_fit
    fit <- fit_step(stratum$pseudobulk, design, fx$filter, TRUE, "Mono")
    expect_named(fit,
        c("fit", "table", "full_library_size", "normalisation_factors"))
    expect_identical(rownames(fit$table), rownames(filtered$dge))
    expect_identical(fit$normalisation_factors,
        filtered$dge$samples$norm.factors)
    ## the factors are TMM factors of the kept genes on the library sizes of
    ## the complete transcriptome
    kept <- stratum$pseudobulk[rownames(fit$table), , drop = FALSE]
    expect_equal(fit$normalisation_factors,
        unname(edgeR::normLibSizes(kept,
            lib.size = colSums(stratum$pseudobulk), method = "TMM")))
    ## an error inside edgeR becomes a warning: no column of this design is
    ## called "exposure", so the test of that coefficient stops
    renamed <- design
    colnames(renamed$design)[2] <- "other"
    expect_warning(
        expect_null(fit_step(stratum$pseudobulk, renamed, fx$filter, TRUE,
            "Mono")),
        "edgeR failed for Mono: ", fixed = TRUE)
    ## a dropped covariate is reported once the cell type has been fitted
    warnings <- character()
    collect <- function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
    }
    dropped <- design
    dropped$dropped_design <- "exposure_copy"
    refit <- withCallingHandlers(
        fit_step(stratum$pseudobulk, dropped, fx$filter, TRUE, "Mono"),
        warning = collect)
    expect_identical(warnings, paste("Mono: covariate column(s)",
        "exposure_copy dropped from the design as collinear with the other",
        "terms; the test is not adjusted for them."))
    expect_identical(refit$table, fit$table)
    ## and not when no gene passes the filter or the fit fails
    renamed$dropped_design <- "exposure_copy"
    warnings <- character()
    withCallingHandlers({
        fit_step(stratum$pseudobulk, dropped, strict, TRUE, "Mono")
        fit_step(stratum$pseudobulk, renamed, fx$filter, TRUE, "Mono")
    }, warning = collect)
    expect_length(warnings, 2L)
    expect_identical(warnings[1],
        "Skipping Mono: no gene passed expression filtering.")
    expect_true(startsWith(warnings[2], "edgeR failed for Mono: "))
})

test_that("the extraction step keeps the tested genes in fit order", {
    fx <- .edger_step_fixture()
    fitted <- .edger_step_fit(fx)
    extract <- exposomeSC:::.edger_extract
    expect_identical(extract(fitted$table, NULL), fitted$table)
    ## G1 goes up with the exposure; the second gene goes down, so that the
    ## sign of the statistic is tested
    down <- rownames(fitted$table)[fitted$table$logFC < 0][1]
    some <- extract(fitted$table, c(down, "G1"))
    expect_identical(rownames(some), c("G1", down))
    as_table <- exposomeSC:::.edger_table
    table <- as_table(some, fitted, "Mono", "exposure", TRUE)
    expect_identical(table$gene, c("G1", down))
    expect_identical(table$log2FC, some$logFC)
    expect_identical(sign(table$log2FC), c(1, -1))
    expect_identical(table$statistic, sign(some$logFC) * sqrt(some$F))
    expect_identical(table$pvalue, some$PValue)
    expect_identical(table$logCPM, some$logCPM)
    expect_true(all(is.na(table$se)))
    expect_identical(unique(table$n_donors), 20L)
    expect_identical(unique(table$method), "edgeR_robust_QL")
    plain <- as_table(some, fitted, "Mono", "exposure", FALSE)
    expect_identical(unique(plain$method), "edgeR_QL")
    ## group sizes and cells per donor that can be told apart
    uneven <- fitted
    uneven$design_info$n_unexposed <- 9L
    uneven$design_info$n_exposed <- 11L
    uneven$donor_cells <- c(12L, 30L, rep(20L, 18))
    odd <- as_table(some, uneven, "Mono", "exposure", TRUE)
    expect_identical(unique(odd$n_unexposed), 9L)
    expect_identical(unique(odd$n_exposed), 11L)
    expect_identical(unique(odd$min_cells), 12L)
    expect_identical(unique(odd$median_cells), 20)
})

test_that("the steps give what run_sc_exwas() returns", {
    fx <- .edger_step_fixture()
    result <- run_sc_exwas(fx$scee, exposure = "exposure",
        celltype_col = "cell_type", celltypes = "Mono")
    fitted <- .edger_step_fit(fx)
    rows <- exposomeSC:::.edger_extract(fitted$table, NULL)
    expect_identical(S4Vectors::metadata(result)$diagnostics$Mono,
        exposomeSC:::.edger_diagnostics(fitted, rows, NULL, fx$filter))
    parameters <- exposomeSC:::.edger_parameters("exposure", NULL, 10L, 5L,
        5L, fx$filter, NULL, TRUE)
    expect_identical(S4Vectors::metadata(result)$parameters, parameters)
    ## settings that are NULL keep their place among the parameters
    expect_true(all(c("covariates", "test_features",
        "n_test_features_requested") %in% names(parameters)))
    table <- exposomeSC:::.edger_table(rows, fitted, "Mono", "exposure",
        TRUE)
    by_gene <- as.data.frame(result)
    by_gene <- by_gene[match(table$gene, by_gene$gene), ]
    expect_identical(by_gene$pvalue, table$pvalue)
    expect_identical(by_gene$log2FC, table$log2FC)
    ## settings other than the defaults reach the steps in their places
    features <- c("G3", "G1")
    focused <- run_sc_exwas(fx$scee, exposure = "exposure",
        celltype_col = "cell_type", celltypes = "Mono", min_cells = 12L,
        min_donors = 15L, min_group_donors = 6L, test_features = features,
        robust = FALSE)
    refit <- exposomeSC:::.edger_fit_celltype(fx$inputs, "Mono", "exposure",
        NULL, fx$formula, 12L, 15L, 6L, fx$filter, FALSE)
    some <- exposomeSC:::.edger_extract(refit$table, features)
    expect_identical(S4Vectors::metadata(focused)$diagnostics$Mono,
        exposomeSC:::.edger_diagnostics(refit, some, features, fx$filter))
    expect_identical(S4Vectors::metadata(focused)$parameters,
        exposomeSC:::.edger_parameters("exposure", NULL, 12L, 15L, 6L,
            fx$filter, features, FALSE))
    ## a cell type that one of the steps leaves out gives NULL
    expect_warning(expect_null(.edger_step_fit(fx, min_cells = 50L)),
        "Skipping Mono: no donors with >= 50 cells.", fixed = TRUE)
})
