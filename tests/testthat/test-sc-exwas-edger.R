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
        "exposomeSC_sc_exwas_edger_v2"
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

test_that("empty edgeR results retain the v2 uncertainty schema", {
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
        "exposomeSC_sc_exwas_edger_v2"
    )
    expect_identical(parameters$se_method, "not_available_edgeR_QL")
    expect_identical(parameters$statistic_type, "signed_sqrt_qlf")
    expect_type(result$se, "double")
})
