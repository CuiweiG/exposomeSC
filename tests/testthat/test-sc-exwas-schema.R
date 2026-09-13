.make_backend_schema_scee <- function() {
    set.seed(20260712)
    donors <- paste0("D", seq_len(8L))
    cell_donors <- rep(donors, each = 10L)
    genes <- paste0("G", seq_len(20L))
    counts <- matrix(
        stats::rnbinom(
            length(genes) * length(cell_donors),
            mu = 20,
            size = 10
        ),
        nrow = length(genes),
        dimnames = list(genes, paste0("cell", seq_along(cell_donors)))
    )
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = counts),
        colData = S4Vectors::DataFrame(
            donor_id = cell_donors,
            cell_type = "Mono"
        )
    )
    exposure <- matrix(
        rep(c(0, 1), each = 4L),
        ncol = 1L,
        dimnames = list(donors, "smoking")
    )
    build_scee(sce, exposure, sample_col = "donor_id")
}

test_that("backend result construction follows the current schema", {
    result <- .sc_exwas_result_frame(
        gene = c("G1", "G2"),
        celltype = "Mono",
        log2FC = c(0.5, -0.2),
        se = c(0.1, 0.2),
        statistic = c(5, -1),
        pvalue = c(0, 0.4),
        exposure = "smoking",
        n_donors = 20L,
        n_unexposed = 10L,
        n_exposed = 10L,
        min_cells = 12L,
        median_cells = 30,
        method = "test_backend"
    )

    expect_identical(
        colnames(result),
        c(
            "gene", "celltype", "log2FC", "se", "statistic", "pvalue",
            "pvalue_underflow_clamped", "padj", "exposure", "n_donors",
            "n_unexposed", "n_exposed", "min_cells", "median_cells",
            "method"
        )
    )
    expect_identical(result$pvalue_underflow_clamped, c(TRUE, FALSE))
    expect_equal(result$pvalue[[1L]], .Machine$double.xmin)
    expect_equal(
        result$padj,
        stats::p.adjust(result$pvalue, method = "BH")
    )
    expect_true(all(result$exposure == "smoking"))
})

test_that("DESeq2 backend emits the current result contract", {
    skip_if_not_installed("DESeq2")
    result <- suppressWarnings(run_sc_exwas(
        .make_backend_schema_scee(),
        exposure = "smoking",
        celltype_col = "cell_type",
        method = "DESeq2",
        min_cells = 5L,
        min_donors = 6L,
        filter_genes = FALSE
    ))
    expect_gt(nrow(result), 0L)
    expect_true(all(c(
        "log2FC", "se", "statistic", "pvalue", "padj", "padj_global",
        "baseMean", "exposure", "method"
    ) %in% colnames(result)))
    expect_false(any(c(
        "log2FoldChange", "lfcSE", "stat"
    ) %in% colnames(result)))
    expect_true(all(result$method == "DESeq2_Wald"))
})

test_that("voom-dream backend emits the current result contract", {
    skip_if_not_installed("variancePartition")
    skip_if_not_installed("limma")
    result <- suppressWarnings(run_sc_exwas(
        .make_backend_schema_scee(),
        exposure = "smoking",
        celltype_col = "cell_type",
        method = "dreamlet",
        min_cells = 5L,
        min_donors = 6L,
        filter_genes = FALSE
    ))
    expect_gt(nrow(result), 0L)
    expect_true(all(c(
        "log2FC", "se", "statistic", "pvalue", "padj", "padj_global",
        "logCPM", "exposure", "method"
    ) %in% colnames(result)))
    expect_false(any(c(
        "log2FoldChange", "lfcSE", "stat"
    ) %in% colnames(result)))
    expect_true(all(result$method == "voom_dream"))
})

test_that("backend result construction rejects malformed vectors", {
    arguments <- list(
        gene = c("G1", "G2"),
        celltype = "Mono",
        log2FC = c(0.5, -0.2),
        se = c(0.1, 0.2),
        statistic = c(5, -1),
        pvalue = c(0.01, 0.4),
        exposure = "smoking",
        n_donors = 20L,
        method = "test_backend"
    )
    bad_length <- arguments
    bad_length$se <- 0.1
    expect_error(
        do.call(.sc_exwas_result_frame, bad_length),
        "incompatible lengths"
    )

    bad_pvalue <- arguments
    bad_pvalue$pvalue <- c(-0.1, 0.4)
    expect_error(
        do.call(.sc_exwas_result_frame, bad_pvalue),
        "p-values"
    )
})

test_that("non-edgeR backends never ignore test_features", {
    empty <- methods::new("SingleCellExposomeExperiment")
    expect_error(
        run_sc_exwas(
            empty,
            exposure = "smoking",
            celltype_col = "cell_type",
            method = "DESeq2",
            test_features = "G1"
        ),
        "supported only by method='edgeR'"
    )
    expect_error(
        run_sc_exwas(
            empty,
            exposure = "smoking",
            celltype_col = "cell_type",
            method = "dreamlet",
            test_features = "G1"
        ),
        "supported only by method='edgeR'"
    )
})

test_that("IERS accepts current and legacy association-statistic fields", {
    current <- data.frame(
        gene = c("G1", "G2", "G3"),
        statistic = c(3, -1, 2),
        stringsAsFactors = FALSE
    )
    legacy <- data.frame(
        gene = current$gene,
        stat = current$statistic,
        stringsAsFactors = FALSE
    )

    current_result <- compute_iers(current)
    legacy_result <- compute_iers(legacy)
    expect_equal(current_result, legacy_result)
    expect_identical(current_result$gene[[1L]], "G1")
})

test_that("result-column resolver prioritises current fields", {
    data <- data.frame(statistic = 1, stat = 2)
    expect_identical(
        .resolve_exwas_result_column(
            data,
            c("statistic", "stat"),
            "association-statistic"
        ),
        "statistic"
    )
    expect_error(
        .resolve_exwas_result_column(
            data.frame(other = 1),
            c("statistic", "stat"),
            "association-statistic"
        ),
        "No supported association-statistic column"
    )
})

test_that("GSEA keeps exposure-specific hypothesis strata separate", {
    skip_if_not_installed("fgsea")
    genes <- paste0("G", seq_len(30))
    result <- rbind(
        data.frame(
            gene = genes,
            celltype = "Mono",
            exposure = "E1",
            statistic = seq(-3, 3, length.out = length(genes))
        ),
        data.frame(
            gene = genes,
            celltype = "Mono",
            exposure = "E2",
            statistic = rev(seq(-3, 3, length.out = length(genes)))
        )
    )
    gene_sets <- list(
        pathway_a = genes[1:10],
        pathway_b = genes[11:20]
    )

    set.seed(20260712)
    enrichment <- suppressWarnings(run_gsea(
        result,
        gene_sets,
        min_size = 5L,
        max_size = 20L
    ))
    expect_gt(nrow(enrichment), 0L)
    expect_setequal(unique(enrichment$exposure), c("E1", "E2"))
    expect_true(all(enrichment$celltype == "Mono"))
})

test_that("compute_iers validates weights and ranks unrounded scores", {
    exwas <- data.frame(gene = paste0("G", 1:4),
                        statistic = c(1, 1 + 1e-9, 3, -2))
    expect_error(compute_iers(exwas, weights = "adaptive"), "weights")
    expect_error(compute_iers(exwas, weights = c(1, -1, 1)), "weights")
    res <- compute_iers(exwas, weights = c(1, 0, 0))
    expect_false(anyDuplicated(res$IERS_rank) > 0)
    expect_error(compute_iers(exwas, celltype = "Mono"), "celltype")
})
