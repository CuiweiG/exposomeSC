.make_interventional_scee <- function() {
    set.seed(20260712)
    donors <- paste0("D", seq_len(12))
    exposure <- rep(c(0, 1), each = 6)
    cell_donor <- character()
    celltype <- character()
    for (i in seq_along(donors)) {
        n_type <- c(
            target = 5L + (i %% 3L),
            context_a = 5L + ((i + 1L) %% 4L),
            context_b = 6L + ((i + 2L) %% 3L)
        )
        cell_donor <- c(cell_donor, rep(donors[[i]], sum(n_type)))
        celltype <- c(celltype, rep(names(n_type), n_type))
    }
    genes <- paste0("G", seq_len(12))
    counts <- matrix(
        stats::rnbinom(length(genes) * length(cell_donor), mu = 10, size = 8),
        nrow = length(genes),
        dimnames = list(genes, paste0("cell", seq_along(cell_donor)))
    )
    signal <- celltype == "target" & exposure[match(cell_donor, donors)] == 1
    counts["G1", signal] <- counts["G1", signal] * 2L
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = counts),
        colData = S4Vectors::DataFrame(
            donor_id = cell_donor,
            cell_type = celltype
        )
    )
    exposure_data <- cbind(
        smoking = exposure,
        age_scaled = as.numeric(scale(seq_along(donors)))
    )
    rownames(exposure_data) <- donors
    build_scee(sce, exposure_data, sample_col = "donor_id")
}

.run_small_erd <- function(scee, target_genes) {
    suppressWarnings(run_erd_interventional(
        scee,
        exposure = "smoking",
        celltype = "target",
        celltype_col = "cell_type",
        covariates = "age_scaled",
        target_genes = target_genes,
        n_mc = 10L,
        n_boot = 20L,
        min_cells = 5L,
        min_donors = 10L,
        min_group_donors = 5L,
        min_boot_valid = 0.5,
        adjust = "none",
        seed = 20260712L,
        BPPARAM = BiocParallel::SerialParam(RNGseed = 20260712L)
    ))
}

test_that("ERD offset is invariant to the target-gene set", {
    scee <- .make_interventional_scee()
    one_gene <- .run_small_erd(scee, "G1")
    several_genes <- .run_small_erd(scee, paste0("G", 1:4))
    matched <- several_genes[several_genes$gene == "G1", ]

    expect_equal(one_gene$IDE_RR, matched$IDE_RR, tolerance = 1e-10)
    expect_equal(one_gene$IIE_RR, matched$IIE_RR, tolerance = 1e-10)
    expect_equal(one_gene$OE_RR, matched$OE_RR, tolerance = 1e-10)
})

test_that("ERD finite-resampling p-values cannot be zero", {
    scee <- .make_interventional_scee()
    result <- .run_small_erd(scee, paste0("G", 1:4))
    finite_p <- c(result$p_IDE, result$p_IIE)
    finite_p <- finite_p[is.finite(finite_p)]

    expect_true(all(finite_p > 0))
    expect_true(all(finite_p <= 1))
    expect_true(all(is.na(result$padj_IDE)))
    expect_true(all(is.na(result$padj_IIE)))
    expect_equal(attr(result, "multiplicity")$scope, "none")
})

test_that("CI E-values equal one when the interval includes the null", {
    expect_equal(
        mediational_evalue_ci(
            rr = c(0.5, 2, 0.5, 2),
            lower = c(0.2, 0.5, 0.2, 1.5),
            upper = c(1.2, 3, 0.8, 4)
        ),
        c(1, 1, mediational_evalue(0.8), mediational_evalue(1.5))
    )
})

test_that("ERD rejects cell types without enough eligible donors", {
    scee <- .make_interventional_scee()
    expect_error(
        run_erd_interventional(
            scee,
            exposure = "smoking",
            celltype = "target",
            celltype_col = "cell_type",
            target_genes = "G1",
            min_cells = 100L,
            min_donors = 10L,
            min_group_donors = 5L,
            n_mc = 2L,
            n_boot = 2L,
            BPPARAM = BiocParallel::SerialParam()
        ),
        "only 0 donor"
    )
})
