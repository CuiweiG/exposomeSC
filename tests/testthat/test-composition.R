library(SingleCellExperiment)
library(S4Vectors)

.make_composition_scee <- function(n_donors = 30L, shifted = FALSE) {
    donor_ids <- paste0("D", seq_len(n_donors))
    exposure <- rep(0:1, length.out = n_donors)
    study <- rep(seq_len(3L), length.out = n_donors)

    cell_rows <- vector("list", n_donors * 2L)
    position <- 0L
    for (donor_index in seq_len(n_donors)) {
        for (replicate_index in seq_len(2L)) {
            position <- position + 1L
            cell_counts <- if (shifted && exposure[[donor_index]] == 1) {
                c(A = 55L, B = 20L, C = 15L)
            } else {
                c(A = 30L, B = 30L, C = 30L)
            }
            delta_a <- (donor_index %% 7L) - 3L
            delta_b <- ((3L * donor_index) %% 7L) - 3L
            cell_counts <- cell_counts + c(
                A = delta_a,
                B = delta_b,
                C = -delta_a - delta_b
            )
            if (replicate_index == 2L) {
                cell_counts <- cell_counts + c(A = 3L, B = -1L, C = -2L)
            }
            cell_type <- rep(names(cell_counts), cell_counts)
            if (donor_index <= 3L && replicate_index == 1L) {
                cell_type <- c(cell_type, "RareLabel")
            }
            cell_rows[[position]] <- data.frame(
                donor_id = donor_ids[[donor_index]],
                sample_id = paste0(
                    donor_ids[[donor_index]],
                    "_S",
                    replicate_index
                ),
                cell_type = cell_type,
                stringsAsFactors = FALSE
            )
        }
    }
    cell_data <- do.call(rbind, cell_rows)
    cell_ids <- paste0("c", seq_len(nrow(cell_data)))
    counts <- matrix(
        1L,
        nrow = 5L,
        ncol = nrow(cell_data),
        dimnames = list(paste0("G", seq_len(5L)), cell_ids)
    )
    sce <- SingleCellExperiment(
        assays = list(counts = counts),
        colData = DataFrame(
            cell_id = cell_ids,
            donor_id = cell_data$donor_id,
            sample_id = cell_data$sample_id,
            cell_type = cell_data$cell_type
        )
    )
    exposure_matrix <- cbind(
        X = exposure,
        age_scaled = as.numeric(scale(40 + (seq_len(n_donors) %% 11L))),
        study_2 = as.numeric(study == 2L),
        study_3 = as.numeric(study == 3L),
        alias = exposure
    )
    rownames(exposure_matrix) <- donor_ids
    object <- build_scee(
        sce,
        exposure_matrix,
        sample_col = "donor_id"
    )
    list(
        object = object,
        study = stats::setNames(paste0("study", study), donor_ids)
    )
}

.run_synthetic_composition <- function(
        object,
        covariates = c("age_scaled", "study_2", "study_3"),
        ...) {
    run_exposure_composition(
        object,
        exposure = "X",
        celltype_col = "cell_type",
        sample_col = "donor_id",
        replicate_col = "sample_id",
        covariates = covariates,
        min_cells = 20L,
        support_min_donors = 5L,
        support_min_total_cells = 50L,
        ...
    )
}

test_that("composition analysis returns donor-level contrasts and omnibus", {
    synthetic <- .make_composition_scee()
    result <- .run_synthetic_composition(synthetic$object)
    result_metadata <- S4Vectors::metadata(result)

    expect_s4_class(result, "DataFrame")
    expect_equal(result$celltype, c("A", "B", "C"))
    expect_true(all(c(
        "log2_ratio_change",
        "se",
        "pvalue",
        "padj_holm",
        "padj_bh",
        "simultaneous_ci_low",
        "simultaneous_ci_high",
        "ratio_change",
        "n_donors",
        "n_replicates"
    ) %in% colnames(result)))
    expect_equal(result$n_donors, rep(30L, 3L))
    expect_equal(result$n_replicates, rep(60L, 3L))
    expect_equal(
        result$padj,
        stats::p.adjust(result$pvalue, method = "holm")
    )
    expect_s4_class(result_metadata$omnibus, "DataFrame")
    expect_equal(result_metadata$omnibus$test, "partial_Pillai_ILR")
    expect_equal(result_metadata$design$full_rank, 5L)
    expect_equal(
        result_metadata$transform$replicate_aggregation,
        "equal_replicate"
    )
})

test_that("Jeffreys expected log-ratios handle observed zeros", {
    counts <- matrix(
        c(0, 10, 30, 5, 0, 20, 0, 0, 8),
        nrow = 3L,
        byrow = TRUE,
        dimnames = list(NULL, c("A", "B", "C"))
    )
    coordinates <- .composition_logratio_coordinates(counts, prior = 0.5)

    expect_true(all(is.finite(coordinates$ilr)))
    expect_true(all(is.finite(coordinates$cell_vs_rest_log2)))
    expect_equal(
        unname(crossprod(coordinates$basis)),
        diag(2L),
        tolerance = 1e-12
    )
    expect_equal(
        unname(colSums(coordinates$basis)),
        c(0, 0),
        tolerance = 1e-12
    )
})

test_that("equal-replicate aggregation gives samples equal donor weight", {
    sample_counts <- matrix(
        c(
            900, 50, 50,
            10, 40, 50,
            30, 40, 30,
            50, 20, 30
        ),
        nrow = 4L,
        byrow = TRUE,
        dimnames = list(paste0("sample", seq_len(4L)), c("A", "B", "C"))
    )
    count_object <- list(
        counts = sample_counts,
        replicate_data = data.frame(
            donor_id = c("D1", "D1", "D2", "D2"),
            stringsAsFactors = FALSE
        )
    )
    equal_replicate <- .composition_donor_coordinates(
        count_object,
        replicate_aggregation = "equal_replicate",
        prior = 0.5
    )
    pooled <- .composition_donor_coordinates(
        count_object,
        replicate_aggregation = "pooled_counts",
        prior = 0.5
    )
    sample_coordinates <- .composition_logratio_coordinates(
        sample_counts,
        prior = 0.5
    )

    expect_equal(
        equal_replicate$ilr["D1", ],
        colMeans(sample_coordinates$ilr[1:2, , drop = FALSE])
    )
    expect_false(isTRUE(all.equal(
        equal_replicate$ilr["D1", ],
        pooled$ilr["D1", ]
    )))
    expect_equal(equal_replicate$n_replicates, c(D1 = 2L, D2 = 2L))
})

test_that("the supported universe is exposure independent", {
    synthetic <- .make_composition_scee()
    result_a <- .run_synthetic_composition(synthetic$object)

    permuted <- synthetic$object
    exposure_data <- exposureData(permuted)
    exposure_data[, "X"] <- rev(exposure_data[, "X"])
    methods::slot(permuted, "exposureData") <- exposure_data
    methods::validObject(permuted)
    result_b <- .run_synthetic_composition(permuted)

    support_a <- as.data.frame(S4Vectors::metadata(result_a)$support)
    support_b <- as.data.frame(S4Vectors::metadata(result_b)$support)
    expect_identical(support_a, support_b)
    expect_equal(result_a$celltype, result_b$celltype)
})

test_that("composition estimates are invariant to cell column order", {
    synthetic <- .make_composition_scee()
    original <- .run_synthetic_composition(synthetic$object)
    reordered_object <- synthetic$object[
        ,
        rev(seq_len(ncol(synthetic$object)))
    ]
    reordered <- .run_synthetic_composition(reordered_object)

    expect_equal(original$coefficient, reordered$coefficient, tolerance = 0)
    expect_equal(original$pvalue, reordered$pvalue, tolerance = 0)
    expect_equal(
        S4Vectors::metadata(original)$omnibus$statistic,
        S4Vectors::metadata(reordered)$omnibus$statistic,
        tolerance = 0
    )
})

test_that("rare and excluded labels are audited rather than silently used", {
    synthetic <- .make_composition_scee()
    result <- .run_synthetic_composition(
        synthetic$object,
        exclude_celltypes = "RareLabel"
    )
    support <- as.data.frame(S4Vectors::metadata(result)$support)
    rare <- support[support$celltype == "RareLabel", , drop = FALSE]

    expect_false(rare$included)
    expect_equal(rare$reason, "pre_specified_exclusion")
    expect_false("RareLabel" %in% result$celltype)
})

test_that("synthetic shifts are detected by reproducible permutation", {
    synthetic <- .make_composition_scee(n_donors = 36L, shifted = TRUE)
    set.seed(2026)
    rng_before <- .Random.seed
    result_a <- .run_synthetic_composition(
        synthetic$object,
        n_permutations = 49L,
        permutation_strata = synthetic$study,
        seed = 17L
    )
    rng_after <- .Random.seed
    result_b <- .run_synthetic_composition(
        synthetic$object,
        n_permutations = 49L,
        permutation_strata = synthetic$study,
        seed = 17L
    )
    omnibus_a <- S4Vectors::metadata(result_a)$omnibus
    omnibus_b <- S4Vectors::metadata(result_b)$omnibus

    expect_identical(rng_before, rng_after)
    expect_lte(omnibus_a$permutation_pvalue, 0.06)
    expect_identical(
        S4Vectors::metadata(result_a)$permutation_statistics,
        S4Vectors::metadata(result_b)$permutation_statistics
    )
    expect_gt(result_a$coefficient[result_a$celltype == "A"], 0)
    expect_lt(result_a$coefficient[result_a$celltype == "B"], 0)
    expect_lt(result_a$coefficient[result_a$celltype == "C"], 0)
})

test_that("pooled-count sensitivity retains donors and universe", {
    synthetic <- .make_composition_scee(shifted = TRUE)
    equal_replicate <- .run_synthetic_composition(synthetic$object)
    supported <- equal_replicate$celltype
    pooled <- .run_synthetic_composition(
        synthetic$object,
        celltype_universe = supported,
        replicate_aggregation = "pooled_counts"
    )

    expect_equal(pooled$celltype, supported)
    expect_equal(pooled$n_donors, equal_replicate$n_donors)
    expect_true(all(is.finite(pooled$coefficient)))
    expect_equal(
        S4Vectors::metadata(pooled)$transform$replicate_aggregation,
        "pooled_counts"
    )
})

test_that("composition inputs and rank are checked strictly", {
    synthetic <- .make_composition_scee()
    expect_error(
        run_exposure_composition(synthetic$object, exposure = "missing"),
        "not found"
    )
    expect_error(
        .run_synthetic_composition(
            synthetic$object,
            covariates = "alias"
        ),
        "aliased"
    )
    expect_error(
        .run_synthetic_composition(
            synthetic$object,
            permutation_strata = c("only", "two"),
            n_permutations = 2L
        ),
        "one value per row"
    )
})

test_that("the omnibus refuses insufficient residual dimensions", {
    synthetic <- .make_composition_scee(n_donors = 3L)
    expect_error(
        run_exposure_composition(
            synthetic$object,
            exposure = "X",
            celltype_col = "cell_type",
            sample_col = "donor_id",
            replicate_col = "sample_id",
            min_cells = 20L,
            support_min_donors = 2L,
            support_min_total_cells = 10L
        ),
        "residual degrees of freedom"
    )
})
