# R/composition.R
# Donor-level exposure--cell-composition analysis

#' @include AllClasses.R
#' @include utils.R
#' @importFrom stats coef p.adjust
NULL

.composition_count_threshold <- function(value, name, minimum = 0L) {
    if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
            !is.finite(value) ||
            value != as.integer(value) || value < minimum) {
        stop(name, " must be one integer >= ", minimum, ".")
    }
    as.integer(value)
}

.composition_replicate_key <- function(donor, replicate) {
    paste0(
        nchar(donor, type = "bytes"), ":", donor, "|",
        nchar(replicate, type = "bytes"), ":", replicate
    )
}

.composition_count_table <- function(cell_types, donors, replicates,
                                      min_cells,
                                      support_min_donors,
                                      support_min_total_cells,
                                      celltype_universe = NULL,
                                      exclude_celltypes = character()) {
    if (!identical(length(cell_types), length(donors)) ||
            !identical(length(cell_types), length(replicates))) {
        stop("Cell-type, donor, and replicate vectors must have equal length.")
    }
    if (!length(cell_types)) {
        stop("No cells are available for composition analysis.")
    }
    if (anyNA(cell_types) || anyNA(donors) || anyNA(replicates) ||
            any(!nzchar(cell_types)) || any(!nzchar(donors)) ||
            any(!nzchar(replicates))) {
        stop(
            "Cell-type, donor, and replicate labels must be non-missing ",
            "and non-empty."
        )
    }

    observed_types <- sort(unique(cell_types))
    exclude_celltypes <- as.character(exclude_celltypes)
    if (anyNA(exclude_celltypes) || any(!nzchar(exclude_celltypes))) {
        stop("exclude_celltypes must contain non-missing, non-empty labels.")
    }
    if (anyDuplicated(exclude_celltypes)) {
        stop("exclude_celltypes must not contain duplicate labels.")
    }

    if (is.null(celltype_universe)) {
        candidate_types <- setdiff(observed_types, exclude_celltypes)
    } else {
        candidate_types <- unique(as.character(celltype_universe))
        if (!length(candidate_types) || anyNA(candidate_types) ||
                any(!nzchar(candidate_types))) {
            stop(
                "celltype_universe must contain unique, non-missing, ",
                "non-empty labels."
            )
        }
        if (anyDuplicated(celltype_universe)) {
            stop("celltype_universe must not contain duplicate labels.")
        }
        candidate_types <- sort(setdiff(
            candidate_types,
            exclude_celltypes
        ))
    }
    if (length(candidate_types) < 2L) {
        stop("Need at least two candidate cell types for composition analysis.")
    }

    replicate_keys <- .composition_replicate_key(donors, replicates)
    replicate_frame <- unique(data.frame(
        key = replicate_keys,
        donor_id = donors,
        replicate_id = replicates,
        stringsAsFactors = FALSE
    ))
    replicate_frame <- replicate_frame[order(
        replicate_frame$donor_id,
        replicate_frame$replicate_id
    ), , drop = FALSE]
    rownames(replicate_frame) <- NULL

    replicate_index <- match(replicate_keys, replicate_frame$key)
    celltype_index <- match(cell_types, candidate_types)
    counted <- !is.na(celltype_index)
    linear_index <- replicate_index[counted] +
        (celltype_index[counted] - 1L) * nrow(replicate_frame)
    counts <- matrix(
        tabulate(
            linear_index,
            nbins = nrow(replicate_frame) * length(candidate_types)
        ),
        nrow = nrow(replicate_frame),
        ncol = length(candidate_types),
        dimnames = list(replicate_frame$key, candidate_types)
    )

    active_types <- rep(TRUE, ncol(counts))
    eligible_replicates <- rep(TRUE, nrow(counts))
    repeat {
        next_replicates <- rowSums(
            counts[, active_types, drop = FALSE]
        ) >= min_cells
        if (!any(next_replicates)) {
            stop(
                "No biological sample retains at least ", min_cells,
                " cells in the supported composition universe."
            )
        }

        donor_counts <- rowsum(
            counts[next_replicates, , drop = FALSE],
            group = replicate_frame$donor_id[next_replicates],
            reorder = TRUE
        )
        donor_presence <- colSums(donor_counts > 0)
        total_cells <- colSums(counts[next_replicates, , drop = FALSE])
        next_types <- active_types &
            donor_presence >= support_min_donors &
            total_cells >= support_min_total_cells

        if (sum(next_types) < 2L) {
            stop(
                "Fewer than two cell types meet the exposure-independent ",
                "support thresholds."
            )
        }
        if (identical(next_types, active_types) &&
                identical(next_replicates, eligible_replicates)) {
            eligible_replicates <- next_replicates
            break
        }
        active_types <- next_types
        eligible_replicates <- next_replicates
    }

    supported_counts <- counts[
        eligible_replicates,
        active_types,
        drop = FALSE
    ]
    supported_replicates <- replicate_frame[
        eligible_replicates,
        ,
        drop = FALSE
    ]
    supported_donor_counts <- rowsum(
        supported_counts,
        group = supported_replicates$donor_id,
        reorder = TRUE
    )

    final_donor_counts_all <- rowsum(
        counts[eligible_replicates, , drop = FALSE],
        group = supported_replicates$donor_id,
        reorder = TRUE
    )
    final_presence <- colSums(final_donor_counts_all > 0)
    final_total <- colSums(counts[eligible_replicates, , drop = FALSE])
    candidate_support <- data.frame(
        celltype = candidate_types,
        n_donors_present = as.integer(final_presence),
        total_cells = as.integer(final_total),
        included = active_types,
        stringsAsFactors = FALSE
    )
    candidate_support$reason <- "included"
    insufficient_donors <- candidate_support$n_donors_present <
        support_min_donors
    insufficient_cells <- candidate_support$total_cells <
        support_min_total_cells
    candidate_support$reason[!candidate_support$included &
        insufficient_donors & insufficient_cells] <-
        "insufficient_donor_and_cell_support"
    candidate_support$reason[!candidate_support$included &
        insufficient_donors & !insufficient_cells] <-
        "insufficient_donor_support"
    candidate_support$reason[!candidate_support$included &
        !insufficient_donors & insufficient_cells] <-
        "insufficient_total_cell_support"

    reported_types <- sort(unique(c(observed_types, candidate_types)))
    support <- data.frame(
        celltype = reported_types,
        n_donors_present = as.integer(vapply(
            reported_types,
            function(celltype) {
                length(unique(donors[cell_types == celltype]))
            },
            integer(1)
        )),
        total_cells = as.integer(vapply(
            reported_types,
            function(celltype) sum(cell_types == celltype),
            integer(1)
        )),
        included = FALSE,
        reason = "outside_pre_specified_universe",
        stringsAsFactors = FALSE
    )
    support_index <- match(candidate_support$celltype, support$celltype)
    support[support_index, c(
        "n_donors_present",
        "total_cells",
        "included",
        "reason"
    )] <- candidate_support[, c(
        "n_donors_present",
        "total_cells",
        "included",
        "reason"
    )]
    excluded_index <- support$celltype %in% exclude_celltypes
    support$included[excluded_index] <- FALSE
    support$reason[excluded_index] <- "pre_specified_exclusion"
    support$count_basis <- ifelse(
        support$celltype %in% candidate_types,
        "eligible_samples_after_filter_closure",
        "all_input_cells"
    )
    support$support_min_donors <- support_min_donors
    support$support_min_total_cells <- support_min_total_cells

    list(
        counts = supported_counts,
        replicate_data = supported_replicates,
        donor_counts = supported_donor_counts,
        support = support
    )
}

.composition_helmert_basis <- function(n_parts) {
    if (length(n_parts) != 1L || n_parts < 2L ||
            n_parts != as.integer(n_parts)) {
        stop("n_parts must be one integer >= 2.")
    }
    basis <- matrix(0, nrow = n_parts, ncol = n_parts - 1L)
    for (column in seq_len(n_parts - 1L)) {
        basis[seq_len(column), column] <-
            1 / sqrt(column * (column + 1))
        basis[column + 1L, column] <-
            -column / sqrt(column * (column + 1))
    }
    colnames(basis) <- paste0("ilr", seq_len(ncol(basis)))
    basis
}

.composition_logratio_coordinates <- function(counts, prior = 0.5) {
    if (!is.matrix(counts) || !is.numeric(counts) || ncol(counts) < 2L ||
            anyNA(counts) || any(!is.finite(counts)) || any(counts < 0)) {
        stop(
            "counts must be a finite, non-negative numeric matrix with ",
            "at least two columns."
        )
    }
    if (!is.numeric(prior) || length(prior) != 1L ||
            !is.finite(prior) || prior <= 0) {
        stop("prior must be one finite number > 0.")
    }

    expected_log <- digamma(counts + prior)
    clr <- expected_log - rowMeans(expected_log)
    basis <- .composition_helmert_basis(ncol(counts))
    ilr <- clr %*% basis

    remainder_mean <- (
        rowSums(expected_log) - expected_log
    ) / (ncol(counts) - 1L)
    cell_vs_rest_log2 <- (expected_log - remainder_mean) / log(2)
    colnames(ilr) <- colnames(basis)
    colnames(cell_vs_rest_log2) <- colnames(counts)
    rownames(ilr) <- rownames(counts)
    rownames(cell_vs_rest_log2) <- rownames(counts)

    list(
        ilr = ilr,
        cell_vs_rest_log2 = cell_vs_rest_log2,
        basis = basis
    )
}

.composition_aggregate_rows <- function(values, donor_ids) {
    donor_levels <- sort(unique(donor_ids))
    donor_factor <- factor(donor_ids, levels = donor_levels)
    sums <- rowsum(values, donor_factor, reorder = FALSE)
    n_replicates <- tabulate(as.integer(donor_factor), length(donor_levels))
    means <- sweep(sums, 1L, n_replicates, "/")
    rownames(means) <- donor_levels
    list(values = means, n_replicates = stats::setNames(
        n_replicates,
        donor_levels
    ))
}

.composition_donor_coordinates <- function(count_object,
                                             replicate_aggregation,
                                             prior) {
    counts <- count_object$counts
    donor_ids <- count_object$replicate_data$donor_id
    pooled_counts <- rowsum(counts, donor_ids, reorder = TRUE)
    observed_proportion <- pooled_counts / rowSums(pooled_counts)

    if (replicate_aggregation == "equal_replicate") {
        replicate_coordinates <- .composition_logratio_coordinates(
            counts,
            prior = prior
        )
        ilr <- .composition_aggregate_rows(
            replicate_coordinates$ilr,
            donor_ids
        )
        log2_ratio <- .composition_aggregate_rows(
            replicate_coordinates$cell_vs_rest_log2,
            donor_ids
        )
        donor_ilr <- ilr$values
        donor_log2_ratio <- log2_ratio$values
        n_replicates <- ilr$n_replicates
        basis <- replicate_coordinates$basis
    } else {
        donor_coordinates <- .composition_logratio_coordinates(
            pooled_counts,
            prior = prior
        )
        donor_ilr <- donor_coordinates$ilr
        donor_log2_ratio <- donor_coordinates$cell_vs_rest_log2
        n_replicates <- table(factor(
            donor_ids,
            levels = rownames(pooled_counts)
        ))
        n_replicates <- stats::setNames(
            as.integer(n_replicates),
            rownames(pooled_counts)
        )
        basis <- donor_coordinates$basis
    }

    list(
        ilr = donor_ilr,
        cell_vs_rest_log2 = donor_log2_ratio,
        observed_proportion = observed_proportion,
        pooled_counts = pooled_counts,
        n_replicates = n_replicates,
        basis = basis
    )
}

.composition_pillai_statistic <- function(response, reduced_qr, full_qr) {
    reduced_residuals <- qr.resid(reduced_qr, response)
    full_residuals <- qr.resid(full_qr, response)
    error_sscp <- crossprod(full_residuals)
    hypothesis_sscp <- crossprod(reduced_residuals) - error_sscp
    total_sscp <- error_sscp + hypothesis_sscp
    solved <- tryCatch(
        qr.solve(total_sscp, hypothesis_sscp, tol = 1e-10),
        error = function(error) {
            stop(
                "The omnibus residual covariance is singular: ",
                conditionMessage(error)
            )
        }
    )
    statistic <- sum(diag(solved))
    if (!is.finite(statistic)) {
        stop("The omnibus Pillai statistic is not finite.")
    }
    max(0, min(1, statistic))
}

.composition_asymptotic_pillai <- function(response,
                                             reduced_design,
                                             full_design) {
    reduced_fit <- stats::lm(response ~ reduced_design - 1)
    full_fit <- stats::lm(response ~ full_design - 1)
    comparison <- tryCatch(
        stats::anova(reduced_fit, full_fit, test = "Pillai"),
        error = function(error) {
            stop(
                "The partial Pillai test could not be estimated: ",
                conditionMessage(error)
            )
        }
    )
    result <- comparison[2L, , drop = FALSE]
    data.frame(
        statistic = as.numeric(result[["Pillai"]]),
        approximate_f = as.numeric(result[["approx F"]]),
        numerator_df = as.numeric(result[["num Df"]]),
        denominator_df = as.numeric(result[["den Df"]]),
        pvalue_asymptotic = as.numeric(result[["Pr(>F)"]]),
        stringsAsFactors = FALSE
    )
}

.composition_preserve_rng <- function(seed, code) {
    if (!is.numeric(seed) || length(seed) != 1L || is.na(seed) ||
            !is.finite(seed) ||
            seed != as.integer(seed)) {
        stop("seed must be one finite integer.")
    }
    .local_rng_scope(as.integer(seed))
    code()
}

.composition_freedman_lane <- function(response, reduced_design,
                                         full_design, strata,
                                         n_permutations, seed) {
    reduced_qr <- qr(reduced_design)
    full_qr <- qr(full_design)
    observed <- .composition_pillai_statistic(
        response,
        reduced_qr,
        full_qr
    )
    fitted_reduced <- qr.fitted(reduced_qr, response)
    residual_reduced <- qr.resid(reduced_qr, response)
    strata_indices <- split(seq_len(nrow(response)), strata, drop = TRUE)

    permutation_statistics <- .composition_preserve_rng(seed, function() {
        vapply(seq_len(n_permutations), function(iteration) {
            permutation <- seq_len(nrow(response))
            for (indices in strata_indices) {
                permutation[indices] <- sample(
                    indices,
                    length(indices),
                    replace = FALSE
                )
            }
            permuted_response <- fitted_reduced +
                residual_reduced[permutation, , drop = FALSE]
            .composition_pillai_statistic(
                permuted_response,
                reduced_qr,
                full_qr
            )
        }, numeric(1))
    })
    pvalue <- (
        1 + sum(permutation_statistics >= observed - 1e-12)
    ) / (n_permutations + 1)

    list(
        statistic = observed,
        pvalue = pvalue,
        monte_carlo_se = sqrt(
            pvalue * (1 - pvalue) / (n_permutations + 1)
        ),
        permutation_statistics = permutation_statistics
    )
}

.composition_hc3 <- function(response, design, coefficient_index) {
    fit <- stats::lm.fit(design, response)
    if (fit$rank != ncol(design)) {
        stop("The cell-type contrast design is rank deficient.")
    }
    coefficients <- fit$coefficients[coefficient_index, ]
    residuals <- fit$residuals
    if (is.null(dim(residuals))) {
        residuals <- matrix(residuals, ncol = 1L)
    }
    q_matrix <- qr.Q(fit$qr, complete = FALSE)
    leverage <- rowSums(q_matrix^2)
    if (any(1 - leverage <= sqrt(.Machine$double.eps))) {
        stop("HC3 standard errors are undefined because leverage is one.")
    }
    bread <- solve(crossprod(design))
    standard_errors <- vapply(seq_len(ncol(response)), function(column) {
        adjusted_squared <- (
            residuals[, column] / (1 - leverage)
        )^2
        meat <- crossprod(design, design * adjusted_squared)
        covariance <- bread %*% meat %*% bread
        sqrt(max(0, covariance[coefficient_index, coefficient_index]))
    }, numeric(1))
    if (any(!is.finite(standard_errors)) || any(standard_errors <= 0)) {
        stop("HC3 standard errors are zero or non-finite.")
    }
    residual_df <- nrow(design) - ncol(design)
    statistic <- coefficients / standard_errors
    pvalue <- 2 * stats::pt(
        abs(statistic),
        df = residual_df,
        lower.tail = FALSE
    )

    list(
        coefficient = as.numeric(coefficients),
        se = standard_errors,
        statistic = as.numeric(statistic),
        pvalue = as.numeric(pvalue),
        df = residual_df
    )
}

.composition_align_strata <- function(permutation_strata,
                                        exposure_data,
                                        donor_ids) {
    if (is.null(permutation_strata)) {
        return(rep("all_donors", length(donor_ids)))
    }
    if (length(permutation_strata) == 1L &&
            is.character(permutation_strata) &&
            permutation_strata %in% colnames(exposure_data)) {
        strata <- exposure_data[, permutation_strata]
        names(strata) <- rownames(exposure_data)
        return(as.character(strata[donor_ids]))
    }
    if (!is.null(names(permutation_strata))) {
        if (anyDuplicated(names(permutation_strata)) ||
                !all(donor_ids %in% names(permutation_strata))) {
            stop(
                "Named permutation_strata must contain exactly one value ",
                "for every model donor."
            )
        }
        return(as.character(permutation_strata[donor_ids]))
    }
    if (length(permutation_strata) != nrow(exposure_data)) {
        stop(
            "Unnamed permutation_strata must have one value per row of ",
            "exposureData."
        )
    }
    as.character(permutation_strata[
        match(donor_ids, rownames(exposure_data))
    ])
}

#' Test exposure associations with donor-level cell composition
#'
#' Tests whether an exposure is associated with the joint cell-composition
#' profile while retaining donors, rather than cells or repeated biological
#' samples, as the independent units of inference.
#'
#' @param x A \code{SingleCellExposomeExperiment}.
#' @param exposure Character scalar naming the exposure variable.
#' @param celltype_col Character scalar naming the cell-type column in
#'   \code{colData(x)}.
#' @param sample_col Character scalar naming the donor-ID column in
#'   \code{colData(x)}. The historical argument name is retained for API
#'   compatibility.
#' @param replicate_col Optional character scalar naming biological samples
#'   nested within donors. If \code{NULL}, each donor is treated as one sample.
#' @param covariates Character vector naming numeric design columns in
#'   \code{exposureData(x)}.
#' @param min_cells Minimum number of cells from the retained universe per
#'   biological sample.
#' @param support_min_donors Minimum number of donors in which a cell type must
#'   be observed. Support is evaluated without reference to exposure values.
#' @param support_min_total_cells Minimum total number of cells required for a
#'   cell type, also evaluated without reference to exposure values.
#' @param celltype_universe Optional pre-specified character vector of cell
#'   types eligible for analysis.
#' @param exclude_celltypes Character vector of annotation labels excluded
#'   before support filtering.
#' @param zero_prior Positive symmetric Dirichlet prior. The default 0.5 is the
#'   Jeffreys prior and replaces the former arbitrary fixed pseudocount.
#' @param replicate_aggregation Either \code{"equal_replicate"}, the primary
#'   approach that averages log-ratio coordinates across samples within each
#'   donor, or \code{"pooled_counts"}, which pools cells within donor and is
#'   intended as a sensitivity analysis.
#' @param n_permutations Number of Freedman--Lane residual permutations for the
#'   omnibus test. Zero reports the asymptotic partial Pillai test only.
#' @param permutation_strata Optional exchangeability strata for permutation.
#'   Supply a named vector indexed by donor ID for study-stratified inference.
#' @param seed Integer random seed used only for permutation. The caller's RNG
#'   state is restored before return.
#'
#' @return A \code{DataFrame} with one HC3-adjusted cell-type-versus-rest
#'   log-ratio contrast per supported cell type. \code{metadata(result)} stores
#'   the partial Pillai omnibus test, support audit, design diagnostics, and the
#'   donor-level model matrices used for inference. \code{padj} is the Holm
#'   family-wise adjusted p-value; \code{padj_bh} is supplied as a secondary
#'   false-discovery-rate summary.
#'
#' @details
#' Counts are formed for every biological sample by supported cell type. Zeros
#' are handled through posterior expected log proportions under a symmetric
#' Dirichlet(1/2) prior. Orthonormal isometric log-ratio coordinates define a
#' basis-invariant partial Pillai omnibus test. Repeated samples are averaged
#' within donor in log-ratio space, giving every donor equal regression weight.
#' Cell-type contrasts are exposure coefficients for
#' \code{log2(cell type / geometric mean of all other supported types)} with
#' HC3 heteroscedasticity-consistent standard errors.
#' Nominal 95 percent confidence intervals and Bonferroni simultaneous
#' 95 percent family-wise intervals are both returned and labelled explicitly.
#'
#' The method estimates relative composition within the supported, retained
#' cell-type universe. It does not identify absolute cell abundance, and the
#' Dirichlet prior handles observed zeros but does not model dissociation,
#' capture, annotation, or tissue-sampling selection mechanisms.
#'
#' @references
#' Aitchison J (1986). The Statistical Analysis of Compositional Data.
#' Chapman and Hall.
#'
#' Freedman D, Lane D (1983). A nonstochastic interpretation of reported
#' significance levels. Journal of Business and Economic Statistics, 1,
#' 292--298.
#'
#' @export
#' @examples
#' library(SingleCellExperiment)
#' library(S4Vectors)
#' set.seed(1)
#' donors <- paste0("D", seq_len(12))
#' donor_type_counts <- t(rmultinom(
#'     length(donors),
#'     size = 30,
#'     prob = c(A = 0.3, B = 0.4, C = 0.3)
#' ))
#' rownames(donor_type_counts) <- donors
#' cell_donor <- rep(donors, each = 30)
#' cell_type <- unlist(
#'     lapply(seq_along(donors), function(i) {
#'         rep(colnames(donor_type_counts), donor_type_counts[i, ])
#'     }),
#'     use.names = FALSE
#' )
#' counts <- matrix(rpois(20 * length(cell_donor), 5), nrow = 20,
#'     dimnames = list(paste0("G", seq_len(20)),
#'         paste0("c", seq_along(cell_donor))))
#' sce <- SingleCellExperiment(
#'     assays = list(counts = counts),
#'     colData = DataFrame(
#'         cell_id = colnames(counts),
#'         donor_id = cell_donor,
#'         cell_type = cell_type))
#' exp_mat <- cbind(X1 = rep(0:1, 6), age = seq(40, 62, by = 2))
#' rownames(exp_mat) <- donors
#' scee <- build_scee(sce, exp_mat, sample_col = "donor_id")
#' run_exposure_composition(
#'     scee,
#'     exposure = "X1",
#'     covariates = "age",
#'     min_cells = 10L,
#'     support_min_donors = 3L,
#'     support_min_total_cells = 10L
#' )
run_exposure_composition <- function(
        x,
        exposure,
        celltype_col = "cell_type",
        sample_col = "donor_id",
        replicate_col = NULL,
        covariates = NULL,
        min_cells = 10L,
        support_min_donors = 10L,
        support_min_total_cells = 100L,
        celltype_universe = NULL,
        exclude_celltypes = character(),
        zero_prior = 0.5,
        replicate_aggregation = c("equal_replicate", "pooled_counts"),
        n_permutations = 0L,
        permutation_strata = NULL,
        seed = 1L) {
    if (!methods::is(x, "SingleCellExposomeExperiment")) {
        stop("x must be a SingleCellExposomeExperiment.")
    }
    if (!is.character(exposure) || length(exposure) != 1L ||
            is.na(exposure) || !nzchar(exposure)) {
        stop("exposure must be one non-missing, non-empty name.")
    }
    name_arguments <- list(
        celltype_col = celltype_col,
        sample_col = sample_col
    )
    if (!is.null(replicate_col)) {
        name_arguments$replicate_col <- replicate_col
    }
    invalid_name_arguments <- names(name_arguments)[vapply(
        name_arguments,
        function(value) {
            !is.character(value) || length(value) != 1L ||
                is.na(value) || !nzchar(value)
        },
        logical(1)
    )]
    if (length(invalid_name_arguments)) {
        stop(
            paste(invalid_name_arguments, collapse = ", "),
            " must each be one non-missing, non-empty column name."
        )
    }
    replicate_aggregation <- match.arg(replicate_aggregation)
    min_cells <- .composition_count_threshold(min_cells, "min_cells", 1L)
    support_min_donors <- .composition_count_threshold(
        support_min_donors,
        "support_min_donors",
        1L
    )
    support_min_total_cells <- .composition_count_threshold(
        support_min_total_cells,
        "support_min_total_cells",
        1L
    )
    n_permutations <- .composition_count_threshold(
        n_permutations,
        "n_permutations",
        0L
    )

    exposure_data <- methods::slot(x, "exposureData")
    if (!exposure %in% colnames(exposure_data)) {
        stop("Exposure '", exposure, "' not found in exposureData.")
    }
    if (!is.null(covariates) && !is.character(covariates)) {
        stop("covariates must be a character vector or NULL.")
    }
    covariates <- as.character(covariates)
    if (anyNA(covariates) || any(!nzchar(covariates))) {
        stop("covariates must contain non-missing, non-empty names.")
    }
    if (anyDuplicated(covariates)) {
        stop("covariates must not contain duplicate names.")
    }
    if (exposure %in% covariates) {
        stop("The exposure must not also be listed as a covariate.")
    }
    missing_covariates <- setdiff(covariates, colnames(exposure_data))
    if (length(missing_covariates)) {
        stop(
            "Covariate(s) not found in exposureData: ",
            paste(missing_covariates, collapse = ", "),
            "."
        )
    }

    col_data <- SummarizedExperiment::colData(x)
    required_columns <- c(celltype_col, sample_col)
    if (!is.null(replicate_col)) {
        required_columns <- c(required_columns, replicate_col)
    }
    missing_columns <- setdiff(required_columns, colnames(col_data))
    if (length(missing_columns)) {
        stop(
            "Column(s) not found in colData(x): ",
            paste(missing_columns, collapse = ", "),
            "."
        )
    }

    cell_types <- as.character(col_data[[celltype_col]])
    donor_ids <- as.character(col_data[[sample_col]])
    replicate_ids <- if (is.null(replicate_col)) {
        donor_ids
    } else {
        as.character(col_data[[replicate_col]])
    }
    if (!all(unique(donor_ids) %in% rownames(exposure_data))) {
        stop("Some cell-level donor IDs are absent from exposureData.")
    }

    count_object <- .composition_count_table(
        cell_types = cell_types,
        donors = donor_ids,
        replicates = replicate_ids,
        min_cells = min_cells,
        support_min_donors = support_min_donors,
        support_min_total_cells = support_min_total_cells,
        celltype_universe = celltype_universe,
        exclude_celltypes = exclude_celltypes
    )
    if (ncol(count_object$counts) == 2L) {
        warning(
            "Only two supported cell types remain; the two reported ",
            "cell-type contrasts are exact sign reversals.",
            call. = FALSE
        )
    }
    donor_coordinates <- .composition_donor_coordinates(
        count_object,
        replicate_aggregation = replicate_aggregation,
        prior = zero_prior
    )

    coordinate_donors <- rownames(donor_coordinates$ilr)
    design_values <- exposure_data[
        coordinate_donors,
        c(exposure, covariates),
        drop = FALSE
    ]
    complete <- stats::complete.cases(design_values)
    model_donors <- coordinate_donors[complete]
    if (!length(model_donors)) {
        stop("No donor has complete exposure and covariate data.")
    }
    design_values <- design_values[complete, , drop = FALSE]
    if (length(unique(design_values[, exposure])) < 2L) {
        stop("The exposure has fewer than two observed values in model donors.")
    }

    reduced_design <- cbind(
        intercept = 1,
        as.matrix(design_values[, covariates, drop = FALSE])
    )
    full_design <- cbind(
        reduced_design,
        exposure = as.numeric(design_values[, exposure])
    )
    storage.mode(reduced_design) <- "double"
    storage.mode(full_design) <- "double"
    if (any(!is.finite(reduced_design)) || any(!is.finite(full_design))) {
        stop("The complete-case design contains non-finite values.")
    }
    reduced_rank <- qr(reduced_design)$rank
    full_rank <- qr(full_design)$rank
    if (reduced_rank != ncol(reduced_design)) {
        stop("The covariate-only composition design is rank deficient.")
    }
    if (full_rank != ncol(full_design) ||
            full_rank != reduced_rank + 1L) {
        stop(
            "The exposure is aliased with the intercept or covariates in ",
            "the composition design."
        )
    }

    donor_ilr <- donor_coordinates$ilr[model_donors, , drop = FALSE]
    donor_log2_ratio <- donor_coordinates$cell_vs_rest_log2[
        model_donors,
        ,
        drop = FALSE
    ]
    residual_df <- nrow(full_design) - full_rank
    if (residual_df < ncol(donor_ilr)) {
        stop(
            "The omnibus model requires at least as many residual degrees ",
            "of freedom as ILR coordinates."
        )
    }

    asymptotic_omnibus <- .composition_asymptotic_pillai(
        donor_ilr,
        reduced_design,
        full_design
    )
    direct_pillai <- .composition_pillai_statistic(
        donor_ilr,
        qr(reduced_design),
        qr(full_design)
    )
    if (!isTRUE(all.equal(
            direct_pillai,
            asymptotic_omnibus$statistic,
            tolerance = 1e-8
        ))) {
        stop("Direct and asymptotic partial Pillai statistics disagree.")
    }
    strata <- .composition_align_strata(
        permutation_strata,
        exposure_data,
        model_donors
    )
    if (anyNA(strata) || any(!nzchar(strata))) {
        stop("Permutation strata must be non-missing and non-empty.")
    }
    if (n_permutations > 0L && all(table(strata) < 2L)) {
        stop("No permutation stratum contains two or more model donors.")
    }

    permutation <- NULL
    if (n_permutations > 0L) {
        permutation <- .composition_freedman_lane(
            donor_ilr,
            reduced_design,
            full_design,
            strata = factor(strata),
            n_permutations = n_permutations,
            seed = seed
        )
    }
    asymptotic_underflow <- asymptotic_omnibus$pvalue_asymptotic == 0
    reported_asymptotic_pvalue <- asymptotic_omnibus$pvalue_asymptotic
    reported_asymptotic_pvalue[asymptotic_underflow] <-
        .Machine$double.xmin
    omnibus_pvalue <- if (is.null(permutation)) {
        reported_asymptotic_pvalue
    } else {
        permutation$pvalue
    }
    omnibus <- data.frame(
        test = "partial_Pillai_ILR",
        statistic = asymptotic_omnibus$statistic,
        approximate_f = asymptotic_omnibus$approximate_f,
        numerator_df = asymptotic_omnibus$numerator_df,
        denominator_df = asymptotic_omnibus$denominator_df,
        pvalue = omnibus_pvalue,
        pvalue_method = if (is.null(permutation)) {
            "asymptotic_partial_Pillai"
        } else {
            "stratified_Freedman_Lane"
        },
        pvalue_asymptotic = reported_asymptotic_pvalue,
        pvalue_asymptotic_underflow_clamped = asymptotic_underflow,
        n_permutations = n_permutations,
        permutation_pvalue = if (is.null(permutation)) {
            NA_real_
        } else {
            permutation$pvalue
        },
        permutation_monte_carlo_se = if (is.null(permutation)) {
            NA_real_
        } else {
            permutation$monte_carlo_se
        },
        stringsAsFactors = FALSE
    )

    contrast_fit <- .composition_hc3(
        donor_log2_ratio,
        full_design,
        coefficient_index = ncol(full_design)
    )
    reported_pvalue <- contrast_fit$pvalue
    pvalue_underflow <- reported_pvalue == 0
    reported_pvalue[pvalue_underflow] <- .Machine$double.xmin
    critical_value <- stats::qt(0.975, df = contrast_fit$df)
    confidence_low <- contrast_fit$coefficient -
        critical_value * contrast_fit$se
    confidence_high <- contrast_fit$coefficient +
        critical_value * contrast_fit$se
    simultaneous_critical_value <- stats::qt(
        1 - 0.05 / (2 * ncol(donor_log2_ratio)),
        df = contrast_fit$df
    )
    simultaneous_confidence_low <- contrast_fit$coefficient -
        simultaneous_critical_value * contrast_fit$se
    simultaneous_confidence_high <- contrast_fit$coefficient +
        simultaneous_critical_value * contrast_fit$se

    model_replicates <- count_object$replicate_data$donor_id %in%
        model_donors
    model_counts <- count_object$counts[model_replicates, , drop = FALSE]
    model_pooled_counts <- donor_coordinates$pooled_counts[
        model_donors,
        ,
        drop = FALSE
    ]
    model_observed_proportion <- donor_coordinates$observed_proportion[
        model_donors,
        ,
        drop = FALSE
    ]
    celltypes <- colnames(donor_log2_ratio)
    support_index <- match(celltypes, count_object$support$celltype)

    output <- S4Vectors::DataFrame(data.frame(
        celltype = celltypes,
        coefficient = contrast_fit$coefficient,
        log2_ratio_change = contrast_fit$coefficient,
        se = contrast_fit$se,
        statistic = contrast_fit$statistic,
        df = rep(contrast_fit$df, length(celltypes)),
        pvalue = reported_pvalue,
        pvalue_underflow_clamped = pvalue_underflow,
        padj = stats::p.adjust(reported_pvalue, method = "holm"),
        padj_holm = stats::p.adjust(reported_pvalue, method = "holm"),
        padj_bh = stats::p.adjust(reported_pvalue, method = "BH"),
        ci_low = confidence_low,
        ci_high = confidence_high,
        ci_level = rep(0.95, length(celltypes)),
        ci_multiplicity = rep("nominal_unadjusted", length(celltypes)),
        simultaneous_ci_low = simultaneous_confidence_low,
        simultaneous_ci_high = simultaneous_confidence_high,
        simultaneous_ci_method = rep(
            "Bonferroni_95_percent_familywise",
            length(celltypes)
        ),
        ratio_change = 2^contrast_fit$coefficient,
        ratio_ci_low = 2^confidence_low,
        ratio_ci_high = 2^confidence_high,
        ratio_simultaneous_ci_low = 2^simultaneous_confidence_low,
        ratio_simultaneous_ci_high = 2^simultaneous_confidence_high,
        exposure = rep(exposure, length(celltypes)),
        n_donors = rep(length(model_donors), length(celltypes)),
        n_replicates = rep(sum(model_replicates), length(celltypes)),
        n_donors_present = colSums(model_pooled_counts > 0),
        zero_replicates = colSums(model_counts == 0),
        total_cells = colSums(model_counts),
        median_observed_proportion = apply(
            model_observed_proportion,
            2L,
            stats::median
        ),
        support_n_donors = count_object$support$n_donors_present[
            support_index
        ],
        support_total_cells = count_object$support$total_cells[
            support_index
        ],
        method = rep(
            "donor_logratio_HC3_with_partial_Pillai",
            length(celltypes)
        ),
        stringsAsFactors = FALSE
    ))

    S4Vectors::metadata(output) <- list(
        omnibus = S4Vectors::DataFrame(omnibus),
        support = S4Vectors::DataFrame(count_object$support),
        design = list(
            exposure = exposure,
            covariates = covariates,
            n_donors_with_eligible_samples = length(coordinate_donors),
            n_donors = length(model_donors),
            n_donors_dropped_for_incomplete_design = sum(!complete),
            n_replicates = sum(model_replicates),
            reduced_rank = reduced_rank,
            full_rank = full_rank,
            residual_df = residual_df,
            permutation_strata_sizes = sort(table(strata), decreasing = TRUE)
        ),
        transform = list(
            zero_handling = "Dirichlet posterior expected log-ratios",
            zero_prior = zero_prior,
            ilr_basis = donor_coordinates$basis,
            replicate_aggregation = replicate_aggregation,
            estimand = paste(
                "Exposure coefficient for log2(cell type / geometric mean",
                "of other supported cell types)"
            )
        ),
        model_data = list(
            donor_id = model_donors,
            design = full_design,
            donor_ilr = donor_ilr,
            donor_cell_vs_rest_log2 = donor_log2_ratio,
            donor_observed_proportion = model_observed_proportion,
            donor_pooled_counts = model_pooled_counts,
            n_replicates = donor_coordinates$n_replicates[model_donors]
        ),
        permutation_statistics = if (is.null(permutation)) {
            numeric()
        } else {
            permutation$permutation_statistics
        }
    )
    output
}
