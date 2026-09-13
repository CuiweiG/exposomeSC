# R/null-calibration.R
# Internal infrastructure for study-stratified empirical null diagnostics.

# The functions in this file deliberately remain internal. The permutation
# diagnostic is an analysis-quality-control procedure, not a general-purpose
# observational conditional-randomisation test.

.null_validate_positive_integer <- function(value, name, minimum = 1L) {
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
            value < minimum || value != as.integer(value)) {
        stop(name, " must be one integer >= ", minimum, ".")
    }
    as.integer(value)
}

.null_replicate_seeds <- function(master_seed, n_replicates) {
    master_seed <- .null_validate_positive_integer(
        master_seed,
        "master_seed"
    )
    n_replicates <- .null_validate_positive_integer(
        n_replicates,
        "n_replicates"
    )
    .local_rng_scope(master_seed, kind = "L'Ecuyer-CMRG")
    sample.int(.Machine$integer.max, n_replicates, replace = FALSE)
}

.null_stratified_permute <- function(exposure, study, seed = NULL) {
    if (!is.numeric(exposure) && !is.logical(exposure)) {
        stop("exposure must be numeric or logical.")
    }
    if (length(exposure) != length(study) || !length(exposure)) {
        stop("exposure and study must have the same positive length.")
    }
    if (anyNA(exposure) || anyNA(study) || any(!is.finite(exposure))) {
        stop("exposure and study must be complete and finite.")
    }
    if (length(unique(exposure)) != 2L) {
        stop("exposure must contain exactly two observed levels.")
    }
    if (!is.null(seed)) {
        seed <- .null_validate_positive_integer(seed, "seed")
        .local_rng_scope(seed, kind = "L'Ecuyer-CMRG")
    }

    study <- as.character(study)
    permuted <- exposure
    for (index in split(seq_along(exposure), study, drop = TRUE)) {
        values <- exposure[index]
        if (length(unique(values)) == 2L) {
            permuted[index] <- sample(values, length(values), replace = FALSE)
        }
        ## Monomorphic studies are deliberately unchanged. They provide no
        ## within-study exposure-label permutation information.
    }
    names(permuted) <- names(exposure)
    permuted
}

.null_within_study_variation <- function(exposure, study) {
    if (length(exposure) != length(study)) {
        stop("exposure and study must have equal length.")
    }
    any(vapply(
        split(exposure, as.character(study), drop = TRUE),
        function(value) length(unique(value)) == 2L,
        logical(1)
    ))
}

.null_draw_eligible_permutation <- function(
    exposure,
    study,
    caches,
    seed,
    min_group_donors = 10L,
    max_attempts = 10000L
) {
    seed <- .null_validate_positive_integer(seed, "seed")
    min_group_donors <- .null_validate_positive_integer(
        min_group_donors,
        "min_group_donors",
        minimum = 2L
    )
    max_attempts <- .null_validate_positive_integer(
        max_attempts,
        "max_attempts"
    )
    if (!is.list(caches) || !length(caches)) {
        stop("caches must be a non-empty list of fixed cell-type caches.")
    }
    if (is.null(names(exposure)) || anyDuplicated(names(exposure)) ||
            any(vapply(
                caches,
                function(cache) {
                    is.null(cache$donors) ||
                        any(!cache$donors %in% names(exposure)) ||
                        is.null(cache$covariates) ||
                        is.null(cache$study)
                },
                logical(1)
            ))) {
        stop("Every cache must map exactly to named exposure donors.")
    }

    .local_rng_scope(seed, kind = "L'Ecuyer-CMRG")
    for (attempt in seq_len(max_attempts)) {
        candidate <- .null_stratified_permute(exposure, study)
        eligible <- vapply(caches, function(cache) {
            design <- .null_build_design(
                exposure = as.numeric(candidate[cache$donors]),
                covariates = cache$covariates,
                study = cache$study,
                min_group_donors = min_group_donors
            )
            isTRUE(design$ok)
        }, logical(1))
        if (all(eligible)) {
            return(list(
                exposure = candidate,
                attempts = as.integer(attempt),
                acceptance_condition = paste(
                    "all fixed cell-type donor subsets satisfy the original",
                    "group-support and design-identifiability criteria"
                )
            ))
        }
    }
    stop(
        "No analysis-eligible study-stratified permutation was found after ",
        max_attempts, " deterministic attempts."
    )
}

.null_fixed_support_cache <- function(counts, host_features,
                                      min_cpm = 1,
                                      min_donors = 10L,
                                      min_total_count = 15L) {
    if (!requireNamespace("edgeR", quietly = TRUE)) {
        stop("Package 'edgeR' is required.")
    }
    counts <- as.matrix(counts)
    if (is.null(rownames(counts)) || is.null(colnames(counts)) ||
            anyDuplicated(rownames(counts)) || anyDuplicated(colnames(counts))) {
        stop("counts must have unique gene and donor dimnames.")
    }
    if (any(!is.finite(counts)) || any(counts < 0)) {
        stop("counts must be finite and non-negative.")
    }
    if (!is.character(host_features) || anyNA(host_features) ||
            anyDuplicated(host_features)) {
        stop("host_features must be a unique, complete character vector.")
    }
    if (!is.numeric(min_cpm) || length(min_cpm) != 1L ||
            !is.finite(min_cpm) || min_cpm <= 0) {
        stop("min_cpm must be one finite number > 0.")
    }
    min_donors <- .null_validate_positive_integer(
        min_donors,
        "min_donors"
    )
    if (!is.numeric(min_total_count) || length(min_total_count) != 1L ||
            !is.finite(min_total_count) || min_total_count < 0) {
        stop("min_total_count must be one finite non-negative number.")
    }

    full_library_size <- colSums(counts)
    if (any(full_library_size <= 0)) {
        stop("Every donor pseudobulk must have a positive library size.")
    }
    dge <- edgeR::DGEList(
        counts = counts,
        lib.size = full_library_size
    )
    donor_cpm <- edgeR::cpm(
        dge,
        normalized.lib.sizes = FALSE,
        log = FALSE
    )
    donor_support <- rowSums(donor_cpm >= min_cpm)
    total_count <- rowSums(counts)
    keep <- donor_support >= min_donors & total_count >= min_total_count
    if (!any(keep)) {
        stop("No feature passed the fixed-support filter.")
    }

    dge <- dge[keep, , keep.lib.sizes = TRUE]
    dge <- edgeR::normLibSizes(dge, method = "TMM")
    retained_features <- rownames(dge$counts)
    retained_host <- retained_features[retained_features %in% host_features]
    if (!length(retained_host)) {
        stop("No host-expression feature passed the fixed-support filter.")
    }

    list(
        counts = dge$counts,
        full_library_size = full_library_size,
        normalisation_factors = dge$samples$norm.factors,
        retained_features = retained_features,
        host_feature_ids = retained_host,
        donor_support = donor_support[keep],
        total_count = total_count[keep],
        n_features_input = nrow(counts),
        n_features_fit = nrow(dge),
        n_host_features_tested = length(retained_host),
        filter = list(
            method = "fixed_support",
            min_cpm = min_cpm,
            min_donors = min_donors,
            min_total_count = min_total_count
        ),
        normalisation = "TMM with complete-transcriptome library sizes"
    )
}

.null_build_design <- function(exposure, covariates, study,
                               min_group_donors = 10L) {
    min_group_donors <- .null_validate_positive_integer(
        min_group_donors,
        "min_group_donors",
        minimum = 2L
    )
    covariates <- as.data.frame(covariates, check.names = FALSE)
    if (length(exposure) != nrow(covariates) ||
            length(study) != nrow(covariates)) {
        stop("exposure, covariates and study must have matching rows.")
    }
    if (anyNA(exposure) || any(!is.finite(exposure)) ||
            anyNA(covariates) || anyNA(study)) {
        return(list(ok = FALSE, reason = "incomplete_design"))
    }
    exposure_levels <- sort(unique(exposure))
    if (length(exposure_levels) != 2L) {
        return(list(ok = FALSE, reason = "non_binary_exposure"))
    }
    group_counts <- table(factor(exposure, levels = exposure_levels))
    if (any(group_counts < min_group_donors)) {
        return(list(
            ok = FALSE,
            reason = "insufficient_permuted_group_size",
            group_counts = group_counts
        ))
    }
    if (!.null_within_study_variation(exposure, study)) {
        return(list(
            ok = FALSE,
            reason = "no_within_study_exposure_variation",
            group_counts = group_counts
        ))
    }

    design_data <- data.frame(
        exposure = as.numeric(exposure),
        covariates,
        check.names = FALSE
    )
    design_formula <- stats::reformulate(colnames(design_data))
    design_full <- stats::model.matrix(design_formula, data = design_data)
    protected <- match(c("(Intercept)", "exposure"), colnames(design_full))
    if (anyNA(protected) ||
            qr(design_full[, protected, drop = FALSE])$rank < 2L) {
        return(list(
            ok = FALSE,
            reason = "exposure_not_identifiable",
            group_counts = group_counts
        ))
    }

    retained <- protected
    nuisance <- setdiff(seq_len(ncol(design_full)), protected)
    for (column in nuisance) {
        candidate <- cbind(
            design_full[, retained, drop = FALSE],
            design_full[, column, drop = FALSE]
        )
        if (qr(candidate)$rank > length(retained)) {
            retained <- c(retained, column)
        }
    }
    design <- design_full[, retained, drop = FALSE]
    if (nrow(design) <= ncol(design) + 2L) {
        return(list(
            ok = FALSE,
            reason = "insufficient_residual_df",
            group_counts = group_counts
        ))
    }

    list(
        ok = TRUE,
        design = design,
        design_columns = colnames(design),
        dropped_design_columns = setdiff(
            colnames(design_full),
            colnames(design)
        ),
        group_counts = group_counts,
        residual_df = nrow(design) - ncol(design)
    )
}

.null_failed_celltype <- function(cache, reason, group_counts = NULL) {
    n_unexposed <- if (length(group_counts) == 2L) {
        unname(group_counts[[1]])
    } else {
        NA_integer_
    }
    n_exposed <- if (length(group_counts) == 2L) {
        unname(group_counts[[2]])
    } else {
        NA_integer_
    }
    list(
        result = data.frame(
            gene = cache$host_feature_ids,
            celltype = cache$celltype,
            pvalue = NA_real_,
            pvalue_underflow_clamped = FALSE,
            log2FC = NA_real_,
            failed = TRUE,
            failure_reason = reason,
            stringsAsFactors = FALSE
        ),
        diagnostic = data.frame(
            celltype = cache$celltype,
            n_donors = length(cache$donors),
            n_unexposed = n_unexposed,
            n_exposed = n_exposed,
            n_hypotheses = length(cache$host_feature_ids),
            n_failed = length(cache$host_feature_ids),
            n_pvalue_underflow_clamped = 0L,
            fit_failed = TRUE,
            failure_reason = reason,
            design_rank = NA_integer_,
            residual_df = NA_integer_,
            stringsAsFactors = FALSE
        )
    )
}

.null_fit_cached_celltype <- function(cache, permuted_exposure,
                                      min_group_donors = 10L,
                                      robust = TRUE) {
    if (!requireNamespace("edgeR", quietly = TRUE)) {
        stop("Package 'edgeR' is required.")
    }
    if (is.null(names(permuted_exposure)) ||
            any(!cache$donors %in% names(permuted_exposure))) {
        stop("permuted_exposure must be named for every cached donor.")
    }
    exposure <- as.numeric(permuted_exposure[cache$donors])
    design_result <- .null_build_design(
        exposure = exposure,
        covariates = cache$covariates,
        study = cache$study,
        min_group_donors = min_group_donors
    )
    if (!isTRUE(design_result$ok)) {
        return(.null_failed_celltype(
            cache,
            design_result$reason,
            design_result$group_counts
        ))
    }

    fit_result <- tryCatch({
        dge <- edgeR::DGEList(
            counts = cache$counts,
            lib.size = cache$full_library_size
        )
        dge$samples$norm.factors <- cache$normalisation_factors
        dge <- edgeR::estimateDisp(
            dge,
            design_result$design,
            robust = robust
        )
        fit <- edgeR::glmQLFit(
            dge,
            design_result$design,
            robust = robust
        )
        test <- edgeR::glmQLFTest(fit, coef = "exposure")
        table <- edgeR::topTags(test, n = Inf, sort.by = "none")$table
        list(fit = fit, table = table)
    }, error = function(error) {
        structure(
            list(message = conditionMessage(error)),
            class = "null_fit_error"
        )
    })
    if (inherits(fit_result, "null_fit_error")) {
        return(.null_failed_celltype(
            cache,
            paste0("edgeR_fit_error: ", fit_result$message),
            design_result$group_counts
        ))
    }

    host_index <- match(cache$host_feature_ids, rownames(fit_result$table))
    if (anyNA(host_index)) {
        return(.null_failed_celltype(
            cache,
            "host_hypothesis_missing_after_fit",
            design_result$group_counts
        ))
    }
    host_table <- fit_result$table[host_index, , drop = FALSE]
    failed <- !is.finite(host_table$PValue) |
        host_table$PValue < 0 | host_table$PValue > 1 |
        !is.finite(host_table$logFC)
    pvalue_underflow <- !failed & host_table$PValue == 0
    reported_pvalue <- host_table$PValue
    reported_pvalue[pvalue_underflow] <- .Machine$double.xmin
    failure_reason <- rep(NA_character_, nrow(host_table))
    failure_reason[failed] <- "non_finite_gene_result"
    result <- data.frame(
        gene = cache$host_feature_ids,
        celltype = cache$celltype,
        pvalue = ifelse(failed, NA_real_, reported_pvalue),
        pvalue_underflow_clamped = pvalue_underflow,
        log2FC = ifelse(failed, NA_real_, host_table$logFC),
        failed = failed,
        failure_reason = failure_reason,
        stringsAsFactors = FALSE
    )
    group_counts <- design_result$group_counts
    residual_df <- min(
        if (!is.null(fit_result$fit$df.residual.zeros)) {
            fit_result$fit$df.residual.zeros
        } else {
            fit_result$fit$df.residual
        },
        na.rm = TRUE
    )
    list(
        result = result,
        diagnostic = data.frame(
            celltype = cache$celltype,
            n_donors = length(cache$donors),
            n_unexposed = unname(group_counts[[1]]),
            n_exposed = unname(group_counts[[2]]),
            n_hypotheses = nrow(result),
            n_failed = sum(failed),
            n_pvalue_underflow_clamped = sum(pvalue_underflow),
            fit_failed = FALSE,
            failure_reason = if (any(failed)) {
                "one_or_more_non_finite_gene_results"
            } else {
                NA_character_
            },
            design_rank = qr(design_result$design)$rank,
            residual_df = residual_df,
            stringsAsFactors = FALSE
        )
    )
}

.null_apply_family_bh <- function(result) {
    required <- c("pvalue", "failed")
    if (!all(required %in% colnames(result))) {
        stop("result must contain pvalue and failed.")
    }
    if (anyNA(result$failed)) {
        stop("result$failed must not contain missing values.")
    }
    valid <- !result$failed & is.finite(result$pvalue) &
        result$pvalue >= 0 & result$pvalue <= 1
    result$padj_global <- NA_real_
    if (any(valid)) {
        ## n is the complete, pre-specified host-expression family. Failed
        ## hypotheses therefore cannot make the multiplicity correction less
        ## stringent by disappearing from the denominator.
        result$padj_global[valid] <- stats::p.adjust(
            result$pvalue[valid],
            method = "BH",
            n = nrow(result)
        )
    }
    result
}

.null_replicate_diagnostics <- function(result, replicate_id, seed,
                                        p_thresholds = c(0.001, 0.01, 0.05),
                                        q_thresholds = c(0.01, 0.05, 0.10),
                                        qq_probabilities = c(
                                            0.001, 0.005, 0.01, 0.05, 0.10,
                                            0.25, 0.50, 0.75, 0.90, 0.95,
                                            0.99, 0.995, 0.999
                                        )) {
    replicate_id <- .null_validate_positive_integer(
        replicate_id,
        "replicate_id"
    )
    seed <- .null_validate_positive_integer(seed, "seed")
    valid <- !result$failed & is.finite(result$pvalue) &
        result$pvalue >= 0 & result$pvalue <= 1
    pvalue <- result$pvalue[valid]
    n_family <- nrow(result)
    n_valid <- length(pvalue)
    lambda <- if (n_valid) {
        bounded <- pmin(
            pmax(pvalue, .Machine$double.xmin),
            1 - .Machine$double.eps
        )
        stats::median(stats::qchisq(
            bounded,
            df = 1,
            lower.tail = FALSE
        )) /
            stats::qchisq(0.5, df = 1)
    } else {
        NA_real_
    }
    summary <- data.frame(
        replicate_id = replicate_id,
        seed = seed,
        n_family = n_family,
        n_valid = n_valid,
        n_failed = n_family - n_valid,
        failed_fraction = (n_family - n_valid) / n_family,
        lambda_gc = lambda,
        stringsAsFactors = FALSE
    )
    marginal <- do.call(rbind, lapply(p_thresholds, function(threshold) {
        n_below <- sum(pvalue <= threshold)
        data.frame(
            replicate_id = replicate_id,
            threshold = threshold,
            n_below = n_below,
            rate_valid = if (n_valid) n_below / n_valid else NA_real_,
            rate_family = n_below / n_family,
            stringsAsFactors = FALSE
        )
    }))
    calls <- do.call(rbind, lapply(q_thresholds, function(threshold) {
        n_called <- sum(
            is.finite(result$padj_global) &
                result$padj_global <= threshold
        )
        data.frame(
            replicate_id = replicate_id,
            threshold = threshold,
            n_called = n_called,
            any_called = n_called > 0,
            stringsAsFactors = FALSE
        )
    }))
    qq <- data.frame(
        replicate_id = replicate_id,
        probability = qq_probabilities,
        observed_p = if (n_valid) {
            as.numeric(stats::quantile(
                pvalue,
                probs = qq_probabilities,
                names = FALSE,
                type = 8
            ))
        } else {
            rep(NA_real_, length(qq_probabilities))
        },
        stringsAsFactors = FALSE
    )
    list(summary = summary, marginal = marginal, calls = calls, qq = qq)
}

.null_wilson_interval <- function(successes, trials, confidence = 0.95) {
    if (!is.numeric(successes) || length(successes) != 1L ||
            !is.numeric(trials) || length(trials) != 1L ||
            !is.finite(successes) || !is.finite(trials) ||
            trials < 1 || successes < 0 || successes > trials) {
        stop("successes and trials must satisfy 0 <= successes <= trials.")
    }
    if (!is.numeric(confidence) || length(confidence) != 1L ||
            !is.finite(confidence) || confidence <= 0 || confidence >= 1) {
        stop("confidence must be one number strictly between zero and one.")
    }
    z <- stats::qnorm(1 - (1 - confidence) / 2)
    proportion <- successes / trials
    denominator <- 1 + z^2 / trials
    centre <- (proportion + z^2 / (2 * trials)) / denominator
    half_width <- z * sqrt(
        proportion * (1 - proportion) / trials +
            z^2 / (4 * trials^2)
    ) / denominator
    c(
        estimate = proportion,
        lower = max(0, centre - half_width),
        upper = min(1, centre + half_width)
    )
}

.null_quantile_distribution <- function(value) {
    value <- value[is.finite(value)]
    if (!length(value)) {
        return(data.frame(
            n_replicates = 0L,
            mean = NA_real_,
            sd = NA_real_,
            median = NA_real_,
            reference_interval_lower = NA_real_,
            reference_interval_upper = NA_real_,
            minimum = NA_real_,
            maximum = NA_real_
        ))
    }
    interval <- stats::quantile(
        value,
        probs = c(0.025, 0.975),
        names = FALSE,
        type = 8
    )
    data.frame(
        n_replicates = length(value),
        mean = mean(value),
        sd = if (length(value) > 1L) stats::sd(value) else NA_real_,
        median = stats::median(value),
        reference_interval_lower = interval[[1]],
        reference_interval_upper = interval[[2]],
        minimum = min(value),
        maximum = max(value)
    )
}

.null_summarise_replicates <- function(receipts) {
    if (!length(receipts)) {
        stop("At least one completed replicate receipt is required.")
    }
    summary <- do.call(rbind, lapply(receipts, `[[`, "summary"))
    marginal <- do.call(rbind, lapply(receipts, `[[`, "marginal"))
    calls <- do.call(rbind, lapply(receipts, `[[`, "calls"))
    qq <- do.call(rbind, lapply(receipts, `[[`, "qq"))

    familywise <- do.call(rbind, lapply(
        sort(unique(calls$threshold)),
        function(threshold) {
            selected <- calls[calls$threshold == threshold, , drop = FALSE]
            interval <- .null_wilson_interval(
                sum(selected$any_called),
                nrow(selected)
            )
            data.frame(
                q_threshold = threshold,
                n_replicates = nrow(selected),
                n_any_discovery = sum(selected$any_called),
                probability_any_discovery = unname(interval[["estimate"]]),
                wilson_lower = unname(interval[["lower"]]),
                wilson_upper = unname(interval[["upper"]]),
                stringsAsFactors = FALSE
            )
        }
    ))
    call_distribution <- do.call(rbind, lapply(
        sort(unique(calls$threshold)),
        function(threshold) {
            selected <- calls$n_called[calls$threshold == threshold]
            cbind(
                q_threshold = threshold,
                .null_quantile_distribution(selected)
            )
        }
    ))
    marginal_distribution <- do.call(rbind, lapply(
        sort(unique(marginal$threshold)),
        function(threshold) {
            selected <- marginal[marginal$threshold == threshold, ]
            valid_summary <- .null_quantile_distribution(selected$rate_valid)
            family_summary <- .null_quantile_distribution(selected$rate_family)
            rbind(
                cbind(
                    p_threshold = threshold,
                    denominator = "successfully_fitted_hypotheses",
                    valid_summary
                ),
                cbind(
                    p_threshold = threshold,
                    denominator = "fixed_host_family",
                    family_summary
                )
            )
        }
    ))
    qq_distribution <- do.call(rbind, lapply(
        sort(unique(qq$probability)),
        function(probability) {
            selected <- qq$observed_p[qq$probability == probability]
            cbind(
                expected_uniform_quantile = probability,
                .null_quantile_distribution(selected)
            )
        }
    ))

    list(
        interpretation = paste(
            "Study-stratified empirical diagnostic only; replicates, not",
            "gene-by-replicate rows, are the independent calibration units."
        ),
        replicate_summary = summary,
        marginal_rejection_distribution = marginal_distribution,
        familywise_error = familywise,
        discovery_count_distribution = call_distribution,
        lambda_distribution = .null_quantile_distribution(summary$lambda_gc),
        permutation_attempt_distribution = if (
            "permutation_attempts" %in% colnames(summary)
        ) {
            .null_quantile_distribution(summary$permutation_attempts)
        } else {
            .null_quantile_distribution(numeric())
        },
        failed_fraction_distribution = .null_quantile_distribution(
            summary$failed_fraction
        ),
        qq_reference_distribution = qq_distribution,
        replicate_marginal = marginal,
        replicate_calls = calls,
        replicate_qq = qq
    )
}

.null_atomic_save_rds <- function(object, path, compress = "gzip",
                                  overwrite = FALSE) {
    directory <- dirname(path)
    dir.create(directory, recursive = TRUE, showWarnings = FALSE)
    temporary <- tempfile(
        pattern = paste0(".", basename(path), ".tmp-"),
        tmpdir = directory
    )
    on.exit(unlink(temporary), add = TRUE)
    saveRDS(object, temporary, compress = compress)
    if (file.exists(path)) {
        if (!isTRUE(overwrite)) {
            stop("Refusing to overwrite existing file: ", path)
        }
        backup <- paste0(
            path,
            ".superseded-",
            format(Sys.time(), "%Y%m%dT%H%M%S", tz = "UTC")
        )
        if (!file.rename(path, backup)) {
            stop("Could not preserve the existing file before replacement.")
        }
    }
    if (!file.rename(temporary, path)) {
        stop("Atomic rename failed for: ", path)
    }
    invisible(path)
}

.null_replicate_path <- function(replicate_dir, replicate_id) {
    replicate_id <- .null_validate_positive_integer(
        replicate_id,
        "replicate_id"
    )
    file.path(
        replicate_dir,
        sprintf("replicate_%06d.rds", replicate_id)
    )
}

.null_receipt_path <- function(replicate_path) {
    sub("[.]rds$", ".receipt.rds", replicate_path)
}

.null_save_replicate <- function(payload, replicate_path) {
    if (is.null(payload$replicate_id) || is.null(payload$seed) ||
            is.null(payload$result) || is.null(payload$summary) ||
            is.null(payload$marginal) || is.null(payload$calls) ||
            is.null(payload$qq)) {
        stop("payload is missing a required replicate component.")
    }
    if (!is.null(payload$replicate_error) ||
            !is.data.frame(payload$summary) || nrow(payload$summary) != 1L ||
            !"n_valid" %in% colnames(payload$summary) ||
            !is.finite(payload$summary$n_valid[[1]]) ||
            payload$summary$n_valid[[1]] < 1L) {
        stop(
            "A failed or zero-valid-hypothesis replicate cannot receive a ",
            "completion receipt."
        )
    }
    .null_atomic_save_rds(payload, replicate_path, compress = "gzip")
    receipt <- list(
        schema_version = "exposomeSC_null_receipt_v1",
        replicate_id = as.integer(payload$replicate_id),
        seed = as.integer(payload$seed),
        result_file = basename(replicate_path),
        bytes = unname(file.info(replicate_path)$size),
        md5 = unname(tools::md5sum(replicate_path)),
        completed_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
        summary = payload$summary,
        marginal = payload$marginal,
        calls = payload$calls,
        qq = payload$qq
    )
    .null_atomic_save_rds(
        receipt,
        .null_receipt_path(replicate_path),
        compress = "gzip"
    )
    receipt
}

.null_create_manifest <- function(contract, n_replicates, master_seed) {
    n_replicates <- .null_validate_positive_integer(
        n_replicates,
        "n_replicates"
    )
    master_seed <- .null_validate_positive_integer(
        master_seed,
        "master_seed"
    )
    list(
        schema_version = "exposomeSC_null_manifest_v1",
        contract = contract,
        n_replicates = n_replicates,
        master_seed = master_seed,
        replicate_seeds = .null_replicate_seeds(
            master_seed,
            n_replicates
        ),
        created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
        updated_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
        completed_replicates = integer(),
        state = "initialised"
    )
}

.null_validate_manifest <- function(manifest, contract, n_replicates,
                                    master_seed) {
    if (!identical(
            manifest$schema_version,
            "exposomeSC_null_manifest_v1"
        ) || !identical(manifest$contract, contract) ||
            !identical(manifest$n_replicates, as.integer(n_replicates)) ||
            !identical(manifest$master_seed, as.integer(master_seed))) {
        stop(
            "Existing null-calibration manifest is incompatible with the ",
            "requested contract, replicate count or master seed."
        )
    }
    expected <- .null_replicate_seeds(master_seed, n_replicates)
    if (!identical(manifest$replicate_seeds, expected)) {
        stop("Existing manifest has an invalid replicate-seed schedule.")
    }
    invisible(TRUE)
}

.null_manifest_pending <- function(manifest, replicate_dir,
                                   replicate_ids = NULL,
                                   check_md5 = FALSE) {
    if (is.null(replicate_ids)) {
        replicate_ids <- seq_len(manifest$n_replicates)
    }
    replicate_ids <- as.integer(replicate_ids)
    if (anyNA(replicate_ids) || any(replicate_ids < 1L) ||
            any(replicate_ids > manifest$n_replicates) ||
            anyDuplicated(replicate_ids)) {
        stop("replicate_ids are invalid for this manifest.")
    }
    complete <- vapply(replicate_ids, function(replicate_id) {
        result_path <- .null_replicate_path(replicate_dir, replicate_id)
        receipt_path <- .null_receipt_path(result_path)
        if (!file.exists(result_path) || !file.exists(receipt_path)) {
            return(FALSE)
        }
        receipt <- tryCatch(readRDS(receipt_path), error = function(error) NULL)
        if (is.null(receipt) ||
                !identical(receipt$schema_version,
                    "exposomeSC_null_receipt_v1") ||
                !identical(receipt$replicate_id, as.integer(replicate_id)) ||
                !identical(
                    receipt$seed,
                    as.integer(manifest$replicate_seeds[[replicate_id]])
                ) || !identical(
                    as.numeric(receipt$bytes),
                    as.numeric(file.info(result_path)$size)
                ) || !is.data.frame(receipt$summary) ||
                nrow(receipt$summary) != 1L ||
                !"n_valid" %in% colnames(receipt$summary) ||
                !is.finite(receipt$summary$n_valid[[1]]) ||
                receipt$summary$n_valid[[1]] < 1L) {
            return(FALSE)
        }
        if (isTRUE(check_md5) && !identical(
            unname(tools::md5sum(result_path)),
            receipt$md5
        )) {
            return(FALSE)
        }
        TRUE
    }, logical(1))
    replicate_ids[!complete]
}
