# Internal helpers for sample-level pseudobulk sensitivity analyses.
#
# These functions deliberately remain unexported. The canonical sc-ExWAS
# estimand uses one pseudobulk column per donor. The helpers below support a
# secondary sample-by-cell-type analysis in which repeated biospecimens are
# correlated within donor.

.sample_voom_lm_fit_backend <- function() {
    # voomLmFit() moved from edgeR to limma during the Bioconductor 3.24 devel
    # cycle. Resolve it at run time so the same source works against either
    # package, and use getExportedValue() rather than :: so that R CMD check
    # does not flag whichever of the two does not currently export it.
    for (package in c("limma", "edgeR")) {
        if (requireNamespace(package, quietly = TRUE) &&
                "voomLmFit" %in% getNamespaceExports(asNamespace(package))) {
            return(getExportedValue(package, "voomLmFit"))
        }
    }
    stop(
        "voomLmFit() is exported by neither limma nor edgeR; ",
        "install limma (>= 3.99.0) or edgeR (< 4.99.0)."
    )
}

.sample_voom_positive_integer <- function(value, name) {
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
            value < 1 || value != as.integer(value)) {
        stop(name, " must be one positive integer.")
    }
    as.integer(value)
}

.sample_voom_site_class <- function(tissue_level_2) {
    tissue <- as.character(tissue_level_2)
    output <- rep("unknown", length(tissue))
    observed <- !is.na(tissue) & nzchar(tissue)
    output[observed & grepl("^parenchyma", tissue)] <- "parenchyma"
    output[observed & grepl(
        "lobular bronchi|segmental bronchi|distal lobular airways",
        tissue
    )] <- "lower_airway"
    output[observed & grepl("trachea|inferior turbinate", tissue)] <-
        "upper_airway"
    factor(
        output,
        levels = c(
            "parenchyma",
            "lower_airway",
            "upper_airway",
            "unknown"
        )
    )
}

.sample_voom_aggregate <- function(counts, donor_id, sample_id, celltype,
                                   target_celltype,
                                   min_donor_cells = 20L,
                                   min_sample_cells = 10L) {
    if (!inherits(counts, "Matrix") && !is.matrix(counts)) {
        stop("counts must be a matrix or Matrix object.")
    }
    if (is.null(rownames(counts)) || anyDuplicated(rownames(counts))) {
        stop("counts must have unique gene row names.")
    }
    n_cells <- ncol(counts)
    fields <- list(
        donor_id = donor_id,
        sample_id = sample_id,
        celltype = celltype
    )
    if (any(vapply(fields, length, integer(1)) != n_cells)) {
        stop("Cell annotations must have length ncol(counts).")
    }
    if (anyNA(donor_id) || anyNA(sample_id) || anyNA(celltype)) {
        stop("Cell annotations must not contain missing values.")
    }
    min_donor_cells <- .sample_voom_positive_integer(
        min_donor_cells,
        "min_donor_cells"
    )
    min_sample_cells <- .sample_voom_positive_integer(
        min_sample_cells,
        "min_sample_cells"
    )
    if (!is.character(target_celltype) || length(target_celltype) != 1L ||
            is.na(target_celltype) || !nzchar(target_celltype)) {
        stop("target_celltype must be one non-empty string.")
    }

    donor_id <- as.character(donor_id)
    sample_id <- as.character(sample_id)
    celltype <- as.character(celltype)
    sample_donor_n <- tapply(
        donor_id,
        sample_id,
        function(value) length(unique(value))
    )
    if (any(sample_donor_n != 1L)) {
        stop("Every sample_id must map to exactly one donor_id.")
    }

    target_index <- which(celltype == target_celltype)
    if (!length(target_index)) {
        return(NULL)
    }
    donor_cell_counts_all <- table(donor_id[target_index])
    eligible_donors <- names(
        donor_cell_counts_all[donor_cell_counts_all >= min_donor_cells]
    )
    if (!length(eligible_donors)) {
        return(NULL)
    }

    donor_eligible_index <- target_index[
        donor_id[target_index] %in% eligible_donors
    ]
    sample_cell_counts_all <- table(sample_id[donor_eligible_index])
    eligible_samples <- sort(names(
        sample_cell_counts_all[
            sample_cell_counts_all >= min_sample_cells
        ]
    ))
    retained_index <- donor_eligible_index[
        sample_id[donor_eligible_index] %in% eligible_samples
    ]
    if (length(eligible_samples)) {
        sample_factor <- factor(
            sample_id[retained_index],
            levels = eligible_samples
        )
        membership <- Matrix::sparseMatrix(
            i = seq_along(retained_index),
            j = as.integer(sample_factor),
            x = 1,
            dims = c(length(retained_index), length(eligible_samples)),
            dimnames = list(NULL, eligible_samples)
        )
        pseudobulk <- counts[, retained_index, drop = FALSE] %*% membership
        rownames(pseudobulk) <- rownames(counts)
        colnames(pseudobulk) <- eligible_samples
        sample_donor <- vapply(eligible_samples, function(sample) {
            unique(donor_id[retained_index][
                sample_id[retained_index] == sample
            ])
        }, character(1))
        sample_cell_counts <- as.integer(
            sample_cell_counts_all[eligible_samples]
        )
        names(sample_cell_counts) <- eligible_samples
    } else {
        pseudobulk <- counts[, integer(), drop = FALSE]
        sample_donor <- setNames(character(), character())
        sample_cell_counts <- setNames(integer(), character())
    }
    donor_cell_counts <- as.integer(donor_cell_counts_all[eligible_donors])
    names(donor_cell_counts) <- eligible_donors

    list(
        counts = pseudobulk,
        sample_ids = eligible_samples,
        sample_donor = sample_donor,
        sample_cell_counts = sample_cell_counts,
        donor_cell_counts = donor_cell_counts,
        eligible_donors = eligible_donors,
        donors_without_eligible_sample = setdiff(
            eligible_donors,
            unname(sample_donor)
        ),
        retained_cell_index = retained_index,
        min_donor_cells = min_donor_cells,
        min_sample_cells = min_sample_cells
    )
}

.sample_voom_overlap <- function(metadata) {
    required <- c("donor_id", "study", "exposure")
    missing_fields <- setdiff(required, colnames(metadata))
    if (length(missing_fields)) {
        stop(
            "metadata is missing: ",
            paste(missing_fields, collapse = ", ")
        )
    }
    if (!is.numeric(metadata$exposure) && !is.logical(metadata$exposure)) {
        stop("overlap_only exposure must be numeric or logical.")
    }
    if (anyNA(metadata$donor_id) || anyNA(metadata$study) ||
            any(!nzchar(as.character(metadata$donor_id))) ||
            any(!nzchar(as.character(metadata$study)))) {
        stop("overlap_only requires complete donor and study identifiers.")
    }
    exposure <- as.numeric(metadata$exposure)
    if (anyNA(exposure) || !all(exposure %in% c(0, 1))) {
        stop("overlap_only requires a complete binary exposure.")
    }
    donor_split <- split(seq_len(nrow(metadata)), metadata$donor_id)
    donor_constant <- vapply(donor_split, function(index) {
        length(unique(exposure[index])) == 1L &&
            length(unique(metadata$study[index])) == 1L
    }, logical(1))
    if (!all(donor_constant)) {
        stop("Exposure and study must be constant within donor.")
    }
    donor_index <- vapply(donor_split, `[[`, integer(1), 1L)
    donor_data <- data.frame(
        donor_id = names(donor_index),
        study = as.character(metadata$study[donor_index]),
        exposure = exposure[donor_index],
        stringsAsFactors = FALSE
    )
    exposure_by_study <- table(
        donor_data$study,
        factor(donor_data$exposure, levels = c(0, 1))
    )
    overlap_studies <- rownames(exposure_by_study)[
        rowSums(exposure_by_study > 0) == 2L
    ]
    retained <- as.character(metadata$study) %in% overlap_studies
    list(
        metadata = metadata[retained, , drop = FALSE],
        overlap_studies = sort(overlap_studies),
        excluded_studies = sort(setdiff(
            unique(as.character(metadata$study)),
            overlap_studies
        )),
        donor_counts = exposure_by_study
    )
}

.sample_voom_factor <- function(value, levels = NULL, reference = NULL,
                                name = "factor") {
    value <- as.character(value)
    if (anyNA(value) || any(!nzchar(value))) {
        stop(name, " must contain complete, non-empty values.")
    }
    if (is.null(levels)) {
        levels <- sort(unique(value))
    }
    if (!is.character(levels) || !length(levels) || anyNA(levels) ||
            any(!nzchar(levels)) || anyDuplicated(levels)) {
        stop(name, " levels must be unique, non-empty strings.")
    }
    unknown <- setdiff(unique(value), levels)
    if (length(unknown)) {
        stop(
            name, " contains value(s) outside the fixed levels: ",
            paste(unknown, collapse = ", ")
        )
    }
    if (is.null(reference)) {
        reference <- levels[[1L]]
    }
    if (!is.character(reference) || length(reference) != 1L ||
            is.na(reference) || !reference %in% levels) {
        stop(name, " reference must be one of its fixed levels.")
    }
    factor(value, levels = c(reference, setdiff(levels, reference)))
}

.sample_voom_constant_within_donor <- function(metadata, fields) {
    donor_split <- split(seq_len(nrow(metadata)), metadata$donor_id)
    for (field in fields) {
        constant <- vapply(donor_split, function(index) {
            length(unique(metadata[[field]][index])) == 1L
        }, logical(1))
        if (!all(constant)) {
            stop(
                field, " must be constant within donor; affected donor(s): ",
                paste(names(constant)[!constant], collapse = ", ")
            )
        }
    }
    invisible(TRUE)
}

.sample_voom_design <- function(metadata,
                                variant = c(
                                    "mirror",
                                    "site",
                                    "dataset_site",
                                    "overlap_only"
                                ),
                                study_levels = NULL,
                                study_reference = NULL,
                                dataset_levels = NULL,
                                dataset_reference = NULL) {
    variant <- match.arg(variant)
    age_terms <- paste0("age_ns", seq_len(3L))
    grouping_term <- if (variant == "dataset_site") "dataset" else "study"
    site_term <- if (variant %in% c(
        "site", "dataset_site", "overlap_only"
    )) {
        "site_class"
    } else {
        character()
    }
    required <- c(
        "donor_id",
        "exposure",
        age_terms,
        "sex_binary",
        grouping_term,
        site_term
    )
    missing_fields <- setdiff(required, colnames(metadata))
    if (length(missing_fields)) {
        stop(
            "metadata is missing design field(s): ",
            paste(missing_fields, collapse = ", ")
        )
    }
    if (is.null(rownames(metadata)) || anyDuplicated(rownames(metadata))) {
        stop("metadata must have unique sample row names.")
    }
    if (anyNA(rownames(metadata)) || any(!nzchar(rownames(metadata)))) {
        stop("metadata sample row names must be complete and non-empty.")
    }
    original_samples <- rownames(metadata)
    overlap_studies <- character()
    excluded_studies <- character()
    if (variant == "overlap_only") {
        overlap <- .sample_voom_overlap(metadata)
        metadata <- overlap$metadata
        overlap_studies <- overlap$overlap_studies
        excluded_studies <- overlap$excluded_studies
        if (!nrow(metadata)) {
            stop("No study contains donors from both exposure groups.")
        }
    }
    complete <- stats::complete.cases(
        metadata[, required, drop = FALSE]
    )
    if (!all(complete)) {
        stop(
            "The sample-level design contains incomplete required fields: ",
            paste(rownames(metadata)[!complete], collapse = ", ")
        )
    }
    metadata <- metadata[complete, , drop = FALSE]
    if (!is.numeric(metadata$exposure) &&
            !is.logical(metadata$exposure)) {
        stop("Sample-blocked exposure must be numeric or logical.")
    }
    metadata$exposure <- as.numeric(metadata$exposure)
    if (!all(metadata$exposure %in% c(0, 1))) {
        stop("Sample-blocked analysis requires a binary exposure.")
    }
    .sample_voom_constant_within_donor(
        metadata,
        unique(c(
            "exposure",
            age_terms,
            "sex_binary",
            "study",
            grouping_term
        ))
    )
    if (grouping_term == "study") {
        metadata$study <- .sample_voom_factor(
            metadata$study,
            levels = study_levels,
            reference = study_reference,
            name = "study"
        )
    } else {
        metadata$dataset <- .sample_voom_factor(
            metadata$dataset,
            levels = dataset_levels,
            reference = dataset_reference,
            name = "dataset"
        )
    }
    if (length(site_term)) {
        metadata$site_class <- .sample_voom_factor(
            metadata$site_class,
            levels = c(
                "parenchyma",
                "lower_airway",
                "upper_airway",
                "unknown"
            ),
            reference = "parenchyma",
            name = "site_class"
        )
    }

    formula_terms <- c(
        "exposure",
        age_terms,
        "sex_binary",
        grouping_term,
        site_term
    )
    design_formula <- stats::reformulate(formula_terms)
    design_full <- stats::model.matrix(design_formula, data = metadata)
    mandatory_names <- c(
        "(Intercept)",
        "exposure",
        age_terms,
        "sex_binary"
    )
    mandatory <- match(mandatory_names, colnames(design_full))
    if (anyNA(mandatory)) {
        stop("The mandatory primary-design columns were not constructed.")
    }
    mandatory_design <- design_full[, mandatory, drop = FALSE]
    if (qr(mandatory_design)$rank < ncol(mandatory_design)) {
        stop(
            "Exposure, the adult age spline, and sex are not jointly ",
            "identifiable; the primary adjustment set cannot be weakened."
        )
    }
    retained_columns <- mandatory
    nuisance <- setdiff(seq_len(ncol(design_full)), mandatory)
    for (column in nuisance) {
        candidate <- cbind(
            design_full[, retained_columns, drop = FALSE],
            design_full[, column, drop = FALSE]
        )
        if (qr(candidate)$rank > length(retained_columns)) {
            retained_columns <- c(retained_columns, column)
        }
    }
    design <- design_full[, retained_columns, drop = FALSE]
    if (qr(design)$rank != ncol(design)) {
        stop("Internal error: the retained design is not full rank.")
    }
    if (nrow(design) <= ncol(design) + 2L) {
        stop("Insufficient residual degrees of freedom.")
    }

    list(
        metadata = metadata,
        design = design,
        formula = design_formula,
        design_columns = colnames(design),
        dropped_design_columns = setdiff(
            colnames(design_full),
            colnames(design)
        ),
        design_rank = qr(design)$rank,
        residual_df = nrow(design) - ncol(design),
        excluded_samples = setdiff(original_samples, rownames(metadata)),
        overlap_studies = overlap_studies,
        excluded_studies = excluded_studies,
        grouping_levels = levels(metadata[[grouping_term]]),
        grouping_reference = levels(metadata[[grouping_term]])[[1L]],
        variant = variant
    )
}

.sample_voom_filter <- function(counts, donor_id, design,
                                min_donor_support = 10L,
                                min_cpm = 1,
                                filter_method = c(
                                    "fixed_support",
                                    "fixed_support_plus_filterByExpr"
                                ),
                                filter_min_count = 10L,
                                filter_min_total_count = 15L) {
    if (!requireNamespace("edgeR", quietly = TRUE)) {
        stop("Package 'edgeR' is required.")
    }
    counts <- as.matrix(counts)
    if (is.null(rownames(counts)) || is.null(colnames(counts)) ||
            anyNA(rownames(counts)) || anyNA(colnames(counts)) ||
            any(!nzchar(rownames(counts))) ||
            any(!nzchar(colnames(counts))) ||
            anyDuplicated(rownames(counts)) ||
            anyDuplicated(colnames(counts))) {
        stop("counts must have unique, complete gene and sample dimnames.")
    }
    if (length(donor_id) != ncol(counts)) {
        stop("donor_id must have length ncol(counts).")
    }
    if (!identical(colnames(counts), rownames(design))) {
        stop("counts columns must match design rows in order.")
    }
    if (!is.matrix(design) || !is.numeric(design) ||
            any(!is.finite(design)) || qr(design)$rank != ncol(design)) {
        stop("design must be a finite, full-rank numeric matrix.")
    }
    if (any(!is.finite(counts)) || any(counts < 0)) {
        stop("counts must be finite and non-negative.")
    }
    if (anyNA(donor_id) || any(!nzchar(as.character(donor_id)))) {
        stop("donor_id must contain complete, non-empty values.")
    }
    if (!is.numeric(min_donor_support) ||
            length(min_donor_support) != 1L ||
            !is.finite(min_donor_support) || min_donor_support < 1 ||
            min_donor_support != as.integer(min_donor_support)) {
        stop("min_donor_support must be one positive integer.")
    }
    if (!is.numeric(min_cpm) || length(min_cpm) != 1L ||
            !is.finite(min_cpm) || min_cpm <= 0) {
        stop("min_cpm must be one finite number > 0.")
    }
    if (!is.numeric(filter_min_count) || length(filter_min_count) != 1L ||
            !is.finite(filter_min_count) || filter_min_count < 0) {
        stop("filter_min_count must be one finite non-negative number.")
    }
    if (!is.numeric(filter_min_total_count) ||
            length(filter_min_total_count) != 1L ||
            !is.finite(filter_min_total_count) ||
            filter_min_total_count < 0) {
        stop(
            "filter_min_total_count must be one finite non-negative number."
        )
    }
    filter_method <- match.arg(filter_method)
    donor_id <- as.character(donor_id)
    donor_factor <- factor(donor_id, levels = unique(donor_id))
    donor_counts <- t(rowsum(
        t(counts),
        group = donor_factor,
        reorder = FALSE
    ))
    donor_library_size <- pmax(colSums(donor_counts), 1)
    donor_dge <- edgeR::DGEList(
        counts = donor_counts,
        lib.size = donor_library_size
    )
    donor_cpm <- edgeR::cpm(
        donor_dge,
        normalized.lib.sizes = FALSE,
        log = FALSE
    )
    donor_support <- rowSums(donor_cpm >= min_cpm)
    support_keep <- donor_support >= min_donor_support

    full_library_size <- colSums(counts)
    if (any(full_library_size <= 0)) {
        stop("Every sample pseudobulk must have a positive library size.")
    }
    dge <- edgeR::DGEList(
        counts = counts,
        lib.size = full_library_size
    )
    total_count_keep <- rowSums(counts) >= filter_min_total_count
    edger_keep <- if (filter_method == "fixed_support_plus_filterByExpr") {
        edgeR::filterByExpr(
            dge,
            design = design,
            min.count = filter_min_count,
            min.total.count = filter_min_total_count
        )
    } else {
        rep(NA, nrow(counts))
    }
    keep <- support_keep & total_count_keep
    if (filter_method == "fixed_support_plus_filterByExpr") {
        keep <- keep & edger_keep
    }
    if (!any(keep)) {
        stop("No gene passed the sample-level expression filter.")
    }

    list(
        dge = dge,
        keep = keep,
        support_keep = support_keep,
        total_count_keep = total_count_keep,
        edger_keep = edger_keep,
        donor_support = donor_support,
        full_library_size = full_library_size,
        donor_library_size = donor_library_size,
        min_donor_support = as.integer(min_donor_support),
        min_cpm = min_cpm,
        filter_method = filter_method,
        filter_min_count = as.integer(filter_min_count),
        filter_min_total_count = as.integer(filter_min_total_count)
    )
}

.sample_voom_fit <- function(counts, metadata, design,
                             robust = TRUE,
                             min_donor_support = 10L,
                             min_cpm = 1,
                             filter_method = c(
                                 "fixed_support",
                                 "fixed_support_plus_filterByExpr"
                             ),
                             filter_min_count = 10L,
                             filter_min_total_count = 15L,
                             test_features = NULL,
                             exposure_name = "smoking_ever") {
    if (!requireNamespace("limma", quietly = TRUE)) {
        stop("Package 'limma' is required.")
    }
    if (!is.logical(robust) || length(robust) != 1L || is.na(robust)) {
        stop("robust must be TRUE or FALSE.")
    }
    if (isTRUE(robust) &&
            !requireNamespace("statmod", quietly = TRUE)) {
        stop("Package 'statmod' is required for robust limma inference.")
    }
    if (!identical(colnames(counts), rownames(metadata)) ||
            !identical(rownames(metadata), rownames(design))) {
        stop("counts, metadata, and design must have identical sample order.")
    }
    if (!"donor_id" %in% colnames(metadata)) {
        stop("metadata must contain donor_id.")
    }
    if (!is.character(exposure_name) || length(exposure_name) != 1L ||
            is.na(exposure_name) || !nzchar(exposure_name)) {
        stop("exposure_name must be one non-empty string.")
    }
    feature_ids <- rownames(counts)
    if (!is.null(test_features)) {
        if (!is.character(test_features) || !length(test_features) ||
                anyNA(test_features) || any(!nzchar(test_features)) ||
                anyDuplicated(test_features)) {
            stop(
                "test_features must be NULL or a non-empty vector of ",
                "unique feature identifiers."
            )
        }
        missing_features <- setdiff(test_features, feature_ids)
        if (length(missing_features)) {
            stop(
                "test_features not found in counts: ",
                paste(missing_features, collapse = ", ")
            )
        }
    }

    filtered <- .sample_voom_filter(
        counts = counts,
        donor_id = metadata$donor_id,
        design = design,
        min_donor_support = min_donor_support,
        min_cpm = min_cpm,
        filter_method = match.arg(filter_method),
        filter_min_count = filter_min_count,
        filter_min_total_count = filter_min_total_count
    )
    requested <- if (is.null(test_features)) feature_ids else test_features
    requested_index <- match(requested, feature_ids)
    requested_support <- filtered$support_keep[requested_index]
    requested_total <- filtered$total_count_keep[requested_index]
    requested_edger <- filtered$edger_keep[requested_index]
    requested_estimable <- filtered$keep[requested_index]
    not_estimable_reason <- rep(NA_character_, length(requested))
    not_estimable_reason[!requested_support & !requested_total] <-
        "below_donor_support_and_total_count"
    not_estimable_reason[!requested_support & requested_total] <-
        "below_donor_support"
    not_estimable_reason[requested_support & !requested_total] <-
        "below_total_count"
    if (filtered$filter_method == "fixed_support_plus_filterByExpr") {
        filter_by_expr_failure <- requested_support & requested_total &
            !requested_edger
        not_estimable_reason[filter_by_expr_failure] <-
            "below_filterByExpr"
    }
    requested_audit <- data.frame(
        gene = requested,
        donor_support = filtered$donor_support[requested_index],
        passed_donor_support = requested_support,
        passed_total_count = requested_total,
        passed_filter_by_expr = requested_edger,
        estimable = requested_estimable,
        not_estimable_reason = not_estimable_reason,
        stringsAsFactors = FALSE
    )
    if (sum(filtered$keep) < 2L) {
        stop("At least two genes must pass filtering for voomLmFit.")
    }
    dge <- filtered$dge[
        filtered$keep,
        ,
        keep.lib.sizes = TRUE
    ]
    dge <- edgeR::calcNormFactors(dge, method = "TMM")
    donor_block <- factor(metadata$donor_id)
    donor_sample_counts <- table(donor_block)
    block_requested <- any(donor_sample_counts > 1L)
    block <- if (block_requested) donor_block else NULL
    voom_lm_fit <- .sample_voom_lm_fit_backend()
    fit <- voom_lm_fit(
        dge,
        design = design,
        block = block,
        sample.weights = FALSE,
        normalize.method = "none",
        keep.EList = FALSE
    )
    fit <- limma::eBayes(
        fit,
        trend = FALSE,
        robust = robust
    )
    coefficient <- match("exposure", colnames(fit$coefficients))
    if (is.na(coefficient)) {
        stop("The fitted model does not contain an exposure coefficient.")
    }
    standard_error <- fit$stdev.unscaled[, coefficient] *
        sqrt(fit$s2.post)
    raw_pvalue <- fit$p.value[, coefficient]
    underflow <- is.finite(raw_pvalue) & raw_pvalue == 0
    reported_pvalue <- raw_pvalue
    reported_pvalue[underflow] <- .Machine$double.xmin
    fitted_gene <- rownames(fit$coefficients)
    report_gene <- intersect(requested, fitted_gene)
    report_index <- match(report_gene, fitted_gene)
    result <- data.frame(
        gene = report_gene,
        log2FC = fit$coefficients[report_index, coefficient],
        se = standard_error[report_index],
        statistic = fit$t[report_index, coefficient],
        df_total = fit$df.total[report_index],
        pvalue = reported_pvalue[report_index],
        pvalue_underflow_clamped = underflow[report_index],
        logCPM = fit$Amean[report_index],
        exposure = exposure_name,
        method = if (block_requested) {
            "voomLmFit_donor_duplicateCorrelation_robust_eBayes"
        } else {
            "voomLmFit_independent_samples_robust_eBayes"
        },
        stringsAsFactors = FALSE
    )
    if (anyNA(result$pvalue) || any(!is.finite(result$pvalue)) ||
            any(result$pvalue <= 0 | result$pvalue > 1) ||
            any(!is.finite(result$log2FC)) ||
            any(!is.finite(result$se)) || any(result$se <= 0) ||
            any(!is.finite(result$statistic)) ||
            any(!is.finite(result$logCPM))) {
        stop("voomLmFit returned an invalid effect, standard error, or p-value.")
    }
    correlation <- if (!is.null(block) && !is.null(fit$correlation)) {
        unname(fit$correlation)
    } else {
        NA_real_
    }
    block_preserved <- if (block_requested) {
        !is.null(fit$block) && identical(
            as.character(fit$block),
            as.character(donor_block)
        )
    } else {
        is.null(fit$block)
    }
    if (!block_preserved) {
        stop("voomLmFit did not preserve the requested donor block.")
    }

    list(
        result = result,
        diagnostics = list(
            n_genes_input = nrow(filtered$dge),
            n_genes_donor_support = sum(filtered$support_keep),
            n_genes_total_count = sum(filtered$total_count_keep),
            n_genes_filter_by_expr = if (all(is.na(filtered$edger_keep))) {
                NA_integer_
            } else {
                sum(filtered$edger_keep)
            },
            n_genes_fit = nrow(fit$coefficients),
            n_genes_requested = length(requested),
            n_genes_estimable = sum(requested_estimable),
            n_genes_reported = nrow(result),
            n_p_underflow_clamped = sum(
                result$pvalue_underflow_clamped
            ),
            min_donor_support = filtered$min_donor_support,
            min_cpm = filtered$min_cpm,
            filter_method = filtered$filter_method,
            filter_min_count = filtered$filter_min_count,
            filter_min_total_count = filtered$filter_min_total_count,
            full_library_size = filtered$full_library_size,
            normalisation_factors = setNames(
                dge$samples$norm.factors,
                colnames(dge)
            ),
            donor_support_tested = filtered$donor_support[filtered$keep],
            requested_feature_audit = requested_audit,
            block_requested = block_requested,
            block_preserved = block_preserved,
            n_repeat_donors = sum(donor_sample_counts > 1L),
            max_samples_per_donor = max(donor_sample_counts),
            consensus_correlation = correlation
        )
    )
}

.sample_voom_adjust_by <- function(pvalue, group) {
    adjusted <- rep(NA_real_, length(pvalue))
    for (index in split(seq_along(pvalue), group, drop = TRUE)) {
        adjusted[index] <- stats::p.adjust(pvalue[index], method = "BH")
    }
    adjusted
}

.sample_voom_primary_hypotheses <- function(payload) {
    required_payload <- c("result", "parameters", "donor_table", "provenance")
    missing_payload <- setdiff(required_payload, names(payload))
    if (length(missing_payload)) {
        stop(
            "Canonical primary payload is missing: ",
            paste(missing_payload, collapse = ", ")
        )
    }
    parameters <- payload$parameters
    expected_scalar <- list(
        label = "primary",
        exposure = "smoking_ever",
        celltype_col = "ann_level_3",
        schema_version = "exposomeSC_hlca_primary_edger_v2",
        se_method = "not_available_edgeR_QL",
        statistic_type = "signed_sqrt_qlf"
    )
    for (field in names(expected_scalar)) {
        if (!identical(parameters[[field]], expected_scalar[[field]])) {
            stop("Canonical primary parameter '", field, "' changed.")
        }
    }
    expected_integer <- c(
        min_cells = 20L,
        min_donors = 30L,
        min_group_donors = 10L
    )
    for (field in names(expected_integer)) {
        if (length(parameters[[field]]) != 1L ||
                is.na(parameters[[field]]) ||
                as.integer(parameters[[field]]) != expected_integer[[field]]) {
            stop("Canonical primary parameter '", field, "' changed.")
        }
    }
    fixed_covariates <- c(
        "age_ns1",
        "age_ns2",
        "age_ns3",
        "sex_binary"
    )
    covariates <- parameters$covariates
    if (!is.character(covariates) || anyNA(covariates) ||
            anyDuplicated(covariates) ||
            !identical(utils::head(covariates, 4L), fixed_covariates) ||
            !length(setdiff(covariates, fixed_covariates)) ||
            any(!grepl("^study_", setdiff(covariates, fixed_covariates)))) {
        stop("Canonical primary adjustment covariates changed.")
    }

    result <- as.data.frame(payload$result)
    required_result <- c(
        "gene",
        "gene_symbol",
        "chromosome",
        "sex_chromosome",
        "feature_type",
        "feature_biotype",
        "annotation_mapping_status",
        "celltype",
        "log2FC",
        "se",
        "statistic",
        "pvalue",
        "pvalue_underflow_clamped",
        "padj",
        "logCPM",
        "exposure",
        "n_donors",
        "n_unexposed",
        "n_exposed",
        "min_cells",
        "median_cells",
        "method",
        "padj_global",
        "feature_family",
        "padj_family_celltype",
        "padj_family_global"
    )
    missing_result <- setdiff(required_result, colnames(result))
    if (length(missing_result)) {
        stop(
            "Canonical primary result is missing: ",
            paste(missing_result, collapse = ", ")
        )
    }
    key <- paste(result$gene, result$celltype, sep = "::")
    if (!nrow(result) || anyNA(result$gene) || anyNA(result$celltype) ||
            anyNA(key) || any(!nzchar(result$gene)) ||
            any(!nzchar(result$celltype)) || anyDuplicated(key)) {
        stop("Canonical primary gene-by-cell-type keys are malformed.")
    }
    if (anyNA(result$pvalue) || any(!is.finite(result$pvalue)) ||
            any(result$pvalue <= 0 | result$pvalue > 1) ||
            any(as.character(result$exposure) != "smoking_ever")) {
        stop("Canonical primary p-values or exposure semantics are invalid.")
    }
    if (!is.numeric(result$se) || any(is.nan(result$se)) ||
            !all(is.na(result$se))) {
        stop(
            "Canonical edgeR QL primary SE values must be numeric NA; ",
            "glmQLFTest does not provide coefficient standard errors."
        )
    }
    if (anyNA(result$log2FC) || any(!is.finite(result$log2FC)) ||
            anyNA(result$statistic) || any(!is.finite(result$statistic)) ||
            anyNA(result$logCPM) || any(!is.finite(result$logCPM)) ||
            any(as.character(result$method) != "edgeR_robust_QL")) {
        stop("Canonical primary edgeR QL estimates are malformed.")
    }
    directional <- result$log2FC != 0 & result$statistic != 0
    if (any(sign(result$log2FC[directional]) !=
            sign(result$statistic[directional]))) {
        stop(
            "Canonical primary signed-root QL statistic disagrees with ",
            "the log2 fold-change direction."
        )
    }
    family_levels <- c(
        "host_expression",
        "repertoire_composition",
        "other_annotation"
    )
    if (anyNA(result$feature_family) ||
            any(!result$feature_family %in% family_levels)) {
        stop("Canonical primary feature families are malformed.")
    }
    family_by_gene <- split(result$feature_family, result$gene)
    if (any(vapply(
            family_by_gene,
            function(value) length(unique(value)) != 1L,
            logical(1)
        ))) {
        stop("A primary feature changes family across cell types.")
    }
    expected_padj <- .sample_voom_adjust_by(result$pvalue, result$celltype)
    expected_global <- stats::p.adjust(result$pvalue, method = "BH")
    expected_family_celltype <- .sample_voom_adjust_by(
        result$pvalue,
        interaction(
            result$feature_family,
            result$celltype,
            drop = TRUE,
            lex.order = TRUE
        )
    )
    expected_family_global <- .sample_voom_adjust_by(
        result$pvalue,
        result$feature_family
    )
    observed_adjustments <- list(
        padj = result$padj,
        padj_global = result$padj_global,
        padj_family_celltype = result$padj_family_celltype,
        padj_family_global = result$padj_family_global
    )
    expected_adjustments <- list(
        padj = expected_padj,
        padj_global = expected_global,
        padj_family_celltype = expected_family_celltype,
        padj_family_global = expected_family_global
    )
    adjustment_ok <- vapply(names(observed_adjustments), function(field) {
        isTRUE(all.equal(
            observed_adjustments[[field]],
            expected_adjustments[[field]],
            tolerance = 1e-12,
            check.attributes = FALSE
        ))
    }, logical(1))
    if (!all(adjustment_ok)) {
        stop(
            "Canonical primary BH field(s) are inconsistent: ",
            paste(names(adjustment_ok)[!adjustment_ok], collapse = ", ")
        )
    }

    result_metadata <- S4Vectors::metadata(payload$result)
    engine_parameters <- result_metadata$parameters
    diagnostics <- result_metadata$diagnostics
    expected_engine <- list(
        exposure = "smoking_ever",
        covariates = covariates,
        min_cells = 20L,
        min_donors = 30L,
        min_group_donors = 10L,
        filter_genes = TRUE,
        filter_method = "fixed_support",
        filter_min_cpm = 1,
        filter_min_donors = 10L,
        filter_min_count = 10L,
        filter_min_total_count = 15L,
        robust = TRUE,
        schema_version = "exposomeSC_sc_exwas_edger_v2",
        se_method = "not_available_edgeR_QL",
        statistic_type = "signed_sqrt_qlf"
    )
    for (field in names(expected_engine)) {
        if (!isTRUE(all.equal(
                engine_parameters[[field]],
                expected_engine[[field]],
                tolerance = 0,
                check.attributes = FALSE
            ))) {
            stop("Canonical primary engine parameter '", field, "' changed.")
        }
    }
    celltypes <- sort(unique(as.character(result$celltype)))
    if (!is.list(diagnostics) || !setequal(names(diagnostics), celltypes)) {
        stop("Canonical primary cell types and diagnostics disagree.")
    }
    donor_table <- payload$donor_table
    if (!is.data.frame(donor_table) ||
            !"donor_id" %in% colnames(donor_table) ||
            anyNA(donor_table$donor_id) ||
            anyDuplicated(donor_table$donor_id)) {
        stop("Canonical primary donor table is malformed.")
    }
    donor_ids_by_celltype <- setNames(vector("list", length(celltypes)), celltypes)
    for (celltype in celltypes) {
        cell_result <- result[result$celltype == celltype, , drop = FALSE]
        diagnostic <- diagnostics[[celltype]]
        required_diagnostic <- c(
            "n_donors",
            "donor_ids",
            "n_genes_tested",
            "test_feature_ids",
            "filter_method",
            "filter_min_cpm",
            "filter_min_donors",
            "filter_min_total_count",
            "design_columns",
            "design_rank",
            "residual_df"
        )
        donor_subset_ok <- all(
            diagnostic$donor_ids %in% donor_table$donor_id
        )
        if (!all(required_diagnostic %in% names(diagnostic)) ||
                anyNA(diagnostic$donor_ids) ||
                anyDuplicated(diagnostic$donor_ids) ||
                !donor_subset_ok ||
                diagnostic$n_donors != length(diagnostic$donor_ids) ||
                diagnostic$n_donors < parameters$min_donors ||
                diagnostic$n_genes_tested != nrow(cell_result) ||
                !setequal(
                    as.character(diagnostic$test_feature_ids),
                    as.character(cell_result$gene)
                ) ||
                any(cell_result$n_donors != diagnostic$n_donors) ||
                any(cell_result$n_unexposed <
                    parameters$min_group_donors) ||
                any(cell_result$n_exposed <
                    parameters$min_group_donors) ||
                diagnostic$filter_method != "fixed_support" ||
                diagnostic$filter_min_cpm != 1 ||
                diagnostic$filter_min_donors != 10L ||
                diagnostic$filter_min_total_count != 15L ||
                diagnostic$design_rank != length(diagnostic$design_columns) ||
                diagnostic$residual_df <= 0) {
            stop("Canonical primary diagnostic changed for ", celltype, ".")
        }
        donor_ids_by_celltype[[celltype]] <- as.character(
            diagnostic$donor_ids
        )
    }

    annotation_columns <- c(
        "gene",
        "gene_symbol",
        "chromosome",
        "sex_chromosome",
        "feature_type",
        "feature_biotype",
        "annotation_mapping_status",
        "celltype",
        "feature_family"
    )
    hypotheses <- result[, annotation_columns, drop = FALSE]
    hypotheses$primary_order <- seq_len(nrow(hypotheses))
    hypotheses$primary_log2FC <- result$log2FC
    hypotheses$primary_se <- result$se
    hypotheses$primary_statistic <- result$statistic
    hypotheses$primary_pvalue <- result$pvalue
    hypotheses$primary_padj <- result$padj
    hypotheses$primary_padj_global <- result$padj_global
    hypotheses$primary_padj_family_celltype <-
        result$padj_family_celltype
    hypotheses$primary_padj_family_global <- result$padj_family_global
    hypotheses$primary_n_donors <- result$n_donors
    hypotheses$primary_n_unexposed <- result$n_unexposed
    hypotheses$primary_n_exposed <- result$n_exposed
    hypotheses$primary_min_cells <- result$min_cells
    hypotheses$primary_median_cells <- result$median_cells

    list(
        hypotheses = hypotheses,
        celltypes = celltypes,
        donor_ids_by_celltype = donor_ids_by_celltype,
        diagnostics = diagnostics,
        parameters = parameters,
        engine_parameters = engine_parameters
    )
}

.sample_voom_complete_hypotheses <- function(primary_hypotheses,
                                              fitted = NULL,
                                              reason = NULL,
                                              exposure_name =
                                                  "smoking_ever") {
    required_primary <- c(
        "gene",
        "celltype",
        "feature_family",
        "primary_order"
    )
    missing_primary <- setdiff(
        required_primary,
        colnames(primary_hypotheses)
    )
    if (length(missing_primary) || !nrow(primary_hypotheses) ||
            length(unique(primary_hypotheses$celltype)) != 1L ||
            anyDuplicated(primary_hypotheses$gene)) {
        stop("primary_hypotheses must be one valid primary cell-type family.")
    }
    output <- primary_hypotheses
    n_hypotheses <- nrow(output)
    output$log2FC <- rep(NA_real_, n_hypotheses)
    output$se <- rep(NA_real_, n_hypotheses)
    output$statistic <- rep(NA_real_, n_hypotheses)
    output$df_total <- rep(NA_real_, n_hypotheses)
    output$pvalue <- rep(NA_real_, n_hypotheses)
    output$pvalue_underflow_clamped <- rep(FALSE, n_hypotheses)
    output$logCPM <- rep(NA_real_, n_hypotheses)
    output$exposure <- rep(exposure_name, n_hypotheses)
    output$method <- rep("not_estimable", n_hypotheses)
    output$estimable <- rep(FALSE, n_hypotheses)
    output$not_estimable_reason <- rep(NA_character_, n_hypotheses)

    if (is.null(fitted)) {
        if (!is.character(reason) || length(reason) != 1L ||
                is.na(reason) || !nzchar(reason)) {
            stop("reason is required when an entire family is not estimable.")
        }
        output$not_estimable_reason <- reason
    } else {
        if (!is.list(fitted) ||
                !all(c("result", "diagnostics") %in% names(fitted))) {
            stop("fitted must be a .sample_voom_fit() result.")
        }
        audit <- fitted$diagnostics$requested_feature_audit
        if (!is.data.frame(audit) ||
                !identical(as.character(audit$gene), output$gene) ||
                anyNA(audit$estimable)) {
            stop("Requested-feature audit does not match the primary family.")
        }
        fitted_result <- fitted$result
        if (anyDuplicated(fitted_result$gene) ||
                !setequal(
                    fitted_result$gene,
                    audit$gene[audit$estimable]
                )) {
            stop("Fitted genes and the requested-feature audit disagree.")
        }
        fitted_index <- match(output$gene, fitted_result$gene)
        estimable <- audit$estimable
        result_columns <- c(
            "log2FC",
            "se",
            "statistic",
            "df_total",
            "pvalue",
            "pvalue_underflow_clamped",
            "logCPM",
            "exposure",
            "method"
        )
        missing_result <- setdiff(result_columns, colnames(fitted_result))
        if (length(missing_result)) {
            stop("Fitted result has an incomplete inference schema.")
        }
        for (column in result_columns) {
            output[[column]][estimable] <-
                fitted_result[[column]][fitted_index[estimable]]
        }
        output$estimable <- estimable
        output$not_estimable_reason <- audit$not_estimable_reason
    }
    output$pvalue_for_multiplicity <- ifelse(
        output$estimable,
        output$pvalue,
        1
    )
    output
}

.sample_voom_adjust_family <- function(result) {
    required <- c(
        "gene",
        "celltype",
        "feature_family",
        "primary_order",
        "estimable",
        "pvalue",
        "pvalue_for_multiplicity"
    )
    missing_fields <- setdiff(required, colnames(result))
    if (length(missing_fields)) {
        stop(
            "result is missing: ",
            paste(missing_fields, collapse = ", ")
        )
    }
    key <- paste(result$gene, result$celltype, sep = "::")
    if (anyNA(key) || anyDuplicated(key) ||
            anyNA(result$primary_order) ||
            anyDuplicated(result$primary_order) ||
            !is.logical(result$estimable) || anyNA(result$estimable)) {
        stop("The fixed primary hypothesis keys are malformed.")
    }
    estimated <- result$estimable
    if (anyNA(result$pvalue[estimated]) ||
            any(!is.finite(result$pvalue[estimated])) ||
            any(result$pvalue[estimated] <= 0 |
                result$pvalue[estimated] > 1) ||
            any(!is.na(result$pvalue[!estimated])) ||
            anyNA(result$pvalue_for_multiplicity) ||
            any(!is.finite(result$pvalue_for_multiplicity)) ||
            any(result$pvalue_for_multiplicity <= 0 |
                result$pvalue_for_multiplicity > 1) ||
            any(result$pvalue_for_multiplicity[!estimated] != 1) ||
            any(result$pvalue_for_multiplicity[estimated] !=
                result$pvalue[estimated])) {
        stop("Actual and fixed-family p-value fields are inconsistent.")
    }
    result$padj <- .sample_voom_adjust_by(
        result$pvalue_for_multiplicity,
        result$celltype
    )
    result$padj_global <- stats::p.adjust(
        result$pvalue_for_multiplicity,
        method = "BH"
    )
    result$padj_family_celltype <- .sample_voom_adjust_by(
        result$pvalue_for_multiplicity,
        interaction(
            result$feature_family,
            result$celltype,
            drop = TRUE,
            lex.order = TRUE
        )
    )
    result$padj_family_global <- .sample_voom_adjust_by(
        result$pvalue_for_multiplicity,
        result$feature_family
    )
    result[order(result$primary_order), , drop = FALSE]
}

.sample_voom_validate_code_continuity <- function(
        primary_code_manifest,
        current_code_manifest,
        primary_script_sha256,
        current_primary_script_sha256) {
    required <- c("relative_path", "bytes", "sha256")
    manifests <- list(
        primary = primary_code_manifest,
        current = current_code_manifest
    )
    for (label in names(manifests)) {
        manifest <- manifests[[label]]
        if (!is.data.frame(manifest) ||
                !all(required %in% colnames(manifest)) ||
                anyNA(manifest[, required, drop = FALSE]) ||
                anyDuplicated(manifest$relative_path) ||
                any(!grepl("^[0-9a-f]{64}$", manifest$sha256))) {
            stop(label, " code manifest is malformed.")
        }
    }
    shared_pattern <- "^(DESCRIPTION|NAMESPACE|R/[^/]+[.]R)$"
    primary_shared <- primary_code_manifest[
        grepl(shared_pattern, primary_code_manifest$relative_path),
        required,
        drop = FALSE
    ]
    current_shared <- current_code_manifest[
        grepl(shared_pattern, current_code_manifest$relative_path),
        required,
        drop = FALSE
    ]
    required_loading_files <- c("DESCRIPTION", "NAMESPACE")
    if (!all(required_loading_files %in% primary_shared$relative_path) ||
            !all(required_loading_files %in% current_shared$relative_path) ||
            !any(grepl("^R/[^/]+[.]R$", primary_shared$relative_path)) ||
            !any(grepl("^R/[^/]+[.]R$", current_shared$relative_path))) {
        stop(
            "Primary and current manifests must include DESCRIPTION, ",
            "NAMESPACE, and R package code."
        )
    }
    if (!setequal(
            primary_shared$relative_path,
            current_shared$relative_path
        )) {
        stop("Primary and current shared package-code paths differ.")
    }
    primary_index <- match(
        current_shared$relative_path,
        primary_shared$relative_path
    )
    if (any(current_shared$bytes != primary_shared$bytes[primary_index]) ||
            any(current_shared$sha256 !=
                primary_shared$sha256[primary_index])) {
        stop(
            "DESCRIPTION, NAMESPACE, or R package code drifted after the ",
            "canonical ",
            "primary analysis."
        )
    }
    primary_script_index <- match(
        "analysis/01_hlca_smoking_edger.R",
        primary_code_manifest$relative_path
    )
    if (is.na(primary_script_index) ||
            !identical(
                primary_code_manifest$sha256[[primary_script_index]],
                primary_script_sha256
            ) ||
            !identical(
                primary_script_sha256,
                current_primary_script_sha256
            )) {
        stop(
            "analysis/01_hlca_smoking_edger.R drifted after the canonical ",
            "primary analysis."
        )
    }
    invisible(TRUE)
}

.sample_voom_atomic_promote <- function(temporary, path) {
    backup <- NULL
    if (file.exists(path)) {
        backup <- paste0(path, ".previous-", Sys.getpid())
        if (file.exists(backup)) {
            unlink(backup)
        }
        if (!file.rename(path, backup)) {
            stop("Could not stage the previous output for replacement: ", path)
        }
    }
    if (!file.rename(temporary, path)) {
        if (!is.null(backup) && file.exists(backup)) {
            file.rename(backup, path)
        }
        stop("Could not atomically promote output: ", path)
    }
    if (!is.null(backup) && file.exists(backup)) {
        unlink(backup)
    }
    invisible(path)
}

.sample_voom_atomic_save_rds <- function(object, path, compress = "xz") {
    temporary <- tempfile(
        paste0(".", basename(path), "-"),
        tmpdir = dirname(path)
    )
    on.exit(unlink(temporary), add = TRUE)
    saveRDS(object, temporary, compress = compress)
    .sample_voom_atomic_promote(temporary, path)
}

.sample_voom_atomic_write_csv <- function(data, path) {
    temporary <- tempfile(
        paste0(".", basename(path), "-"),
        tmpdir = dirname(path)
    )
    on.exit(unlink(temporary), add = TRUE)
    utils::write.csv(data, temporary, row.names = FALSE)
    .sample_voom_atomic_promote(temporary, path)
}
