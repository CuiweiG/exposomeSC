# R/utils.R
# Internal utility functions (not exported)

#' @include AllClasses.R
#' @importFrom SummarizedExperiment assay colData
#' @importFrom methods is
#' @importFrom Matrix rowSums
NULL

# Pseudobulk aggregation (internal, not exported)
#
# Sums the raw counts of one cell type within each donor that has at least
# `min_cells` cells of that type. The summation is done by
# scrapper::aggregateAcrossCells(), the successor of the deprecated
# scuttle::aggregateAcrossCells(). Donors keep the order of
# names(table(samples)), so results do not depend on scrapper's own ordering.
# Returns NULL when no donor qualifies.
# @keywords internal
# @noRd
.pseudobulk_aggregate <- function(counts_mat, samples,
                                   cell_types, celltype,
                                   min_cells = 10L) {
    ct_index <- which(!is.na(cell_types) & cell_types == celltype)
    if (!length(ct_index)) return(NULL)

    ct_samples <- samples[ct_index]
    donor_ncells <- table(ct_samples)
    valid <- names(donor_ncells)[donor_ncells >= min_cells]
    if (!length(valid)) return(NULL)

    keep <- ct_samples %in% valid
    aggregated <- scrapper::aggregateAcrossCells(
        counts_mat[, ct_index[keep], drop = FALSE],
        list(donor = factor(ct_samples[keep], levels = valid)),
        compute.sum = TRUE,
        compute.detected = FALSE,
        compute.median = FALSE,
        num.threads = 1L
    )
    donor_order <- match(valid, as.character(aggregated$combinations$donor))
    pb_mat <- aggregated$sums[, donor_order, drop = FALSE]
    dimnames(pb_mat) <- list(rownames(counts_mat), valid)

    list(
        pb_mat = pb_mat,
        valid_donors = valid,
        n_cells = as.integer(donor_ncells[valid])
    )
}

# Log-CPM transformation (internal, not exported)
# @keywords internal
# @noRd
.log_cpm <- function(mat) {
    lib_sizes <- colSums(mat)
    log2(t(t(mat) / lib_sizes * 1e6) + 1)
}

# Construct the backend-independent sc-ExWAS result contract for one cell
# type. Every backend builds its per-cell-type table here and the tables are
# combined by .sc_exwas_finalise().
# @keywords internal
# @noRd
.sc_exwas_result_frame <- function(gene, celltype, log2FC, se,
                                    statistic, pvalue, exposure,
                                    n_donors, method,
                                    logCPM = rep(NA_real_, length(gene)),
                                    n_unexposed = NA_integer_,
                                    n_exposed = NA_integer_,
                                    min_cells = NA_integer_,
                                    median_cells = NA_real_) {
    n <- length(gene)
    vector_arguments <- list(
        log2FC = log2FC,
        se = se,
        statistic = statistic,
        pvalue = pvalue,
        logCPM = logCPM
    )
    wrong_length <- names(vector_arguments)[
        vapply(vector_arguments, length, integer(1)) != n
    ]
    if (length(wrong_length)) {
        stop(
            "Internal sc-ExWAS result vectors have incompatible lengths: ",
            paste(wrong_length, collapse = ", "),
            "."
        )
    }
    if (any(!is.na(pvalue) &
            (!is.finite(pvalue) | pvalue < 0 | pvalue > 1))) {
        stop("Internal sc-ExWAS p-values must be NA or finite in [0, 1].")
    }

    underflow <- !is.na(pvalue) & pvalue == 0
    reported_pvalue <- pvalue
    reported_pvalue[underflow] <- .Machine$double.xmin

    data.frame(
        gene = as.character(gene),
        celltype = rep(as.character(celltype), length.out = n),
        log2FC = as.numeric(log2FC),
        se = as.numeric(se),
        statistic = as.numeric(statistic),
        pvalue = as.numeric(reported_pvalue),
        pvalue_underflow_clamped = as.logical(underflow),
        padj = stats::p.adjust(reported_pvalue, method = "BH"),
        logCPM = as.numeric(logCPM),
        exposure = rep(as.character(exposure), length.out = n),
        n_donors = rep(as.integer(n_donors), length.out = n),
        n_unexposed = rep(as.integer(n_unexposed), length.out = n),
        n_exposed = rep(as.integer(n_exposed), length.out = n),
        min_cells = rep(as.integer(min_cells), length.out = n),
        median_cells = rep(as.numeric(median_cells), length.out = n),
        method = rep(as.character(method), length.out = n),
        stringsAsFactors = FALSE
    )
}

# Columns of every run_sc_exwas() result, in order.
.sc_exwas_columns <- c(
    "gene", "celltype", "log2FC", "se", "statistic", "pvalue",
    "pvalue_underflow_clamped", "padj", "padj_global", "logCPM",
    "exposure", "n_donors", "n_unexposed", "n_exposed", "min_cells",
    "median_cells", "method"
)

# A zero-row result with the full column set and types.
.sc_exwas_empty_frame <- function() {
    S4Vectors::DataFrame(
        gene = character(),
        celltype = character(),
        log2FC = numeric(),
        se = numeric(),
        statistic = numeric(),
        pvalue = numeric(),
        pvalue_underflow_clamped = logical(),
        padj = numeric(),
        padj_global = numeric(),
        logCPM = numeric(),
        exposure = character(),
        n_donors = integer(),
        n_unexposed = integer(),
        n_exposed = integer(),
        min_cells = integer(),
        median_cells = numeric(),
        method = character()
    )
}

# Combine per-cell-type tables from any backend into the run_sc_exwas()
# result: BH across all gene-by-cell-type tests (padj_global), rows ordered
# by padj_global then pvalue, the same columns for every backend, and the
# fitting diagnostics and parameters in metadata().
.sc_exwas_finalise <- function(tables, diagnostics, parameters) {
    tables <- Filter(Negate(is.null), tables)
    if (!length(tables)) {
        result <- .sc_exwas_empty_frame()
    } else {
        output <- do.call(rbind, unname(tables))
        output$padj_global <- stats::p.adjust(output$pvalue, method = "BH")
        output <- output[order(output$padj_global, output$pvalue),
                         .sc_exwas_columns, drop = FALSE]
        rownames(output) <- NULL
        result <- S4Vectors::DataFrame(output, check.names = FALSE)
    }
    S4Vectors::metadata(result)$diagnostics <- diagnostics
    S4Vectors::metadata(result)$parameters <- parameters
    result
}

# Inputs shared by the run_sc_exwas() backends: checks that the exposure,
# covariates and cell-type column exist, and returns the count matrix with
# each cell's donor and cell type.
.sc_exwas_inputs <- function(x, exposure, covariates, celltype_col,
                             celltypes) {
    exp_data <- exposureData(x)
    missing_fields <- setdiff(c(exposure, covariates), colnames(exp_data))
    if (length(missing_fields)) {
        stop(
            "Exposure/covariate column(s) not found in exposureData: ",
            paste(missing_fields, collapse = ", ")
        )
    }
    if ("exposure" %in% covariates) {
        stop(
            "A covariate cannot be called 'exposure': that name is used for ",
            "the exposure term of the model. Rename the variable in ",
            "exposureData(x)."
        )
    }
    if (exposure %in% covariates) {
        stop("The exposure must not also be listed as a covariate.")
    }
    .check_formula_names(covariates, "covariate")
    cd <- SummarizedExperiment::colData(x)
    if (!celltype_col %in% colnames(cd)) {
        stop("celltype_col '", celltype_col, "' is not a column of colData(x).")
    }
    cell_types <- as.character(cd[[celltype_col]])
    if (is.null(celltypes)) {
        celltypes <- sort(unique(cell_types))
    }
    list(
        exp_data = exp_data,
        counts = SummarizedExperiment::assay(x, "counts"),
        samples = .cell_samples(x),
        cell_types = cell_types,
        celltypes = celltypes
    )
}

# Donor-level design data (exposure first, then covariates) for the donors
# of one pseudobulk, restricted to donors with complete values. Also returns
# the group sizes when the exposure takes exactly two values.
.sc_exwas_design_data <- function(exp_data, donors, exposure, covariates) {
    design_data <- data.frame(
        exposure = as.numeric(exp_data[donors, exposure]),
        row.names = donors,
        check.names = FALSE
    )
    for (covariate in covariates) {
        design_data[[covariate]] <- as.numeric(exp_data[donors, covariate])
    }
    complete <- stats::complete.cases(design_data)
    design_data <- design_data[complete, , drop = FALSE]
    exposure_levels <- sort(unique(design_data$exposure))
    group_counts <- if (length(exposure_levels) == 2L) {
        table(factor(design_data$exposure, levels = exposure_levels))
    } else {
        integer()
    }
    list(
        data = design_data,
        complete = complete,
        n_unexposed = if (length(group_counts)) {
            unname(group_counts[[1L]])
        } else {
            NA_integer_
        },
        n_exposed = if (length(group_counts)) {
            unname(group_counts[[2L]])
        } else {
            NA_integer_
        },
        group_counts = group_counts
    )
}

# Seed the random number generator for the rest of the calling function, then
# put the caller's generator back exactly as it was when that function exits.
#
# `withr::local_seed()` does the seeding. On its own it restores `.Random.seed`
# but, when the caller had no seed yet, leaves any generator kind it set in
# place, so a later draw by the user would come from L'Ecuyer-CMRG instead of
# Mersenne-Twister. The kind and seed are therefore recorded here and restored
# after withr's own clean-up, in the same order the package used before.
#
# `envir` is the frame whose exit ends the scope. A helper whose caller keeps
# drawing from the seeded stream passes its caller's frame.
# @keywords internal
# @noRd
.local_rng_scope <- function(seed, kind = NULL, envir = parent.frame()) {
    had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
    old_seed <- if (had_seed) get(".Random.seed", envir = globalenv())
    old_kind <- RNGkind()
    withr::defer({
        do.call(RNGkind, as.list(old_kind))
        if (had_seed) {
            assign(".Random.seed", old_seed, envir = globalenv())
        } else if (exists(".Random.seed", envir = globalenv(),
                          inherits = FALSE)) {
            rm(".Random.seed", envir = globalenv())
        }
    }, envir = envir)
    withr::local_seed(seed, .local_envir = envir, .rng_kind = kind)
    invisible(NULL)
}

# Experimental interfaces warn once per session. Tests reset the record with
# .reset_experimental_warnings() before asserting the warning.
.experimental_warnings <- new.env(parent = emptyenv())

.warn_experimental <- function(key, message) {
    if (!isTRUE(.experimental_warnings[[key]])) {
        assign(key, TRUE, envir = .experimental_warnings)
        warning(message, call. = FALSE)
    }
    invisible(NULL)
}

.reset_experimental_warnings <- function() {
    rm(list = ls(.experimental_warnings, all.names = TRUE),
       envir = .experimental_warnings)
    invisible(NULL)
}

# A misspelled argument absorbed by ... must not silently change an analysis.
.stop_on_unused_dots <- function(function_name, ...) {
    n_extra <- ...length()
    if (!n_extra) return(invisible(NULL))
    extra_names <- ...names()
    if (is.null(extra_names)) extra_names <- rep("", n_extra)
    extra_names[is.na(extra_names) | !nzchar(extra_names)] <- "<unnamed>"
    stop("Unused argument(s) in ", function_name, "(): ",
         paste(extra_names, collapse = ", "), ".", call. = FALSE)
}

# Requested genes that are not in the data are left out with a warning that
# names them, so that the tested family does not shrink without a word.
# @keywords internal
# @noRd
.match_target_genes <- function(target_genes, available) {
    absent <- setdiff(target_genes, available)
    if (length(absent)) {
        shown <- paste(utils::head(absent, 10L), collapse = ", ")
        more <- if (length(absent) > 10L)
            sprintf(" and %d more", length(absent) - 10L) else ""
        warning("target_genes not found in the data are not tested: ",
                shown, more, ".", call. = FALSE)
    }
    intersect(target_genes, available)
}

# Arguments that name one cell type or one exposure. A longer vector would be
# recycled in comparisons such as cell_types == celltype and analyse several
# cell types as one. A factor or a number is taken as its label: the columns
# are compared as text, and a factor would index exposureData by its code.
# @keywords internal
# @noRd
.single_name <- function(value, argument) {
    if (is.factor(value) || is.numeric(value)) value <- as.character(value)
    if (!is.character(value) || length(value) != 1L || is.na(value) ||
            !nzchar(value)) {
        stop("'", argument, "' must be a single name.", call. = FALSE)
    }
    value
}

# The exposures of a mixture must be distinct, and none may also be a
# covariate: the covariate would replace the exposure's quantile score with
# its raw values, or the adjustment would be dropped without a word.
# @keywords internal
# @noRd
.check_mixture_exposures <- function(exposures, covariates) {
    if (anyDuplicated(exposures)) {
        stop("'exposures' must not contain duplicated names.", call. = FALSE)
    }
    overlap <- intersect(exposures, covariates)
    if (length(overlap)) {
        overlap_names <- paste(overlap, collapse = ", ")
        stop("An exposure must not also be listed as a covariate: ",
            overlap_names, ".", call. = FALSE)
    }
    invisible(TRUE)
}

# Exposure, covariate and moderator names are pasted into model formulas, so
# each must be a syntactic R name: 'age-years' would be read as age minus
# years, and a ratio such as 'MEHP/DEHP' as a nested term. Three kinds of
# name that make.names() leaves as they are cannot be used either: in a
# formula '.' stands for every column of the model data other than the
# response, and R reads '...' and '..1', '..2' and so on as arguments passed
# down from a calling function.
# @keywords internal
# @noRd
.check_formula_names <- function(names, role) {
    names <- as.character(names)
    bad <- names[make.names(names) != names]
    if (length(bad)) {
        bad_names <- paste0("'", bad, "'", collapse = ", ")
        suggested <- paste0("'", make.names(bad), "'", collapse = ", ")
        stop("Non-syntactic ", role, " name(s) cannot be used in a model ",
            "formula: ", bad_names, ". Rename the variable(s) in ",
            "exposureData(x), for example to ", suggested, ".",
            call. = FALSE)
    }
    special <- grep("^(\\.|\\.\\.\\.|\\.\\.[0-9]+)$", names, value = TRUE)
    if (length(special)) {
        stop("The ", role, " name(s) ",
            paste0("'", special, "'", collapse = ", "), " cannot be used in ",
            "a model formula, where R reads '.', '...' and names such as ",
            "'..1' specially. Rename the variable(s) in exposureData(x).",
            call. = FALSE)
    }
    invisible(TRUE)
}

# The analyses put the user's variables next to variables of their own (the
# response, the exposure term, an offset, a grouping factor, the columns of a
# meta-regression or a partial correlation), and some packages they call add
# columns to the model data. A user's variable with one of those names
# replaced the other variable, was replaced by it or was dropped, mostly
# without an error that named the clash, so such names are refused. `role`
# says what kind of variable `names` holds and `location` where it can be
# renamed.
# @keywords internal
# @noRd
.check_reserved_names <- function(names, role, reserved,
                                  location = "exposureData(x)") {
    clash <- intersect(names, reserved)
    if (length(clash)) {
        article <- if (grepl("^[aeiou]", role)) "An " else "A "
        clash_names <- paste0("'", clash, "'", collapse = ", ")
        used_for <- if (length(reserved) > 1L) {
            paste("the names", paste(reserved, collapse = ", "),
                "for variables")
        } else {
            paste("the name", reserved, "for a variable")
        }
        stop(article, role, " cannot be called ", clash_names, ": the ",
            "analysis uses ", used_for, " of its own. Rename the variable ",
            "in ", location, ".", call. = FALSE)
    }
    invisible(TRUE)
}
