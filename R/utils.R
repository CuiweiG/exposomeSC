# R/utils.R
# Internal utility functions (not exported)

#' @include AllClasses.R
#' @importFrom SummarizedExperiment assay colData
#' @importFrom methods is
#' @importFrom Matrix rowSums
NULL

# Pseudobulk aggregation (internal, not exported)
# @keywords internal
# @noRd
.pseudobulk_aggregate <- function(counts_mat, samples,
                                   cell_types, celltype,
                                   min_cells = 10L) {
    ct_mask <- !is.na(cell_types) & cell_types == celltype
    if (!any(ct_mask)) return(NULL)

    ct_index <- which(ct_mask)
    ct_samples <- samples[ct_index]
    donor_ncells <- table(ct_samples)
    valid <- names(donor_ncells[donor_ncells >= min_cells])
    if (length(valid) == 0) return(NULL)

    keep <- ct_samples %in% valid
    selected_index <- ct_index[keep]
    donor_index <- match(ct_samples[keep], valid)
    grouping <- Matrix::sparseMatrix(
        i = seq_along(donor_index),
        j = donor_index,
        x = 1,
        dims = c(length(donor_index), length(valid)),
        dimnames = list(NULL, valid)
    )
    ## A single sparse matrix multiplication replaces one full gene-matrix
    ## slice per donor. The resulting dense pseudobulk matrix is small
    ## (features x donors) and is the representation expected by edgeR/DESeq2.
    pb_mat <- as.matrix(
        counts_mat[, selected_index, drop = FALSE] %*% grouping
    )
    rownames(pb_mat) <- rownames(counts_mat)
    colnames(pb_mat) <- valid

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

# Construct the backend-independent sc-ExWAS result contract.
# @keywords internal
# @noRd
.sc_exwas_result_frame <- function(gene, celltype, log2FC, se,
                                    statistic, pvalue, exposure,
                                    n_donors, method,
                                    n_unexposed = NA_integer_,
                                    n_exposed = NA_integer_,
                                    min_cells = NA_integer_,
                                    median_cells = NA_real_) {
    n <- length(gene)
    vector_arguments <- list(
        log2FC = log2FC,
        se = se,
        statistic = statistic,
        pvalue = pvalue
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

# Resolve current sc-ExWAS fields while accepting explicit legacy aliases.
# @keywords internal
# @noRd
.resolve_exwas_result_column <- function(data, candidates, role) {
    available <- candidates[candidates %in% colnames(data)]
    if (!length(available)) {
        stop(
            "No supported ", role, " column found. Tried: ",
            paste(candidates, collapse = ", "),
            "."
        )
    }
    available[[1L]]
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
