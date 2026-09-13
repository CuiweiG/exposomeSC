# R/bridge.R
# Interoperability bridges: rexposome, Seurat, etc.

#' @include AllClasses.R
#' @include build-scee.R
#' @importFrom S4Vectors DataFrame
#' @importFrom methods is
NULL

#' Convert rexposome ExposomeSet to exposomeSC
#'
#' Bridges the bulk exposome ecosystem with single-cell
#' analysis by extracting the exposure matrix from an
#' \code{ExposomeSet} (rexposome) and combining it with
#' a \code{SingleCellExperiment}.
#'
#' @param exposome_set An \code{ExposomeSet} from rexposome.
#' @param sce A \code{SingleCellExperiment} with matching
#'   donor/sample IDs.
#' @param sample_col Character. Column in \code{colData(sce)}
#'   with sample IDs matching \code{sampleNames(exposome_set)}.
#' @param exposures Character vector (optional). Subset of
#'   exposures to include. Default: all.
#'
#' @return A \code{\linkS4class{SingleCellExposomeExperiment}}.
#'
#' @details
#' This function enables a common workflow: researchers have
#' existing bulk exposome data in rexposome format and want
#' to add single-cell resolution. The exposure matrix and
#' metadata are extracted from the \code{ExposomeSet} and
#' combined with the SCE to create an SCEE.
#'
#' Sample IDs in the ExposomeSet must match entries in
#' \code{colData(sce)[[sample_col]]}.
#'
#' @references
#' Hernandez-Ferrer C et al. (2019). Comprehensive study of
#' the exposome and omic data using rexposome Bioconductor
#' packages. \emph{Bioinformatics} 35:5344-5345.
#'
#' @export
#' @examples
#' # Requires rexposome package
#' # library(rexposome)
#' # es <- loadExposome(...)
#' # scee <- as_scee(es, sce, sample_col = "donor_id")
#' cat("See rexposome vignette for ExposomeSet creation\n")
as_scee <- function(exposome_set, sce, sample_col,
                     exposures = NULL) {
    ## rexposome is not required: an ExposomeSet is a Biobase eSet whose
    ## assayData element "exp" holds exposures (rows) by samples (columns)
    if (!requireNamespace("Biobase", quietly = TRUE))
        stop("Package 'Biobase' required. ",
             "BiocManager::install('Biobase')")

    cls <- class(exposome_set)
    if (!any(grepl("ExposomeSet|ExposomeClust", cls)))
        warning("exposome_set does not appear to be an ",
                "ExposomeSet (class: ", cls[1], ")",
                call. = FALSE)

    stopifnot(is(sce, "SingleCellExperiment"))

    exposure_by_sample <- tryCatch(
        Biobase::assayDataElement(exposome_set, "exp"),
        error = function(e) NULL)
    if (is.null(exposure_by_sample))
        stop("Could not extract exposure data from ExposomeSet: ",
             "no assayData element named 'exp'.")
    exp_mat <- t(as.matrix(exposure_by_sample))

    if (!is.numeric(exp_mat)) {
        converted <- suppressWarnings(matrix(
            as.numeric(exp_mat), nrow = nrow(exp_mat),
            dimnames = dimnames(exp_mat)))
        if (any(is.na(converted) & !is.na(exp_mat)))
            stop("Some exposures are not numeric; encode categorical ",
                 "exposures numerically before conversion.")
        exp_mat <- converted
    }

    if (!is.null(exposures)) {
        exposures <- intersect(exposures, colnames(exp_mat))
        if (length(exposures) == 0)
            stop("No matching exposures found")
        exp_mat <- exp_mat[, exposures, drop = FALSE]
    }

    ## Exposure metadata, aligned with the retained exposures
    exp_info <- tryCatch({
        fi <- Biobase::fData(exposome_set)
        if (all(colnames(exp_mat) %in% rownames(fi))) {
            fi <- fi[colnames(exp_mat), setdiff(colnames(fi), "exposure"),
                     drop = FALSE]
            S4Vectors::DataFrame(exposure = colnames(exp_mat), fi,
                                 check.names = FALSE)
        } else NULL
    }, error = function(e) NULL)

    build_scee(sce, exp_mat, sample_col = sample_col,
        exposure_info = exp_info)
}

#' Convert Seurat metadata to exposure matrix
#'
#' Extracts donor-level metadata from a Seurat object's
#' \code{meta.data} and formats it as an exposure matrix
#' for \code{\link{build_scee}}.
#'
#' @param seurat_meta A data.frame (from
#'   \code{seurat_obj@@meta.data}).
#' @param sample_col Character. Column with donor IDs.
#' @param exposure_cols Character vector. Columns to use as
#'   exposures.
#'
#' @return A numeric matrix (donors x exposures). Each value is the
#'   mean over the donor's cells; a warning names columns whose values
#'   vary within a donor, and non-numeric columns are an error.
#'
#' @export
#' @examples
#' meta <- data.frame(
#'     donor = rep(c("D1","D2","D3"), each = 10),
#'     PM2.5 = rep(c(12, 25, 8), each = 10),
#'     age = rep(c(45, 32, 67), each = 10))
#' exp_mat <- seurat_to_exposure(meta,
#'     sample_col = "donor",
#'     exposure_cols = c("PM2.5", "age"))
#' dim(exp_mat)
seurat_to_exposure <- function(seurat_meta, sample_col,
                                exposure_cols) {
    stopifnot(is.data.frame(seurat_meta))
    stopifnot(sample_col %in% colnames(seurat_meta))

    missing <- setdiff(exposure_cols, colnames(seurat_meta))
    if (length(missing) > 0)
        stop("Columns not found: ",
             paste(missing, collapse = ", "))

    non_numeric <- exposure_cols[!vapply(seurat_meta[exposure_cols],
                                         is.numeric, logical(1))]
    if (length(non_numeric))
        stop("Non-numeric exposure column(s): ",
             paste(non_numeric, collapse = ", "),
             ". Encode them numerically first.")

    donors <- unique(as.character(seurat_meta[[sample_col]]))
    varying <- character()
    exp_mat <- do.call(rbind, lapply(donors, function(d) {
        rows <- seurat_meta[as.character(seurat_meta[[sample_col]]) == d, ,
                            drop = FALSE]
        vapply(exposure_cols, function(col) {
            vals <- rows[[col]]
            observed <- vals[!is.na(vals)]
            if (length(unique(observed)) > 1L)
                varying <<- union(varying, col)
            if (length(observed)) mean(observed) else NA_real_
        }, numeric(1))
    }))
    if (length(varying))
        warning("Values vary within a donor for: ",
                paste(varying, collapse = ", "),
                "; the donor mean is used.", call. = FALSE)
    rownames(exp_mat) <- donors
    colnames(exp_mat) <- exposure_cols
    exp_mat
}
