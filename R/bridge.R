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
    ## Check that the input looks like an ExposomeSet
    ## (rexposome). We don't requireNamespace("rexposome")
    ## because rexposome depends on pryr which is retired
    ## from CRAN. Instead we check for the class and use
    ## Biobase accessors directly.
    if (!requireNamespace("Biobase", quietly = TRUE))
        stop("Package 'Biobase' required. ",
             "BiocManager::install('Biobase')")

    cls <- class(exposome_set)
    if (!any(grepl("ExposomeSet|ExposomeClust", cls)))
        warning("exposome_set does not appear to be an ",
                "ExposomeSet (class: ", cls[1], ")",
                call. = FALSE)

    stopifnot(is(sce, "SingleCellExperiment"))

    if (!requireNamespace("Biobase", quietly = TRUE))
        stop("Package 'Biobase' required. ",
             "BiocManager::install('Biobase')")

    ## Extract exposure matrix from ExposomeSet
    exp_df <- tryCatch(
        Biobase::pData(
            Biobase::assayData(exposome_set)[["exp"]]),
        error = function(e) NULL)
    ## If that fails, try exprs() accessor
    if (is.null(exp_df)) {
        exp_df <- tryCatch(
            as.data.frame(t(Biobase::exprs(exposome_set))),
            error = function(e) NULL)
    }
    if (is.null(exp_df))
        stop("Could not extract exposure data from ",
             "ExposomeSet")

    exp_mat <- as.matrix(exp_df)
    if (!is.numeric(exp_mat))
        exp_mat <- apply(exp_mat, 2, as.numeric)

    if (!is.null(exposures)) {
        exposures <- intersect(exposures, colnames(exp_mat))
        if (length(exposures) == 0)
            stop("No matching exposures found")
        exp_mat <- exp_mat[, exposures, drop = FALSE]
    }

    ## Extract exposure metadata if available
    exp_info <- tryCatch({
        if (requireNamespace("Biobase", quietly = TRUE)) {
            fi <- Biobase::fData(exposome_set)
            S4Vectors::DataFrame(fi)
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
#' @return A numeric matrix (donors x exposures).
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

    donors <- unique(as.character(seurat_meta[[sample_col]]))
    exp_mat <- do.call(rbind, lapply(donors, function(d) {
        rows <- seurat_meta[seurat_meta[[sample_col]] == d, ]
        vapply(exposure_cols, function(col) {
            vals <- rows[[col]]
            if (is.numeric(vals)) mean(vals, na.rm = TRUE)
            else as.numeric(vals[1])
        }, numeric(1))
    }))
    rownames(exp_mat) <- donors
    colnames(exp_mat) <- exposure_cols
    exp_mat
}
