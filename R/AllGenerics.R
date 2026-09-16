# R/AllGenerics.R

#' @include AllClasses.R
#' @importFrom methods setGeneric
NULL

#' Access exposure data matrix
#' @param x A \code{\linkS4class{SingleCellExposomeExperiment}}.
#' @return Numeric matrix (samples x exposures).
#' @export
#' @rdname SCEE-accessors
#' @examples
#' library(SingleCellExperiment)
#' counts <- matrix(1L, nrow = 5, ncol = 12,
#'     dimnames = list(paste0("G", 1:5), paste0("c", 1:12)))
#' sce <- SingleCellExperiment(assays = list(counts = counts),
#'     colData = S4Vectors::DataFrame(
#'         cell_id = paste0("c", 1:12),
#'         donor_id = rep(paste0("D", 1:4), each = 3)))
#' exp_mat <- matrix(c(12, 25, 8, 30, 45, 32, 67, 51), nrow = 4,
#'     dimnames = list(paste0("D", 1:4), c("PM2.5", "age")))
#' scee <- build_scee(sce, exp_mat, sample_col = "donor_id")
#'
#' exposureVariables(scee)
#' exposureData(scee)
#' exposureInfo(scee)
#' head(cellSampleMap(scee))
#'
#' ## Replace the exposure matrix, here rescaling PM2.5
#' exp_new <- exposureData(scee)
#' exp_new[, "PM2.5"] <- exp_new[, "PM2.5"] / 10
#' exposureData(scee) <- exp_new
#' exposureData(scee)[, "PM2.5"]
setGeneric("exposureData", function(x)
    standardGeneric("exposureData"))

#' Access exposure metadata
#' @param x A \code{\linkS4class{SingleCellExposomeExperiment}}.
#' @return A \code{DataFrame}.
#' @export
#' @rdname SCEE-accessors
setGeneric("exposureInfo", function(x)
    standardGeneric("exposureInfo"))

#' Access sample-cell mapping
#' @param x A \code{\linkS4class{SingleCellExposomeExperiment}}.
#' @return A \code{DataFrame}.
#' @export
#' @rdname SCEE-accessors
setGeneric("cellSampleMap", function(x)
    standardGeneric("cellSampleMap"))

#' Get exposure variable names
#' @param x A \code{\linkS4class{SingleCellExposomeExperiment}}.
#' @return Character vector.
#' @export
#' @rdname SCEE-accessors
setGeneric("exposureVariables", function(x)
    standardGeneric("exposureVariables"))

#' Replace exposure data matrix
#' @param x A \code{\linkS4class{SingleCellExposomeExperiment}}.
#' @param value Numeric matrix with same row names (samples).
#' @export
#' @rdname SCEE-accessors
setGeneric("exposureData<-", function(x, value)
    standardGeneric("exposureData<-"))

#' Replace exposure metadata
#' @param x A \code{\linkS4class{SingleCellExposomeExperiment}}.
#' @param value A \code{DataFrame}.
#' @export
#' @rdname SCEE-accessors
setGeneric("exposureInfo<-", function(x, value)
    standardGeneric("exposureInfo<-"))

#' Run cell-type-specific ExWAS
#' @param x A \code{\linkS4class{SingleCellExposomeExperiment}}.
#' @param ... Additional arguments.
#' @return A \code{DataFrame} of results.
#' @export
#' @rdname run_sc_exwas
setGeneric("run_sc_exwas", function(x, ...)
    standardGeneric("run_sc_exwas"))
