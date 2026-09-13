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
#' showClass("SingleCellExposomeExperiment")
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
setGeneric("sampleMap", function(x)
    standardGeneric("sampleMap"))

#' Get exposure variable names
#' @param x A \code{\linkS4class{SingleCellExposomeExperiment}}.
#' @return Character vector.
#' @export
#' @rdname SCEE-accessors
setGeneric("exposureNames", function(x)
    standardGeneric("exposureNames"))

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
