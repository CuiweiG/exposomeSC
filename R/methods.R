# R/methods.R
#' @include AllClasses.R
#' @include AllGenerics.R
#' @importFrom methods setMethod show slot validObject new
#' @importFrom utils head
NULL

#' @rdname SCEE-accessors
#' @export
setMethod("exposureData", "SingleCellExposomeExperiment",
    function(x) slot(x, "exposureData"))

#' @rdname SCEE-accessors
#' @export
setMethod("exposureInfo", "SingleCellExposomeExperiment",
    function(x) slot(x, "exposureInfo"))

#' @rdname SCEE-accessors
#' @export
setMethod("cellSampleMap", "SingleCellExposomeExperiment",
    function(x) slot(x, "sampleMap"))

#' @rdname SCEE-accessors
#' @export
setMethod("exposureVariables", "SingleCellExposomeExperiment",
    function(x) colnames(slot(x, "exposureData")))

#' @rdname SCEE-accessors
#' @export
setMethod("exposureData<-",
    "SingleCellExposomeExperiment",
    function(x, value) {
        if (!is.matrix(value) || !is.numeric(value)) {
            stop("replacement exposureData must be a numeric matrix.")
        }
        slot(x, "exposureData") <- value
        validObject(x)
        x
    })

#' @rdname SCEE-accessors
#' @export
setMethod("exposureInfo<-",
    "SingleCellExposomeExperiment",
    function(x, value) {
        stopifnot(is(value, "DataFrame"))
        slot(x, "exposureInfo") <- value
        validObject(x)
        x
    })

#' @rdname SingleCellExposomeExperiment-class
#' @aliases [,SingleCellExposomeExperiment,ANY,ANY,ANY-method
#' @param x A SingleCellExposomeExperiment.
#' @param i Row (gene) subscript.
#' @param j Column (cell) subscript.
#' @param ... Additional arguments.
#' @param drop Logical. Ignored (kept for S4 compatibility).
#' @export
setMethod("[", "SingleCellExposomeExperiment",
    function(x, i, j, ..., drop = TRUE) {
        sce_sub <- callNextMethod()
        sm <- slot(x, "sampleMap")
        ed <- slot(x, "exposureData")
        ei <- slot(x, "exposureInfo")

        if (!missing(j)) {
            ## Update sampleMap for retained cells
            retained <- colnames(sce_sub)
            if (anyDuplicated(retained)) {
                stop(
                    "Duplicated cell selection is not supported for ",
                    "SingleCellExposomeExperiment."
                )
            }
            sm_idx <- match(retained, as.character(sm$cell_id))
            if (anyNA(sm_idx)) {
                stop("Retained cells are missing from sampleMap.")
            }
            sm <- sm[sm_idx, , drop = FALSE]

            ## Subset exposureData to retained donors only
            retained_donors <- unique(as.character(
                sm$sample_id))
            if (!is.null(rownames(ed))) {
                if (length(retained_donors) > 0) {
                    keep <- match(retained_donors, rownames(ed))
                    if (anyNA(keep)) {
                        stop("Retained sample IDs are missing from exposureData.")
                    }
                    ed <- ed[keep, , drop = FALSE]
                } else {
                    ed <- ed[FALSE, , drop = FALSE]
                }
            }
        }

        new("SingleCellExposomeExperiment",
            sce_sub,
            exposureData = ed,
            exposureInfo = ei,
            sampleMap = sm)
    })

#' @rdname SingleCellExposomeExperiment-class
#' @param object A SingleCellExposomeExperiment.
#' @export
setMethod("show", "SingleCellExposomeExperiment",
    function(object) {
    callNextMethod()
    ed <- slot(object, "exposureData")
    message("exposureData(", nrow(ed), " samples x ",
        ncol(ed), " exposures): ",
        paste(utils::head(colnames(ed), 5), collapse = ", "),
        if (ncol(ed) > 5) " ..." else "")
    sm <- slot(object, "sampleMap")
    message("sampleMap: ", nrow(sm), " cell-sample links")
})
