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
    function(x) {
        S4Vectors::DataFrame(
            cell_id = colnames(x),
            sample_id = .cell_samples(x)
        )
    })

#' @rdname SCEE-accessors
#' @export
setMethod("sampleColumn", "SingleCellExposomeExperiment",
    function(x) slot(x, "sampleColumn"))

#' @rdname SCEE-accessors
#' @export
setReplaceMethod("sampleColumn", "SingleCellExposomeExperiment",
    function(x, value) {
        if (!is.character(value) || length(value) != 1L || is.na(value) ||
                !nzchar(value)) {
            stop("the sample column must be a single column name.")
        }
        slot(x, "sampleColumn") <- value
        validObject(x)
        x
    })

## Sample (donor) ID of every cell, read from colData
.cell_samples <- function(x) {
    as.character(SummarizedExperiment::colData(x)[[slot(x, "sampleColumn")]])
}

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
        ed <- slot(x, "exposureData")

        if (!missing(j)) {
            if (anyDuplicated(colnames(sce_sub))) {
                stop(
                    "Duplicated cell selection is not supported for ",
                    "SingleCellExposomeExperiment."
                )
            }
            ## Keep the exposure rows of the donors that still have cells
            retained_donors <- unique(.cell_samples(sce_sub))
            keep <- match(retained_donors, rownames(ed))
            if (anyNA(keep)) {
                stop("Retained sample IDs are missing from exposureData.")
            }
            ed <- ed[keep, , drop = FALSE]
        }

        new("SingleCellExposomeExperiment",
            sce_sub,
            exposureData = ed)
    })

#' @rdname SingleCellExposomeExperiment-class
#' @aliases cbind,SingleCellExposomeExperiment-method
#' @param deparse.level Ignored.
#' @importFrom BiocGenerics cbind
#' @export
setMethod("cbind", "SingleCellExposomeExperiment",
    function(..., deparse.level = 1) {
        objects <- list(...)
        if (!all(vapply(objects, is, logical(1),
                "SingleCellExposomeExperiment"))) {
            stop("cbind() combines SingleCellExposomeExperiment objects ",
                 "only; build the others with build_scee() first.")
        }
        first <- objects[[1L]]
        for (other in objects[-1L]) {
            if (!identical(slot(other, "sampleColumn"),
                    slot(first, "sampleColumn"))) {
                stop("All objects must use the same sample column.")
            }
            if (!identical(colnames(slot(other, "exposureData")),
                    colnames(slot(first, "exposureData")))) {
                stop("All objects must have the same exposure variables ",
                     "in the same order.")
            }
            if (!isTRUE(all.equal(as.data.frame(slot(other, "exposureInfo")),
                    as.data.frame(slot(first, "exposureInfo")),
                    check.attributes = FALSE))) {
                stop("All objects must have the same exposureInfo.")
            }
        }

        ## Stack the exposure rows; a donor in several objects must have the
        ## same values in each
        ed <- do.call(rbind, lapply(objects, slot, "exposureData"))
        donors <- rownames(ed)
        repeated <- unique(donors[duplicated(donors)])
        for (donor in repeated) {
            rows <- ed[donors == donor, , drop = FALSE]
            same <- apply(rows, 2L, function(v) {
                all(is.na(v)) || (!anyNA(v) && all(v == v[1L]))
            })
            if (!all(same)) {
                stop("Sample '", donor, "' has different exposure values ",
                     "in the objects being combined.")
            }
        }
        ed <- ed[!duplicated(donors), , drop = FALSE]

        sces <- lapply(objects, methods::as, "SingleCellExperiment")
        combined <- do.call(cbind, sces)
        new("SingleCellExposomeExperiment",
            combined,
            exposureData = ed,
            exposureInfo = slot(first, "exposureInfo"),
            sampleColumn = slot(first, "sampleColumn"))
    })

## colData<- (also reached through $<- and [[<-) and dimnames<- (reached
## through colnames<- and rownames<-) check validity, so that an edit cannot
## leave the donor column out of step with exposureData.

#' @rdname SingleCellExposomeExperiment-class
#' @aliases colData<-,SingleCellExposomeExperiment,DataFrame-method
#' @param value The replacement \code{colData} or \code{dimnames}.
#' @importFrom SummarizedExperiment colData<-
#' @export
setReplaceMethod("colData", c("SingleCellExposomeExperiment", "DataFrame"),
    function(x, ..., value) {
        x <- callNextMethod()
        validObject(x)
        x
    })

#' @rdname SingleCellExposomeExperiment-class
#' @aliases colData<-,SingleCellExposomeExperiment,NULL-method
#' @export
setReplaceMethod("colData", c("SingleCellExposomeExperiment", "NULL"),
    function(x, ..., value) {
        x <- callNextMethod()
        validObject(x)
        x
    })

#' @rdname SingleCellExposomeExperiment-class
#' @aliases dimnames<-,SingleCellExposomeExperiment,list-method
#' @export
setReplaceMethod("dimnames", c("SingleCellExposomeExperiment", "list"),
    function(x, value) {
        x <- callNextMethod()
        validObject(x)
        x
    })

#' @rdname SingleCellExposomeExperiment-class
#' @aliases dimnames<-,SingleCellExposomeExperiment,NULL-method
#' @export
setReplaceMethod("dimnames", c("SingleCellExposomeExperiment", "NULL"),
    function(x, value) {
        x <- callNextMethod()
        validObject(x)
        x
    })

#' @rdname SingleCellExposomeExperiment-class
#' @param object A SingleCellExposomeExperiment.
#' @export
setMethod("show", "SingleCellExposomeExperiment",
    function(object) {
    callNextMethod()
    ed <- slot(object, "exposureData")
    cat("exposureData(", nrow(ed), " samples x ",
        ncol(ed), " exposures): ",
        paste(utils::head(colnames(ed), 5), collapse = ", "),
        if (ncol(ed) > 5) " ..." else "", "\n", sep = "")
    cat("sampleColumn: ", slot(object, "sampleColumn"), "\n", sep = "")
})
