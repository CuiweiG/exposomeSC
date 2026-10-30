# R/AllClasses.R
# Core S4 classes: SingleCellExposomeExperiment + Network classes

#' @import methods
#' @importFrom S4Vectors DataFrame
#' @importClassesFrom S4Vectors DataFrame
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
NULL

#' SingleCellExposomeExperiment: Container for single-cell
#' exposome data
#'
#' Extends \code{\link[SingleCellExperiment]{SingleCellExperiment}} with
#' sample-level environmental exposure measurements. This enables
#' cell-type-specific exposome-wide association studies (sc-ExWAS).
#'
#' The donor of each cell is read from one column of \code{colData}, whose
#' name is stored in the object (see \code{\link{sampleColumn}}). The
#' cell-to-donor mapping is not stored a second time, so renaming cells,
#' combining objects and editing \code{colData} cannot leave two copies out
#' of step. The object is valid when every donor in that column has exactly
#' one row in \code{exposureData} and every row of \code{exposureData}
#' belongs to at least one cell.
#'
#' Replacing \code{colData}, or a column of it with \code{$<-} or
#' \code{[[<-}, checks validity again, so an edit that removes the donor
#' column or introduces a donor without exposure data is an error; so does
#' replacing the cell names with \code{colnames<-}. To change donor
#' identifiers, edit \code{exposureData} and \code{colData} together and
#' rebuild the object with \code{\link{build_scee}}.
#'
#' @slot exposureData Numeric matrix. Rows = samples (donors),
#'   columns = exposure variables.
#' @slot exposureInfo \code{DataFrame}. Metadata for each
#'   exposure (family, unit, LOD).
#' @slot sampleColumn Character scalar. Name of the \code{colData} column
#'   holding each cell's sample (donor) identifier.
#' @return The subsetting method returns a valid
#'   \code{SingleCellExposomeExperiment}. When cells are subset, its
#'   \code{exposureData} is restricted to the donors that keep at least one
#'   cell; feature subsetting leaves the donor-level slots unchanged.
#'   \code{cbind} combines objects with the same exposure variables, exposure
#'   metadata and sample column; donors present in more than one object must
#'   have identical exposure values. The replacement methods return the
#'   modified object. The \code{show} method is called for its display side
#'   effect and returns \code{NULL} invisibly.
#'
#' @section Data hierarchy:
#' \preformatted{
#' Donor (n=50-1000)
#'   |-- Exposures (PM2.5, Pb, BPA, ...)   <- sample level
#'   +-- Cells x Genes matrix
#'       |-- T cells (1000 cells)
#'       |-- NK cells (500 cells)            <- cell level
#'       +-- Monocytes (800 cells)
#' }
#'
#' @references
#' van der Wijst MGP et al. (2020). The single-cell eQTLGen
#' consortium. \emph{eLife} 9:e52155. \doi{10.7554/eLife.52155}
#'
#' @name SingleCellExposomeExperiment-class
#' @exportClass SingleCellExposomeExperiment
#' @examples
#' showClass("SingleCellExposomeExperiment")
.SCEE <- setClass("SingleCellExposomeExperiment",
    contains = "SingleCellExperiment",
    slots = list(
        exposureData   = "matrix",
        exposureInfo   = "DataFrame",
        sampleColumn   = "character"
    ),
    prototype = list(
        exposureData   = matrix(nrow = 0, ncol = 0),
        exposureInfo   = S4Vectors::DataFrame(),
        sampleColumn   = character()
    )
)

setValidity("SingleCellExposomeExperiment", function(object) {
    msg <- character()
    ed <- slot(object, "exposureData")
    ei <- slot(object, "exposureInfo")
    sc <- slot(object, "sampleColumn")

    if (!is.numeric(ed)) {
        msg <- c(msg, "exposureData must be a numeric matrix")
    } else {
        if (any(is.infinite(ed)) || any(is.nan(ed))) {
            msg <- c(
                msg,
                "exposureData may contain NA, but not NaN or infinite values"
            )
        }
    }

    if (nrow(ed) > 0) {
        if (is.null(rownames(ed)) || anyNA(rownames(ed)) ||
                any(!nzchar(rownames(ed)))) {
            msg <- c(msg, "exposureData must have non-missing sample row names")
        } else if (anyDuplicated(rownames(ed))) {
            msg <- c(msg, "exposureData sample row names must be unique")
        }
        if (is.null(colnames(ed)) || anyNA(colnames(ed)) ||
                any(!nzchar(colnames(ed)))) {
            msg <- c(msg, "exposureData must have non-missing column names")
        } else if (anyDuplicated(colnames(ed))) {
            msg <- c(msg, "exposureData column names must be unique")
        }
    }

    if (ncol(ed) > 0) {
        if (nrow(ei) != ncol(ed)) {
            msg <- c(msg, paste0(
                "nrow(exposureInfo) [", nrow(ei),
                "] must equal ncol(exposureData) [",
                ncol(ed), "]"))
        }
        if (!"exposure" %in% colnames(ei)) {
            msg <- c(msg, "exposureInfo must contain an exposure column")
        } else if (!identical(
                as.character(ei$exposure),
                colnames(ed)
            )) {
            msg <- c(msg, paste0(
                "exposureInfo$exposure must match ",
                "colnames(exposureData) in order"))
        }
    }

    ## An empty object (no cells and no exposure rows) needs no sample column
    if (ncol(object) == 0L && nrow(ed) == 0L && length(sc) == 0L) {
        return(if (length(msg) == 0L) TRUE else msg)
    }

    cells <- colnames(object)
    if (ncol(object) > 0L && (is.null(cells) || anyNA(cells) ||
            any(!nzchar(cells)))) {
        msg <- c(msg, "cells (colnames) must have non-missing names")
    } else if (anyDuplicated(cells)) {
        msg <- c(msg, "cell names (colnames) must be unique")
    }

    cd <- SummarizedExperiment::colData(object)
    if (length(sc) != 1L || is.na(sc) || !nzchar(sc)) {
        msg <- c(msg, "sampleColumn must be a single column name")
    } else if (!sc %in% colnames(cd)) {
        msg <- c(msg, paste0(
            "sampleColumn '", sc, "' is not a column of colData"))
    } else {
        donors <- as.character(cd[[sc]])
        if (anyNA(donors) || any(!nzchar(donors))) {
            msg <- c(msg, paste0(
                "colData column '", sc, "' must hold non-missing, ",
                "non-empty sample IDs"))
        } else {
            ## exposureData row names are checked above
            ed_samples <- as.character(rownames(ed))
            missing <- setdiff(unique(donors), ed_samples)
            if (length(missing) > 0L)
                msg <- c(msg, paste0(
                    "sample(s) in colData column '", sc,
                    "' not in exposureData: ",
                    paste(head(missing, 3), collapse = ", ")))
            extra <- setdiff(ed_samples, donors)
            if (length(extra) > 0L)
                msg <- c(msg, paste0(
                    "exposureData sample row(s) without cells: ",
                    paste(head(extra, 3), collapse = ", ")))
        }
    }

    if (length(msg) == 0L) TRUE else msg
})

# -------------------------------------------------------
# Network S4 classes for cross-omic network inference
# -------------------------------------------------------

#' CelltypeNetwork: virtual base class for cell-type networks
#'
#' @name CelltypeNetwork-class
#' @examples
#' showClass("CelltypeNetwork")
#' @exportClass CelltypeNetwork
setClass("CelltypeNetwork", contains = "VIRTUAL")

#' CelltypeNetworkResult: cell-type-resolved cross-omic network
#'
#' Stores the estimated precision matrix, adjacency matrix,
#' stability scores, and node metadata from
#' \code{\link{run_celltype_network}} and
#' \code{\link{run_exposure_network}}. Use the accessors in
#' \code{\link{network-accessors}} to extract them, and
#' \code{\link{CelltypeNetworkResult}()} to wrap a network estimated
#' elsewhere.
#'
#' @slot precision_matrix Numeric matrix. Estimated block
#'   precision (inverse covariance) matrix.
#' @slot adjacency_matrix Numeric matrix. Binary edge indicators
#'   (1 = edge present, 0 = absent).
#' @slot stability_scores Numeric matrix. Edge-level selection
#'   frequencies (0-1) across the XStARS subsamples at the selected
#'   penalties, for coglasso networks built with \code{stability = TRUE}.
#'   Empty otherwise.
#' @slot node_info \code{DataFrame}. Maps each node (row/col
#'   index) to its omic layer ("transcript" or "metabolite"),
#'   feature name, and optional module/community assignment.
#' @slot celltype Character. The cell type this network
#'   represents.
#' @slot method Character. Estimation method: \code{"coglasso"} or
#'   \code{"block_glasso"} from \code{run_celltype_network()},
#'   \code{"exposure_driven"} from \code{run_exposure_network()}, or
#'   \code{"precomputed"} for a network supplied by the user.
#' @slot metadata List. Additional parameters: n_donors,
#'   n_features, lambda values, exposure, covariates.
#'
#' @return The \code{show} method prints a summary and returns
#'   \code{NULL} invisibly.
#' @seealso \code{\link{network-accessors}},
#'   \code{\link{CelltypeNetworkResult}}
#' @name CelltypeNetworkResult-class
#' @exportClass CelltypeNetworkResult
#' @examples
#' showClass("CelltypeNetworkResult")
setClass("CelltypeNetworkResult",
    contains = "CelltypeNetwork",
    slots = c(
        precision_matrix = "matrix",
        adjacency_matrix = "matrix",
        stability_scores = "matrix",
        node_info        = "DataFrame",
        celltype         = "character",
        method           = "character",
        metadata         = "list"
    ),
    prototype = list(
        precision_matrix = matrix(nrow = 0, ncol = 0),
        adjacency_matrix = matrix(nrow = 0, ncol = 0),
        stability_scores = matrix(nrow = 0, ncol = 0),
        node_info        = S4Vectors::DataFrame(),
        celltype         = NA_character_,
        method           = "coglasso",
        metadata         = list()
    )
)

setValidity("CelltypeNetworkResult", function(object) {
    msg <- character()
    pm <- object@precision_matrix
    am <- object@adjacency_matrix

    if (nrow(pm) > 0 && ncol(pm) > 0) {
        if (nrow(pm) != ncol(pm))
            msg <- c(msg,
                "precision_matrix must be square")
        if (nrow(am) != nrow(pm) || ncol(am) != ncol(pm))
            msg <- c(msg, paste0(
                "adjacency_matrix dimensions must match ",
                "precision_matrix"))
    }

    sm <- object@stability_scores
    if (length(sm) > 0L &&
            (nrow(sm) != nrow(am) || ncol(sm) != ncol(am)))
        msg <- c(msg, paste0(
            "stability_scores must be empty or have the dimensions of ",
            "adjacency_matrix"))

    ni <- object@node_info
    if (nrow(ni) > 0 && nrow(pm) > 0) {
        if (nrow(ni) != nrow(pm))
            msg <- c(msg, paste0(
                "nrow(node_info) [", nrow(ni),
                "] must equal nrow(precision_matrix) [",
                nrow(pm), "]"))
    }
    if (nrow(ni) > 0 && !"feature" %in% colnames(ni))
        msg <- c(msg, "node_info must have a feature column")

    if (length(object@celltype) != 1L)
        msg <- c(msg, "celltype must be a single character string")
    if (length(object@method) != 1L || is.na(object@method))
        msg <- c(msg, "method must be a single non-missing string")

    if (length(msg) == 0L) TRUE else msg
})

#' NetworkComparison: differential network analysis result
#'
#' Stores results from \code{\link{run_comparative_network}},
#' including differential, shared, and cell-type-specific edges. Use the
#' accessors in \code{\link{network-accessors}} to extract them.
#'
#' @slot diff_edges \code{DataFrame}. Edges present in some but not all
#'   of the compared networks (including edges unique to one cell type),
#'   with the range of their partial correlations and, for
#'   \code{method = "fisher_z"}, p-values and adjusted p-values.
#' @slot shared_edges \code{DataFrame}. Edges present in all
#'   compared cell types.
#' @slot celltype_specific Named list of \code{DataFrame}s. Edges
#'   unique to each cell type.
#' @slot summary List. Summary statistics: Jaccard similarity,
#'   edge counts, cell types compared.
#'
#' @return The \code{show} method prints a summary and returns
#'   \code{NULL} invisibly.
#' @seealso \code{\link{network-accessors}},
#'   \code{\link{NetworkComparison}}
#' @name NetworkComparison-class
#' @exportClass NetworkComparison
#' @examples
#' showClass("NetworkComparison")
setClass("NetworkComparison",
    slots = c(
        diff_edges        = "DataFrame",
        shared_edges      = "DataFrame",
        celltype_specific = "list",
        summary           = "list"
    ),
    prototype = list(
        diff_edges        = S4Vectors::DataFrame(),
        shared_edges      = S4Vectors::DataFrame(),
        celltype_specific = list(),
        summary           = list()
    )
)

setValidity("NetworkComparison", function(object) {
    msg <- character()
    for (slot_name in c("diff_edges", "shared_edges")) {
        edges <- slot(object, slot_name)
        if (nrow(edges) > 0L &&
                !all(c("node_i", "node_j") %in% colnames(edges)))
            msg <- c(msg, paste0(
                slot_name, " must have node_i and node_j columns"))
    }
    shared <- object@shared_edges
    if ("category" %in% colnames(shared) &&
            any(as.character(shared$category) != "shared"))
        msg <- c(msg, "every edge in shared_edges must have category 'shared'")
    diff <- object@diff_edges
    if ("category" %in% colnames(diff) &&
            any(as.character(diff$category) == "shared"))
        msg <- c(msg, "diff_edges must not contain shared edges")

    specific <- object@celltype_specific
    if (length(specific) > 0L) {
        if (!all(vapply(specific, is, logical(1), "DataFrame")))
            msg <- c(msg, "celltype_specific must be a list of DataFrames")
        ct_names <- names(specific)
        if (is.null(ct_names) || anyNA(ct_names) || any(!nzchar(ct_names)) ||
                anyDuplicated(ct_names)) {
            msg <- c(msg, paste0(
                "celltype_specific must be named by cell type, ",
                "with unique names"))
        } else if (!is.null(object@summary$celltypes) &&
                !all(ct_names %in% object@summary$celltypes)) {
            msg <- c(msg, paste0(
                "names(celltype_specific) must be among ",
                "summary$celltypes"))
        }
    }
    celltypes <- object@summary$celltypes
    if (!is.null(celltypes) && !is.character(celltypes))
        msg <- c(msg, "summary$celltypes must be a character vector")

    if (length(msg) == 0L) TRUE else msg
})

# Show methods for network classes

#' @rdname CelltypeNetworkResult-class
#' @aliases show,CelltypeNetworkResult-method
#' @param object A \code{CelltypeNetworkResult} object.
#' @export
setMethod("show", "CelltypeNetworkResult", function(object) {
    p <- nrow(object@precision_matrix)
    n_edges <- sum(object@adjacency_matrix[
        upper.tri(object@adjacency_matrix)] != 0)
    ni <- object@node_info
    n_tx <- if (nrow(ni) > 0 && "omic_layer" %in% colnames(ni))
        sum(ni$omic_layer == "transcript") else 0L
    n_met <- if (nrow(ni) > 0 && "omic_layer" %in% colnames(ni))
        sum(ni$omic_layer == "metabolite") else 0L
    n_cross <- 0L
    if (n_tx > 0 && n_met > 0 && p > 0) {
        am <- object@adjacency_matrix
        layers <- ni$omic_layer
        for (i in seq_len(p - 1L)) {
            for (j in (i + 1L):p) {
                if (am[i, j] != 0 && layers[i] != layers[j])
                    n_cross <- n_cross + 1L
            }
        }
    }
    cat(sprintf(
        "CelltypeNetworkResult for '%s'\n", object@celltype))
    cat(sprintf("  Method: %s\n", object@method))
    cat(sprintf("  Nodes: %d (%d transcripts, %d metabolites)\n",
        p, n_tx, n_met))
    cat(sprintf("  Edges: %d (%d cross-omic)\n",
        n_edges, n_cross))
    if (nrow(object@stability_scores) > 0)
        cat("  Stability selection: yes\n")
    md <- object@metadata
    if (!is.null(md$n_donors))
        cat(sprintf("  Donors: %d\n", md$n_donors))
})

#' @rdname NetworkComparison-class
#' @aliases show,NetworkComparison-method
#' @param object A \code{NetworkComparison} object.
#' @export
setMethod("show", "NetworkComparison", function(object) {
    cat("NetworkComparison\n")
    cat(sprintf("  Cell types compared: %s\n",
        paste(object@summary$celltypes, collapse = ", ")))
    cat(sprintf("  Shared edges: %d\n",
        nrow(object@shared_edges)))
    cat(sprintf("  Differential edges: %d\n",
        nrow(object@diff_edges)))
    for (ct in names(object@celltype_specific)) {
        cat(sprintf("  Unique to %s: %d\n",
            ct, nrow(object@celltype_specific[[ct]])))
    }
})
