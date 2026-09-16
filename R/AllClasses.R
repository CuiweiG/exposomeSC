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
#' sample-level environmental exposure measurements and a
#' sample-to-cell mapping. This enables cell-type-specific
#' exposome-wide association studies (sc-ExWAS).
#'
#' @slot exposureData Numeric matrix. Rows = samples (donors),
#'   columns = exposure variables.
#' @slot exposureInfo \code{DataFrame}. Metadata for each
#'   exposure (family, unit, LOD).
#' @slot sampleMap \code{DataFrame}. Maps cell barcodes to
#'   sample/donor IDs.
#' @return The subsetting method returns a valid
#'   \code{SingleCellExposomeExperiment}. When cells are subset, its
#'   \code{sampleMap} is reordered to the retained cells and its
#'   \code{exposureData} is restricted to represented donors; feature
#'   subsetting leaves the donor-level slots unchanged. The \code{show} method
#'   is called for its display side effect and returns \code{NULL} invisibly.
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
        sampleMap      = "DataFrame"
    ),
    prototype = list(
        exposureData   = matrix(nrow = 0, ncol = 0),
        exposureInfo   = S4Vectors::DataFrame(),
        sampleMap      = S4Vectors::DataFrame()
    )
)

setValidity("SingleCellExposomeExperiment", function(object) {
    msg <- character()
    ed <- slot(object, "exposureData")
    ei <- slot(object, "exposureInfo")
    sm <- slot(object, "sampleMap")

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

    required_map_columns <- c("cell_id", "sample_id")
    missing_map_columns <- setdiff(required_map_columns, colnames(sm))
    if (nrow(sm) != ncol(object)) {
        msg <- c(msg, paste0(
            "nrow(sampleMap) [", nrow(sm),
            "] must equal ncol(object) [", ncol(object), "]"))
    }
    if (length(missing_map_columns)) {
        msg <- c(msg, paste0(
            "sampleMap missing required column(s): ",
            paste(missing_map_columns, collapse = ", ")))
    } else {
        map_cells <- as.character(sm$cell_id)
        map_samples <- as.character(sm$sample_id)
        if (anyNA(map_cells) || any(!nzchar(map_cells)))
            msg <- c(msg, "sampleMap$cell_id must be non-missing")
        if (anyNA(map_samples) || any(!nzchar(map_samples)))
            msg <- c(msg, "sampleMap$sample_id must be non-missing")
        if (anyDuplicated(map_cells))
            msg <- c(msg, "sampleMap$cell_id must be unique")
        if (nrow(ed) > 0 && !is.null(rownames(ed))) {
            missing <- setdiff(unique(sm$sample_id),
                               rownames(ed))
            if (length(missing) > 0)
                msg <- c(msg, paste0(
                    "sampleMap sample_id(s) not in ",
                    "exposureData: ",
                    paste(head(missing, 3), collapse = ", ")))
            extra <- setdiff(rownames(ed), unique(map_samples))
            if (length(extra) > 0)
                msg <- c(msg, paste0(
                    "exposureData sample row(s) not represented in sampleMap: ",
                    paste(head(extra, 3), collapse = ", ")))
        }
        sce_cells <- colnames(object)
        if (ncol(object) > 0 && is.null(sce_cells)) {
            msg <- c(msg, "SingleCellExposomeExperiment cells must have names")
        } else if (!is.null(sce_cells) &&
                !identical(map_cells, sce_cells)) {
            msg <- c(msg,
                "sampleMap$cell_id must match colnames(object) in order")
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
#' \code{\link{run_celltype_network}}.
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
#' @slot method Character. Estimation method used
#'   ("coglasso" or "block_glasso").
#' @slot metadata List. Additional parameters: n_donors,
#'   n_features, lambda values, exposure, covariates.
#'
#' @return A \code{CelltypeNetworkResult} object.
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
            msg <- c(msg,
                "adjacency_matrix dimensions must match ",
                "precision_matrix")
    }

    ni <- object@node_info
    if (nrow(ni) > 0 && nrow(pm) > 0) {
        if (nrow(ni) != nrow(pm))
            msg <- c(msg, paste0(
                "nrow(node_info) [", nrow(ni),
                "] must equal nrow(precision_matrix) [",
                nrow(pm), "]"))
    }

    if (length(msg) == 0L) TRUE else msg
})

#' NetworkComparison: differential network analysis result
#'
#' Stores results from \code{\link{run_comparative_network}},
#' including differential, shared, and cell-type-specific edges.
#'
#' @slot diff_edges \code{DataFrame}. Edges with significantly
#'   different partial correlations across cell types, with
#'   p-values and adjusted p-values.
#' @slot shared_edges \code{DataFrame}. Edges present in all
#'   compared cell types.
#' @slot celltype_specific List of \code{DataFrame}s. Edges
#'   unique to each cell type.
#' @slot summary List. Summary statistics: Jaccard similarity,
#'   edge counts, cell types compared.
#'
#' @return A \code{NetworkComparison} object.
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

#' TemporalNetwork: temporal network evolution result
#'
#' Stores results from \code{\link{run_temporal_network}},
#' tracking edge dynamics across time points.
#'
#' @slot networks List of \code{CelltypeNetworkResult} objects,
#'   one per time point.
#' @slot edge_dynamics \code{DataFrame}. Edge-level temporal
#'   trajectory: node_i, node_j, status at each time point.
#' @slot timepoints Character vector. Ordered time point labels.
#'
#' @return A \code{TemporalNetwork} object.
#' @name TemporalNetwork-class
#' @exportClass TemporalNetwork
#' @examples
#' showClass("TemporalNetwork")
setClass("TemporalNetwork",
    slots = c(
        networks       = "list",
        edge_dynamics  = "DataFrame",
        timepoints     = "character"
    ),
    prototype = list(
        networks       = list(),
        edge_dynamics  = S4Vectors::DataFrame(),
        timepoints     = character()
    )
)

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

#' @rdname TemporalNetwork-class
#' @aliases show,TemporalNetwork-method
#' @param object A \code{TemporalNetwork} object.
#' @export
setMethod("show", "TemporalNetwork", function(object) {
    cat("TemporalNetwork\n")
    cat(sprintf("  Time points: %s\n",
        paste(object@timepoints, collapse = " -> ")))
    cat(sprintf("  Networks: %d\n", length(object@networks)))
    if (nrow(object@edge_dynamics) > 0) {
        ed <- object@edge_dynamics
        if ("status" %in% colnames(ed)) {
            cat(sprintf("  Edge dynamics tracked: %d edges\n",
                nrow(ed)))
        }
    }
})
