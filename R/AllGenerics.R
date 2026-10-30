# R/AllGenerics.R

#' @include AllClasses.R
#' @importFrom methods setGeneric
NULL

#' Accessors for SingleCellExposomeExperiment objects
#'
#' Get and set the donor-level exposure data of a
#' \code{\linkS4class{SingleCellExposomeExperiment}} and the column of
#' \code{colData} that links cells to donors.
#'
#' @param x A \code{\linkS4class{SingleCellExposomeExperiment}}.
#' @param value For \code{exposureData<-}, a numeric matrix with one row per
#'   sample, named by sample ID; for \code{exposureInfo<-}, a
#'   \code{DataFrame} whose \code{exposure} column equals
#'   \code{colnames(exposureData(x))}; for \code{sampleColumn<-}, the name of
#'   a \code{colData} column whose sample IDs match the rows of
#'   \code{exposureData(x)}.
#'
#' @return
#' \describe{
#'   \item{\code{exposureData()}}{Numeric matrix, samples by exposures.}
#'   \item{\code{exposureInfo()}}{\code{DataFrame} with one row per exposure.}
#'   \item{\code{exposureVariables()}}{Character vector of exposure names.}
#'   \item{\code{sampleColumn()}}{Character scalar: the \code{colData}
#'     column holding each cell's sample (donor) ID. The analysis functions
#'     read each cell's donor from this column.}
#'   \item{\code{cellSampleMap()}}{\code{DataFrame} with columns
#'     \code{cell_id} (the cell names) and \code{sample_id}, derived from
#'     that column.}
#'   \item{Replacement functions}{The modified object, after a validity
#'     check.}
#' }
#'
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
#' sampleColumn(scee)
#' head(cellSampleMap(scee))
#'
#' ## Replace the exposure matrix, here rescaling PM2.5
#' exp_new <- exposureData(scee)
#' exp_new[, "PM2.5"] <- exp_new[, "PM2.5"] / 10
#' exposureData(scee) <- exp_new
#' exposureData(scee)[, "PM2.5"]
setGeneric("exposureData", function(x)
    standardGeneric("exposureData"))

#' @export
#' @rdname SCEE-accessors
setGeneric("exposureInfo", function(x)
    standardGeneric("exposureInfo"))

#' @export
#' @rdname SCEE-accessors
setGeneric("cellSampleMap", function(x)
    standardGeneric("cellSampleMap"))

#' @export
#' @rdname SCEE-accessors
setGeneric("sampleColumn", function(x)
    standardGeneric("sampleColumn"))

#' @export
#' @rdname SCEE-accessors
setGeneric("sampleColumn<-", function(x, value)
    standardGeneric("sampleColumn<-"))

#' @export
#' @rdname SCEE-accessors
setGeneric("exposureVariables", function(x)
    standardGeneric("exposureVariables"))

#' @export
#' @rdname SCEE-accessors
setGeneric("exposureData<-", function(x, value)
    standardGeneric("exposureData<-"))

#' @export
#' @rdname SCEE-accessors
setGeneric("exposureInfo<-", function(x, value)
    standardGeneric("exposureInfo<-"))

#' Accessors for network results
#'
#' Extract the parts of a \code{\linkS4class{CelltypeNetworkResult}}
#' (returned by \code{\link{run_celltype_network}} and
#' \code{\link{run_exposure_network}}) and of a
#' \code{\linkS4class{NetworkComparison}} (returned by
#' \code{\link{run_comparative_network}}).
#'
#' @param x A \code{CelltypeNetworkResult} or, for \code{diffEdges()},
#'   \code{sharedEdges()}, \code{celltypeSpecificEdges()} and
#'   \code{comparisonSummary()}, a \code{NetworkComparison}.
#'
#' @return
#' \describe{
#'   \item{\code{networkPrecision()}}{Numeric feature-by-feature precision
#'     matrix.}
#'   \item{\code{networkAdjacency()}}{Feature-by-feature matrix of edge
#'     indicators (1 = edge, 0 = no edge).}
#'   \item{\code{networkStability()}}{Feature-by-feature matrix of edge
#'     selection frequencies, or a 0 by 0 matrix when none were stored.}
#'   \item{\code{networkNodes()}}{\code{DataFrame} with one row per feature,
#'     in the order of the matrix rows.}
#'   \item{\code{networkCelltype()}}{Character scalar: the cell type.}
#'   \item{\code{networkMethod()}}{Character scalar: the estimation method.}
#'   \item{\code{metadata()}}{List of parameters and summaries, for example
#'     \code{n_donors}, \code{transform} and \code{edge_counts}.}
#'   \item{\code{diffEdges()}}{\code{DataFrame} of edges that are unique to
#'     one cell type or present in some but not all networks.}
#'   \item{\code{sharedEdges()}}{\code{DataFrame} of edges present in every
#'     network.}
#'   \item{\code{celltypeSpecificEdges()}}{Named list of \code{DataFrame}s,
#'     one per cell type with edges found only in that cell type.}
#'   \item{\code{comparisonSummary()}}{List with the compared cell types,
#'     edge counts, pairwise Jaccard overlap and the comparison method.}
#' }
#'
#' @name network-accessors
#' @aliases networkPrecision
#' @export
#' @examples
#' features <- c("G1", "G2", "M1")
#' make_network <- function(celltype, cross_weight) {
#'     precision <- diag(3)
#'     precision[1, 3] <- precision[3, 1] <- cross_weight
#'     adjacency <- (precision != 0) * 1L
#'     diag(adjacency) <- 0L
#'     dimnames(precision) <- dimnames(adjacency) <- list(features, features)
#'     CelltypeNetworkResult(
#'         precision = precision,
#'         adjacency = adjacency,
#'         node_info = S4Vectors::DataFrame(
#'             feature = features,
#'             omic_layer = c("transcript", "transcript", "metabolite")),
#'         celltype = celltype,
#'         metadata = list(n_donors = 20L))
#' }
#' mono <- make_network("Monocyte", -0.45)
#' networkCelltype(mono)
#' networkAdjacency(mono)
#' networkNodes(mono)
#' metadata(mono)$n_donors
#'
#' comparison <- run_comparative_network(mono, make_network("B cell", 0))
#' diffEdges(comparison)
#' celltypeSpecificEdges(comparison)
#' comparisonSummary(comparison)$jaccard
setGeneric("networkPrecision", function(x)
    standardGeneric("networkPrecision"))

#' @rdname network-accessors
#' @export
setGeneric("networkAdjacency", function(x)
    standardGeneric("networkAdjacency"))

#' @rdname network-accessors
#' @export
setGeneric("networkStability", function(x)
    standardGeneric("networkStability"))

#' @rdname network-accessors
#' @export
setGeneric("networkNodes", function(x)
    standardGeneric("networkNodes"))

#' @rdname network-accessors
#' @export
setGeneric("networkCelltype", function(x)
    standardGeneric("networkCelltype"))

#' @rdname network-accessors
#' @export
setGeneric("networkMethod", function(x)
    standardGeneric("networkMethod"))

#' @rdname network-accessors
#' @export
setGeneric("diffEdges", function(x)
    standardGeneric("diffEdges"))

#' @rdname network-accessors
#' @export
setGeneric("sharedEdges", function(x)
    standardGeneric("sharedEdges"))

#' @rdname network-accessors
#' @export
setGeneric("celltypeSpecificEdges", function(x)
    standardGeneric("celltypeSpecificEdges"))

#' @rdname network-accessors
#' @export
setGeneric("comparisonSummary", function(x)
    standardGeneric("comparisonSummary"))

#' Run cell-type-specific ExWAS
#' @param x A \code{\linkS4class{SingleCellExposomeExperiment}}.
#' @param ... Additional arguments.
#' @export
#' @rdname run_sc_exwas
setGeneric("run_sc_exwas", function(x, ...)
    standardGeneric("run_sc_exwas"))
