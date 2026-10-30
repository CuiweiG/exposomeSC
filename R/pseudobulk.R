# R/pseudobulk.R
# Donor-by-cell-type pseudobulk counts as SummarizedExperiment objects

#' @include AllClasses.R
#' @include utils.R
NULL

#' Donor pseudobulk counts for each cell type
#'
#' Sums the raw counts of each cell type within each donor and returns one
#' \code{SummarizedExperiment} per cell type, with the donors' exposure data
#' as column data. These are the pseudobulk counts that
#' \code{\link{run_sc_exwas}} and the other pseudobulk functions of the
#' package aggregate internally, so they can be passed to other
#' Bioconductor tools, inspected, or plotted.
#'
#' @param x A \code{\linkS4class{SingleCellExposomeExperiment}} with a
#'   \code{"counts"} assay. Each cell's donor is read from the
#'   \code{colData} column named by \code{\link{sampleColumn}(x)}.
#' @param celltype_col Character. Column in \code{colData} containing cell
#'   type annotations. Default \code{"cell_type"}.
#' @param celltypes Character vector (optional). Cell types to aggregate.
#'   Default: all non-missing labels, sorted.
#' @param min_cells Integer. A donor is included for a cell type only when
#'   it has at least this many cells of that type. Default 10, the default
#'   of \code{\link{run_sc_exwas}}.
#'
#' @return A \code{\link[S4Vectors]{SimpleList}} named by cell type. Each
#'   element is a \code{\link[SummarizedExperiment]{SummarizedExperiment}}
#'   with one column per included donor and
#'   \itemize{
#'     \item an assay \code{"counts"} with the summed raw counts (features
#'       by donors);
#'     \item \code{rowData} copied from \code{x};
#'     \item \code{colData} with the donor identifier (\code{donor_id}),
#'       the cell type (\code{cell_type}), the number of cells summed
#'       (\code{n_cells}) and one column per variable of
#'       \code{exposureData(x)};
#'     \item \code{metadata} recording the cell type, \code{min_cells} and
#'       the \code{colData} columns used for donors and cell types.
#'   }
#'   Cell types in which no donor reaches \code{min_cells} are left out and
#'   listed in \code{metadata(result)$skipped}.
#'
#' @details
#' The summation is done by \code{scrapper::aggregateAcrossCells()}. Donors
#' with missing exposure values are kept here; the analysis functions drop
#' them, and filter genes, when they fit a model, so a fitted model may use
#' fewer donors and genes than the object returned here. Model diagnostics
#' (\code{metadata()} of the \code{run_sc_exwas} result) list the donors and
#' gene counts actually used.
#'
#' @references
#' Lun ATL, Kancherla J (2023). Powering single-cell analyses in the browser
#' with WebAssembly. \emph{Journal of Open Source Software} 8(89):5603.
#' \doi{10.21105/joss.05603}
#'
#' @seealso \code{\link{run_sc_exwas}}
#' @export
#' @examples
#' library(SingleCellExperiment)
#' library(S4Vectors)
#' set.seed(1)
#' counts <- matrix(rpois(2000, 5), nrow = 20,
#'     dimnames = list(paste0("G", 1:20), paste0("c", 1:100)))
#' sce <- SingleCellExperiment(assays = list(counts = counts),
#'     colData = DataFrame(donor_id = rep(paste0("D", 1:5), each = 20),
#'         cell_type = rep(c("Mono", "NK"), 50)))
#' exp_mat <- matrix(rnorm(5), ncol = 1,
#'     dimnames = list(paste0("D", 1:5), "E1"))
#' scee <- build_scee(sce, exp_mat, sample_col = "donor_id")
#' pb <- donor_pseudobulk(scee, min_cells = 5L)
#' pb
#' SummarizedExperiment::colData(pb[["Mono"]])
donor_pseudobulk <- function(x, celltype_col = "cell_type",
                             celltypes = NULL, min_cells = 10L) {
    if (!methods::is(x, "SingleCellExposomeExperiment")) {
        stop("x must be a SingleCellExposomeExperiment.")
    }
    if (!is.character(celltype_col) || length(celltype_col) != 1L ||
            is.na(celltype_col)) {
        stop("celltype_col must be one column name.")
    }
    if (!is.numeric(min_cells) || length(min_cells) != 1L ||
            is.na(min_cells) || min_cells < 1) {
        stop("min_cells must be one positive integer.")
    }
    cd <- SummarizedExperiment::colData(x)
    if (!celltype_col %in% colnames(cd)) {
        stop("celltype_col '", celltype_col, "' is not a column of colData(x).")
    }
    cell_types <- as.character(cd[[celltype_col]])
    if (is.null(celltypes)) {
        celltypes <- sort(unique(cell_types))
    }
    counts <- SummarizedExperiment::assay(x, "counts")
    samples <- .cell_samples(x)
    exp_data <- exposureData(x)

    out <- stats::setNames(vector("list", length(celltypes)), celltypes)
    for (ct in celltypes) {
        aggregated <- .pseudobulk_aggregate(
            counts, samples, cell_types, ct, min_cells = min_cells
        )
        if (is.null(aggregated)) next
        donors <- aggregated$valid_donors
        donor_data <- S4Vectors::DataFrame(
            donor_id = donors,
            cell_type = rep(ct, length(donors)),
            n_cells = aggregated$n_cells,
            row.names = donors
        )
        if (ncol(exp_data)) {
            donor_data <- cbind(
                donor_data,
                S4Vectors::DataFrame(
                    exp_data[donors, , drop = FALSE],
                    check.names = FALSE
                )
            )
        }
        out[[ct]] <- SummarizedExperiment::SummarizedExperiment(
            assays = list(counts = aggregated$pb_mat),
            rowData = SummarizedExperiment::rowData(x),
            colData = donor_data,
            metadata = list(
                celltype = ct,
                min_cells = as.integer(min_cells),
                sample_column = sampleColumn(x),
                celltype_column = celltype_col
            )
        )
    }
    skipped <- celltypes[vapply(out, is.null, logical(1))]
    result <- S4Vectors::SimpleList(out[!vapply(out, is.null, logical(1))])
    S4Vectors::metadata(result)$skipped <- skipped
    result
}
