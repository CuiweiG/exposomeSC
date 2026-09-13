# R/build-scee.R
# Construct SingleCellExposomeExperiment from components

#' @include AllClasses.R
#' @importFrom S4Vectors DataFrame
#' @importFrom SummarizedExperiment colData
#' @importFrom SingleCellExperiment SingleCellExperiment
NULL

#' Build a SingleCellExposomeExperiment
#'
#' Combines a \code{SingleCellExperiment} with sample-level
#' exposure measurements into a unified container for sc-ExWAS.
#'
#' @param sce A \code{\link[SingleCellExperiment:SingleCellExperiment-class]{SingleCellExperiment}}.
#' @param exposure_matrix Numeric matrix. Rows = samples,
#'   columns = exposure variables, with unique row and column names.
#'   The row names must be exactly the sample IDs in
#'   \code{colData(sce)[[sample_col]]}: every sample of the SCE must
#'   have a row, and rows for samples absent from the SCE are an
#'   error. \code{NA} is allowed; \code{NaN} and infinite values are
#'   not.
#' @param sample_col Character. Column in \code{colData(sce)}
#'   containing the sample/donor ID.
#' @param exposure_info \code{DataFrame} (optional). Metadata
#'   for exposures (for example family, unit and LOD). It must contain
#'   an \code{exposure} column equal to
#'   \code{colnames(exposure_matrix)}, in the same order.
#'
#' @return A \code{\linkS4class{SingleCellExposomeExperiment}}.
#'
#' @export
#' @examples
#' library(SingleCellExperiment)
#' library(S4Vectors)
#' set.seed(1)
#' counts <- matrix(rpois(2000, 5), nrow = 50,
#'     dimnames = list(paste0("G", 1:50), paste0("c", 1:40)))
#' sce <- SingleCellExperiment(
#'     assays = list(counts = counts),
#'     colData = DataFrame(
#'         cell_id = paste0("c", seq_len(40)),
#'         donor_id = rep(paste0("D", 1:4), each = 10),
#'         cell_type = rep(c("Mono", "NK"), 20)))
#' exp_mat <- matrix(rnorm(12), nrow = 4,
#'     dimnames = list(paste0("D", 1:4),
#'                     c("PM2.5", "Pb", "BPA")))
#' scee <- build_scee(sce, exp_mat, sample_col = "donor_id")
#' scee
build_scee <- function(sce, exposure_matrix, sample_col,
                        exposure_info = NULL) {
    if (!is(sce, "SingleCellExperiment")) {
        stop("sce must inherit from SingleCellExperiment.")
    }
    if (!is.matrix(exposure_matrix) || !is.numeric(exposure_matrix)) {
        stop("exposure_matrix must be a numeric matrix.")
    }
    if (!nrow(exposure_matrix) || !ncol(exposure_matrix)) {
        stop("exposure_matrix must have at least one row and one column.")
    }
    if (any(is.infinite(exposure_matrix)) || any(is.nan(exposure_matrix))) {
        stop("exposure_matrix may contain NA, but not NaN or infinite values.")
    }
    cd <- colData(sce)
    if (!sample_col %in% colnames(cd)) {
        stop("'", sample_col, "' not found in colData. ",
             "Available: ",
             paste(colnames(cd), collapse = ", "))
    }

    cell_samples <- as.character(cd[[sample_col]])
    if (anyNA(cell_samples) || any(!nzchar(cell_samples))) {
        stop("Sample IDs in colData(sce) must be non-missing and non-empty.")
    }
    exp_samples <- rownames(exposure_matrix)

    if (is.null(exp_samples) || anyNA(exp_samples) ||
            any(!nzchar(exp_samples)) || anyDuplicated(exp_samples)) {
        stop(
            "exposure_matrix must have unique, non-missing sample row names."
        )
    }
    exposure_names <- colnames(exposure_matrix)
    if (is.null(exposure_names) || anyNA(exposure_names) ||
            any(!nzchar(exposure_names)) || anyDuplicated(exposure_names)) {
        stop(
            "exposure_matrix must have unique, non-missing exposure ",
            "column names."
        )
    }

    missing <- setdiff(unique(cell_samples), exp_samples)
    if (length(missing) > 0) {
        stop("Samples in SCE not found in exposure_matrix: ",
             paste(missing, collapse = ", "))
    }
    extra <- setdiff(exp_samples, unique(cell_samples))
    if (length(extra) > 0) {
        stop(
            "Samples in exposure_matrix not represented in SCE: ",
            paste(extra, collapse = ", ")
        )
    }
    if (is.null(colnames(sce)) || anyNA(colnames(sce)) ||
            any(!nzchar(colnames(sce))) || anyDuplicated(colnames(sce))) {
        stop("sce must have unique, non-missing cell column names.")
    }

    sample_map <- S4Vectors::DataFrame(
        cell_id = colnames(sce),
        sample_id = cell_samples
    )

    if (is.null(exposure_info)) {
        exposure_info <- S4Vectors::DataFrame(
            exposure = colnames(exposure_matrix),
            family = rep(NA_character_,
                         ncol(exposure_matrix)),
            unit = rep(NA_character_,
                       ncol(exposure_matrix))
        )
    } else {
        if (!is(exposure_info, "DataFrame")) {
            stop("exposure_info must be an S4Vectors::DataFrame.")
        }
        if (!"exposure" %in% colnames(exposure_info) ||
                !identical(
                    as.character(exposure_info$exposure),
                    exposure_names
                )) {
            stop(
                "exposure_info$exposure must match exposure_matrix ",
                "column names in order."
            )
        }
    }

    new("SingleCellExposomeExperiment",
        sce,
        exposureData = exposure_matrix,
        exposureInfo = exposure_info,
        sampleMap = sample_map)
}
