# R/multi-exwas.R
# Batch ExWAS across multiple exposures

#' @include AllClasses.R
#' @include AllGenerics.R
#' @include sc-exwas.R
#' @importFrom stats p.adjust
NULL

#' Run sc-ExWAS across multiple exposures
#'
#' Convenience wrapper that runs \code{\link{run_sc_exwas}}
#' for each exposure and combines results with a final
#' global FDR correction across all exposures and cell types.
#'
#' @param x A \code{\linkS4class{SingleCellExposomeExperiment}}.
#' @param exposures Character vector. Exposure names to test.
#'   Default: all exposures in \code{exposureData}.
#' @param celltype_col Character. Column with cell type labels.
#' @param sample_col Character. Column with donor IDs.
#' @param covariates Character vector (optional).
#' @param min_cells Integer. Default 10.
#' @param min_donors Integer. Default 5.
#' @param filter_genes Logical. Default TRUE.
#' @param ... Further arguments of \code{run_sc_exwas}, such as
#'   \code{method} or \code{celltypes}. An argument that
#'   \code{run_sc_exwas} does not accept is an error.
#'
#' @return A \code{DataFrame} with all results combined and
#'   \code{padj_all} for global FDR across all tests. The \code{se} and
#'   \code{statistic} fields retain the selected backend's semantics: edgeR QL
#'   has numeric-\code{NA} coefficient standard errors and a signed square-root
#'   QL F statistic, whereas DESeq2 and voom-dream return their coefficient
#'   standard errors with Wald and moderated t statistics, respectively.
#'
#' @export
#' @examples
#' library(SingleCellExperiment); library(S4Vectors)
#' set.seed(1)
#' counts <- matrix(rpois(5000, 8), nrow = 50,
#'     dimnames = list(paste0("G", 1:50), paste0("c", 1:100)))
#' sce <- SingleCellExperiment(assays = list(counts = counts),
#'     colData = DataFrame(cell_id = paste0("c", 1:100),
#'         donor_id = rep(paste0("D", 1:5), each = 20),
#'         cell_type = rep(c("Mono", "NK"), 50)))
#' exp_mat <- matrix(rnorm(15), nrow = 5,
#'     dimnames = list(paste0("D", 1:5), c("E1", "E2", "E3")))
#' scee <- build_scee(sce, exp_mat, sample_col = "donor_id")
#' multi <- run_multi_exwas(scee,
#'     exposures = c("E1", "E2"),
#'     celltype_col = "cell_type",
#'     sample_col = "donor_id",
#'     min_cells = 5L, min_donors = 3L)
#' head(multi)
run_multi_exwas <- function(x, exposures = NULL,
                             celltype_col = "cell_type",
                             sample_col = "donor_id",
                             covariates = NULL,
                             min_cells = 10L,
                             min_donors = 5L,
                             filter_genes = TRUE, ...) {
    stopifnot(is(x, "SingleCellExposomeExperiment"))
    extra_names <- ...names()
    if (...length() && (is.null(extra_names) || any(!nzchar(extra_names)))) {
        stop("Arguments passed through ... to run_sc_exwas() must be named.")
    }
    unknown <- setdiff(extra_names, .run_sc_exwas_passthrough)
    if (length(unknown)) {
        stop("Unused argument(s) in run_multi_exwas(): ",
             paste(unknown, collapse = ", "), ".")
    }

    if (is.null(exposures)) {
        exposures <- exposureVariables(x)
    }
    if (!is.character(exposures) || !length(exposures) ||
            anyNA(exposures) || any(!nzchar(exposures)) ||
            anyDuplicated(exposures)) {
        stop(
            "exposures must be a non-empty character vector of unique, ",
            "non-missing exposure names."
        )
    }

    exp_data <- slot(x, "exposureData")
    missing <- setdiff(exposures, colnames(exp_data))
    if (length(missing) > 0) {
        stop("Exposures not found: ",
             paste(missing, collapse = ", "))
    }

    all_results <- list()
    for (exp_name in exposures) {
        message("Running sc-ExWAS for: ", exp_name)
        res <- tryCatch({
            result <- run_sc_exwas(x,
                exposure = exp_name,
                celltype_col = celltype_col,
                sample_col = sample_col,
                covariates = covariates,
                min_cells = min_cells,
                min_donors = min_donors,
                filter_genes = filter_genes, ...)
            required <- c(
                "gene", "celltype", "log2FC", "se", "statistic",
                "pvalue", "exposure"
            )
            missing_result_fields <- setdiff(required, colnames(result))
            if (length(missing_result_fields)) {
                stop(
                    "run_sc_exwas returned an incompatible result schema; ",
                    "missing: ",
                    paste(missing_result_fields, collapse = ", "),
                    "."
                )
            }
            if (nrow(result) &&
                    !all(as.character(result$exposure) == exp_name)) {
                stop("run_sc_exwas returned an incorrect exposure label.")
            }
            result
        },
            error = function(e) {
                warning("Failed for '", exp_name, "': ",
                        conditionMessage(e), call. = FALSE)
                NULL
            })
        if (!is.null(res) && nrow(res) > 0) {
            all_results <- c(all_results, list(
                as.data.frame(res)))
        }
    }

    if (length(all_results) == 0) {
        return(S4Vectors::DataFrame(
            gene = character(),
            celltype = character(),
            log2FC = numeric(),
            se = numeric(),
            statistic = numeric(),
            pvalue = numeric(),
            pvalue_underflow_clamped = logical(),
            padj = numeric(),
            padj_global = numeric(),
            exposure = character(),
            n_donors = integer(),
            n_unexposed = integer(),
            n_exposed = integer(),
            min_cells = integer(),
            median_cells = numeric(),
            method = character(),
            padj_all = numeric()))
    }

    out <- do.call(rbind, all_results)
    out$padj_all <- p.adjust(out$pvalue, method = "BH")
    S4Vectors::DataFrame(out)
}

# Arguments of the run_sc_exwas() method that run_multi_exwas() does not set
# itself. A test keeps this list equal to the method's formals.
.run_sc_exwas_passthrough <- c(
    "celltypes", "method", "min_group_donors", "filter_min_count",
    "filter_min_total_count", "filter_method", "filter_min_cpm",
    "filter_min_donors", "test_features", "robust"
)
