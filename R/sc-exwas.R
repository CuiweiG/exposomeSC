# R/sc-exwas.R
# Core: cell-type-specific exposome-wide association

#' @include AllClasses.R
#' @include AllGenerics.R
#' @include utils.R
#' @importFrom S4Vectors DataFrame
#' @importFrom SummarizedExperiment assay colData
#' @importFrom stats p.adjust as.formula
NULL

#' Cell-type-specific exposome-wide association study
#'
#' Performs pseudobulk differential expression analysis with
#' a continuous exposure variable as covariate. Aggregates
#' counts within each donor-celltype stratum to avoid
#' pseudoreplication, then uses an empirical-Bayes pseudobulk model for
#' inference.
#'
#' @param x A \code{\linkS4class{SingleCellExposomeExperiment}}.
#' @param exposure Character. Name of the exposure variable
#'   (column in \code{exposureData}).
#' @param celltype_col Character. Column in \code{colData}
#'   containing cell type annotations.
#' @param celltypes Character vector (optional). Which cell
#'   types to test. Default: all.
#' @param sample_col Character. Column containing donor IDs.
#' @param covariates Character vector (optional). Additional
#'   covariates from \code{exposureData}.
#' @param min_cells Integer. Minimum cells per donor-celltype
#'   stratum. Default 10.
#' @param min_donors Integer. Minimum donors per cell type.
#'   Default 5. Cell types with fewer valid donors are skipped
#'   with a warning.
#' @param filter_genes Logical. If \code{TRUE} (default), filter
#'   low-expression genes before model fitting. The edgeR backend defaults to
#'   exposure-label-independent donor-support filtering; the DESeq2 and
#'   dreamlet backends retain genes detected in at least 20\% of donors.
#' @param method Character. Backend for differential expression.
#'   \code{"edgeR"} (default) uses robust quasi-likelihood inference with
#'   empirical-Bayes dispersion shrinkage. \code{"DESeq2"} (and its retained
#'   alias \code{"pseudobulk"}) uses DESeq2 on aggregated counts.
#'   \code{"dreamlet"} uses variancePartition's
#'   \code{voomWithDreamWeights()} and \code{dream()} with precision
#'   weights from the mean-variance trend. Pseudobulk data have one
#'   observation per donor and the model has no random effects, so
#'   \code{dream()} fits it with limma and the moderated t statistic uses
#'   the residual degrees of freedom; it does not call the
#'   \pkg{dreamlet} package.
#' @param min_group_donors Integer. For a binary exposure, minimum complete
#'   donors required in each group for the edgeR backend. Default 5.
#' @param filter_min_count Numeric. Minimum count passed to
#'   \code{edgeR::filterByExpr} when \code{filter_method = "filterByExpr"}.
#' @param filter_min_total_count Numeric. Minimum total pseudobulk count for
#'   either edgeR filtering method.
#' @param filter_method Character. edgeR expression-filtering rule.
#'   \code{"fixed_support"} (default) requires a fixed CPM threshold in a
#'   fixed number of independent donors and does not use exposure labels.
#'   \code{"filterByExpr"} retains the edgeR compatibility option.
#' @param filter_min_cpm Numeric. CPM threshold for
#'   \code{filter_method = "fixed_support"}. Default 1.
#' @param filter_min_donors Integer. Minimum number of independent donors
#'   meeting \code{filter_min_cpm}. Default 5.
#' @param test_features Optional character vector of feature identifiers that
#'   defines the reported hypothesis family. Library-size calculation,
#'   expression filtering, normalisation and dispersion estimation still use
#'   all available features; subsetting occurs only after model fitting. Every
#'   requested feature must pass the pre-specified expression filter in every
#'   tested cell type, otherwise the call fails rather than silently changing
#'   the requested hypothesis family. Currently supported only by the edgeR
#'   backend; other backends fail explicitly rather than ignoring it.
#' @param robust Logical. Use robust empirical-Bayes dispersion and
#'   quasi-likelihood fitting in edgeR. Default \code{TRUE}.
#' @param ... Must be empty. Any further argument is an error, so that a
#'   misspelled argument name cannot silently change the analysis.
#'
#' @return A \code{DataFrame} with common columns including gene, celltype,
#'   log2FC, se, statistic, pvalue, padj (within-cell-type BH),
#'   padj_global (cross-cell-type BH), exposure, donor-support diagnostics and
#'   method. Backend-specific abundance summaries are reported as logCPM for
#'   edgeR/voom-dream and baseMean for DESeq2. Uncertainty fields are
#'   backend-specific: edgeR QL reports \code{se} as numeric \code{NA} and
#'   \code{statistic} as the directional signed square root of the QL F
#'   statistic; DESeq2 reports its coefficient standard error and Wald
#'   statistic; voom-dream reports its coefficient standard error and
#'   moderated t statistic.
#'
#' @details
#' The pseudobulk approach aggregates single-cell counts
#' within each (donor x cell type) stratum, creating one
#' observation per donor per cell type. This correctly
#' treats cells from the same donor as non-independent,
#' avoiding the inflated Type I error rates that plague
#' cell-level analyses (Squair et al. 2021 *Nat Commun*).
#'
#' The exposure variable and specified covariates enter one donor-level design
#' matrix. The default edgeR backend estimates robust negative-binomial
#' dispersions and uses quasi-likelihood F-tests. Complete-transcriptome library
#' sizes are retained when expression filtering or \code{test_features}
#' restricts the reported hypothesis family.
#'
#' \code{edgeR::glmQLFTest} does not return a coefficient standard error.
#' Consequently, the edgeR backend never reconstructs \code{se} from the ratio
#' of the log-fold change to the QL statistic. Its \code{statistic} field is a
#' directional ranking representation of the QL F statistic, not a Wald z or
#' t statistic, and must not be used to derive a coefficient confidence
#' interval. The edgeR result metadata records
#' \code{se_method = "not_available_edgeR_QL"} and
#' \code{statistic_type = "signed_sqrt_qlf"}.
#'
#' @references
#' Squair JW et al. (2021). Confronting false discoveries in
#' single-cell differential expression. \emph{Nat Commun}
#' 12:5692. \doi{10.1038/s41467-021-25960-2}
#'
#' Chen Y, Lun ATL, Smyth GK (2016). From reads to genes to pathways:
#' differential expression analysis of RNA-Seq experiments using Rsubread and
#' the edgeR quasi-likelihood pipeline. \emph{F1000Research} 5:1438.
#' \doi{10.12688/f1000research.8987.2}
#'
#' @export
#' @rdname run_sc_exwas
#' @examples
#' library(SingleCellExperiment)
#' library(S4Vectors)
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
#' \donttest{
#' if (requireNamespace("edgeR", quietly = TRUE)) {
#'     result <- run_sc_exwas(scee, exposure = "E1",
#'         celltype_col = "cell_type", sample_col = "donor_id",
#'         min_cells = 5L)
#'     head(result)
#' }
#' }
setMethod("run_sc_exwas",
    "SingleCellExposomeExperiment",
    function(x, exposure, celltype_col, celltypes = NULL,
             sample_col = "donor_id",
             covariates = NULL,
             min_cells = 10L,
             min_donors = 5L,
             filter_genes = TRUE,
             method = c("edgeR", "DESeq2", "pseudobulk", "dreamlet"),
             min_group_donors = 5L,
             filter_min_count = 10L,
             filter_min_total_count = 15L,
             filter_method = c("fixed_support", "filterByExpr"),
             filter_min_cpm = 1,
             filter_min_donors = 5L,
             test_features = NULL,
             robust = TRUE,
             ...) {

    method <- match.arg(method)
    .stop_on_unused_dots("run_sc_exwas", ...)

    if (method != "edgeR" && !is.null(test_features)) {
        stop(
            "test_features is currently supported only by method='edgeR'; ",
            "refusing to ignore the requested hypothesis family."
        )
    }

    if (method == "edgeR") {
        filter_method <- match.arg(filter_method)
        return(.run_sc_exwas_edger(
            x = x,
            exposure = exposure,
            celltype_col = celltype_col,
            celltypes = celltypes,
            sample_col = sample_col,
            covariates = covariates,
            min_cells = min_cells,
            min_donors = min_donors,
            min_group_donors = min_group_donors,
            filter_genes = filter_genes,
            filter_min_count = filter_min_count,
            filter_min_total_count = filter_min_total_count,
            filter_method = filter_method,
            filter_min_cpm = filter_min_cpm,
            filter_min_donors = filter_min_donors,
            test_features = test_features,
            robust = robust
        ))
    }

    if (method == "dreamlet") {
        return(.run_sc_exwas_dreamlet(
            x, exposure, celltype_col, celltypes,
            sample_col, covariates, min_cells, min_donors,
            filter_genes))
    }

    if (!requireNamespace("DESeq2", quietly = TRUE)) {
        stop("Package 'DESeq2' is required. ",
             "BiocManager::install('DESeq2')")
    }

    exp_data <- exposureData(x)
    required_fields <- c(exposure, covariates)
    missing_fields <- setdiff(required_fields, colnames(exp_data))
    if (length(missing_fields)) {
        stop(
            "Exposure/covariate column(s) not found in exposureData: ",
            paste(missing_fields, collapse = ", ")
        )
    }

    cd <- SummarizedExperiment::colData(x)
    if (!all(c(sample_col, celltype_col) %in% colnames(cd))) {
        stop("sample_col and celltype_col must be present in colData(x).")
    }
    counts_mat <- SummarizedExperiment::assay(x, "counts")
    samples <- as.character(cd[[sample_col]])
    cell_types <- as.character(cd[[celltype_col]])

    if (is.null(celltypes)) {
        celltypes <- sort(unique(cell_types))
    }

    all_results <- list()

    for (ct in celltypes) {
        ## Pseudobulk aggregation via shared utility
        pb_result <- .pseudobulk_aggregate(
            counts_mat, samples, cell_types, ct,
            min_cells = min_cells)

        if (is.null(pb_result)) {
            warning("Skipping ", ct, ": no donors with >= ",
                    min_cells, " cells", call. = FALSE)
            next
        }

        valid_donors <- pb_result$valid_donors
        pb_mat <- pb_result$pb_mat
        donor_cells <- pb_result$n_cells

        if (length(valid_donors) < min_donors) {
            warning("Skipping ", ct, ": only ",
                    length(valid_donors), " donors with >= ",
                    min_cells, " cells (need >= ", min_donors,
                    ")", call. = FALSE)
            next
        }
        if (length(valid_donors) < 8L) {
            warning("Cell type '", ct, "': only ",
                    length(valid_donors),
                    " donors -- DESeq2 dispersion estimates ",
                    "may be unreliable. Consider ",
                    "interpreting results cautiously.",
                    call. = FALSE)
        }

        ## Build design data
        design_df <- data.frame(
            exposure = as.numeric(exp_data[valid_donors, exposure]),
            row.names = valid_donors,
            check.names = FALSE
        )
        for (cov in covariates) {
            design_df[[cov]] <- as.numeric(exp_data[valid_donors, cov])
        }
        complete <- stats::complete.cases(design_df)
        valid_donors <- valid_donors[complete]
        design_df <- design_df[complete, , drop = FALSE]
        pb_mat <- pb_mat[, complete, drop = FALSE]
        donor_cells <- donor_cells[complete]

        if (length(valid_donors) < min_donors) {
            warning(
                "Skipping ", ct, ": only ", length(valid_donors),
                " complete donors with >= ", min_cells, " cells (need >= ",
                min_donors, ").",
                call. = FALSE
            )
            next
        }
        exposure_levels <- sort(unique(design_df$exposure))
        group_counts <- if (length(exposure_levels) == 2L) {
            table(factor(design_df$exposure, levels = exposure_levels))
        } else {
            integer()
        }

        ## Gene filtering
        if (filter_genes) {
            ## Keep genes detected (count > 0) in >= 20% of
            ## donors -- reduces multiple testing burden and
            ## improves dispersion estimation
            min_samples <- max(2L,
                ceiling(0.2 * ncol(pb_mat)))
            keep <- rowSums(pb_mat > 0) >= min_samples
        } else {
            keep <- rowSums(pb_mat) > 0
        }
        pb_mat <- pb_mat[keep, , drop = FALSE]

        if (nrow(pb_mat) == 0) {
            warning("Skipping ", ct, ": no gene passed expression filtering.",
                    call. = FALSE)
            next
        }

        ## DESeq2
        design_formula <- stats::reformulate(c("exposure", covariates))

        tryCatch({
            dds <- DESeq2::DESeqDataSetFromMatrix(
                countData = pb_mat,
                colData = design_df,
                design = design_formula)
            dds <- DESeq2::DESeq(dds, quiet = TRUE)
            res <- DESeq2::results(dds, name = "exposure")
            res_df <- .sc_exwas_result_frame(
                gene = rownames(res),
                celltype = ct,
                log2FC = res$log2FoldChange,
                se = res$lfcSE,
                statistic = res$stat,
                pvalue = res$pvalue,
                exposure = exposure,
                n_donors = length(valid_donors),
                n_unexposed = if (length(group_counts)) {
                    unname(group_counts[[1L]])
                } else NA_integer_,
                n_exposed = if (length(group_counts)) {
                    unname(group_counts[[2L]])
                } else NA_integer_,
                min_cells = min(donor_cells),
                median_cells = stats::median(donor_cells),
                method = "DESeq2_Wald"
            )
            res_df$baseMean <- as.numeric(res$baseMean)
            all_results <- c(all_results, list(res_df))
        }, error = function(e) {
            warning("DESeq2 failed for ", ct, ": ",
                    conditionMessage(e), call. = FALSE)
        })
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
            baseMean = numeric(),
            exposure = character(),
            n_donors = integer(),
            n_unexposed = integer(),
            n_exposed = integer(),
            min_cells = integer(),
            median_cells = numeric(),
            method = character()))
    }

    out <- do.call(rbind, all_results)
    out$padj_global <- p.adjust(out$pvalue, method = "BH")
    S4Vectors::DataFrame(out)
})


# Robust edgeR quasi-likelihood backend (internal)
#
# Library sizes are computed from the complete pseudobulk transcriptome before
# expression filtering. Filtering therefore changes only the tested hypothesis
# universe, never the offset for a retained gene.
.run_sc_exwas_edger <- function(x, exposure, celltype_col, celltypes,
                                sample_col, covariates, min_cells,
                                min_donors, min_group_donors,
                                filter_genes, filter_min_count,
                                filter_min_total_count, filter_method,
                                filter_min_cpm, filter_min_donors,
                                test_features, robust) {
    if (!requireNamespace("edgeR", quietly = TRUE)) {
        stop("Package 'edgeR' is required for method='edgeR'.")
    }
    if (!is.numeric(min_cells) || length(min_cells) != 1L || min_cells < 1) {
        stop("min_cells must be one positive integer.")
    }
    if (!is.numeric(min_donors) || length(min_donors) != 1L ||
            min_donors < 3) {
        stop("min_donors must be one integer >= 3.")
    }
    if (!is.numeric(min_group_donors) || length(min_group_donors) != 1L ||
            min_group_donors < 2) {
        stop("min_group_donors must be one integer >= 2.")
    }
    if (!is.logical(filter_genes) || length(filter_genes) != 1L ||
            is.na(filter_genes)) {
        stop("filter_genes must be TRUE or FALSE.")
    }
    filter_method <- match.arg(
        filter_method,
        choices = c("fixed_support", "filterByExpr")
    )
    if (!is.numeric(filter_min_cpm) || length(filter_min_cpm) != 1L ||
            !is.finite(filter_min_cpm) || filter_min_cpm <= 0) {
        stop("filter_min_cpm must be one finite number > 0.")
    }
    if (!is.numeric(filter_min_donors) ||
            length(filter_min_donors) != 1L ||
            !is.finite(filter_min_donors) || filter_min_donors < 1 ||
            filter_min_donors != as.integer(filter_min_donors)) {
        stop("filter_min_donors must be one positive integer.")
    }
    if (!is.numeric(filter_min_count) || length(filter_min_count) != 1L ||
            !is.finite(filter_min_count) || filter_min_count < 0) {
        stop("filter_min_count must be one finite non-negative number.")
    }
    if (!is.numeric(filter_min_total_count) ||
            length(filter_min_total_count) != 1L ||
            !is.finite(filter_min_total_count) ||
            filter_min_total_count < 0) {
        stop(
            "filter_min_total_count must be one finite non-negative number."
        )
    }
    if (!is.logical(robust) || length(robust) != 1L || is.na(robust)) {
        stop("robust must be TRUE or FALSE.")
    }

    exp_data <- exposureData(x)
    required <- c(exposure, covariates)
    missing_fields <- setdiff(required, colnames(exp_data))
    if (length(missing_fields)) {
        stop(
            "Exposure/covariate column(s) not found in exposureData: ",
            paste(missing_fields, collapse = ", ")
        )
    }

    cd <- SummarizedExperiment::colData(x)
    if (!all(c(sample_col, celltype_col) %in% colnames(cd))) {
        stop("sample_col and celltype_col must be present in colData(x).")
    }
    counts_mat <- SummarizedExperiment::assay(x, "counts")
    feature_ids <- rownames(counts_mat)
    if (!is.null(test_features)) {
        if (!is.character(test_features) || !length(test_features) ||
                anyNA(test_features) || any(!nzchar(test_features)) ||
                anyDuplicated(test_features)) {
            stop(
                "test_features must be NULL or a non-empty character vector ",
                "of unique, non-missing feature identifiers."
            )
        }
        if (is.null(feature_ids) || anyDuplicated(feature_ids)) {
            stop(
                "Unique feature row names are required when test_features ",
                "is supplied."
            )
        }
        missing_features <- setdiff(test_features, feature_ids)
        if (length(missing_features)) {
            stop(
                "test_features not found in rownames(x): ",
                paste(missing_features, collapse = ", ")
            )
        }
    }
    samples <- as.character(cd[[sample_col]])
    cell_types <- as.character(cd[[celltype_col]])
    if (is.null(celltypes)) {
        celltypes <- sort(unique(cell_types))
    }

    formula_terms <- c("exposure", covariates)
    design_formula <- stats::reformulate(formula_terms)
    all_results <- list()
    diagnostics <- list()

    for (celltype in celltypes) {
        aggregated <- .pseudobulk_aggregate(
            counts_mat,
            samples,
            cell_types,
            celltype,
            min_cells = as.integer(min_cells)
        )
        if (is.null(aggregated)) {
            warning(
                "Skipping ", celltype, ": no donors with >= ",
                min_cells, " cells.",
                call. = FALSE
            )
            next
        }

        donors <- aggregated$valid_donors
        design_data <- data.frame(
            exposure = as.numeric(exp_data[donors, exposure]),
            row.names = donors,
            check.names = FALSE
        )
        for (covariate in covariates) {
            design_data[[covariate]] <- as.numeric(
                exp_data[donors, covariate]
            )
        }
        complete <- stats::complete.cases(design_data)
        donors <- donors[complete]
        design_data <- design_data[complete, , drop = FALSE]
        pseudobulk <- aggregated$pb_mat[, complete, drop = FALSE]
        donor_cells <- aggregated$n_cells[complete]

        if (length(donors) < min_donors) {
            warning(
                "Skipping ", celltype, ": only ", length(donors),
                " complete donors with >= ", min_cells, " cells.",
                call. = FALSE
            )
            next
        }
        exposure_levels <- sort(unique(design_data$exposure))
        group_counts <- if (length(exposure_levels) == 2L) {
            table(factor(design_data$exposure, levels = exposure_levels))
        } else {
            integer()
        }
        if (length(group_counts) && any(group_counts < min_group_donors)) {
            warning(
                "Skipping ", celltype, ": binary exposure group has fewer than ",
                min_group_donors, " complete donors.",
                call. = FALSE
            )
            next
        }

        design_full <- stats::model.matrix(design_formula, data = design_data)
        protected <- match(c("(Intercept)", "exposure"), colnames(design_full))
        if (anyNA(protected) || qr(design_full[, protected, drop = FALSE])$rank < 2L) {
            warning(
                "Skipping ", celltype,
                ": the exposure is not identifiable.",
                call. = FALSE
            )
            next
        }
        keep_design <- protected
        nuisance <- setdiff(seq_len(ncol(design_full)), protected)
        for (column in nuisance) {
            candidate <- cbind(
                design_full[, keep_design, drop = FALSE],
                design_full[, column, drop = FALSE]
            )
            if (qr(candidate)$rank > length(keep_design)) {
                keep_design <- c(keep_design, column)
            }
        }
        design <- design_full[, keep_design, drop = FALSE]
        dropped_design <- setdiff(colnames(design_full), colnames(design))
        if (nrow(design) <= ncol(design) + 2L) {
            warning(
                "Skipping ", celltype,
                ": insufficient residual degrees of freedom.",
                call. = FALSE
            )
            next
        }

        full_library_size <- colSums(pseudobulk)
        dge <- edgeR::DGEList(
            counts = pseudobulk,
            lib.size = pmax(full_library_size, 1)
        )
        keep <- if (isTRUE(filter_genes)) {
            if (identical(filter_method, "fixed_support")) {
                donor_cpm <- edgeR::cpm(
                    dge,
                    normalized.lib.sizes = FALSE,
                    log = FALSE
                )
                rowSums(donor_cpm >= filter_min_cpm) >=
                    as.integer(filter_min_donors) &
                    rowSums(pseudobulk) >= filter_min_total_count
            } else {
                edgeR::filterByExpr(
                    dge,
                    design = design,
                    min.count = filter_min_count,
                    min.total.count = filter_min_total_count
                )
            }
        } else {
            rowSums(pseudobulk) > 0
        }
        if (!any(keep)) {
            warning(
                "Skipping ", celltype, ": no gene passed expression filtering.",
                call. = FALSE
            )
            next
        }
        dge <- dge[keep, , keep.lib.sizes = TRUE]
        dge <- edgeR::normLibSizes(dge, method = "TMM")

        fit_result <- tryCatch({
            dge <- edgeR::estimateDisp(
                dge,
                design,
                robust = robust
            )
            fit <- edgeR::glmQLFit(
                dge,
                design,
                robust = robust
            )
            test <- edgeR::glmQLFTest(fit, coef = "exposure")
            list(
                fit = fit,
                table = edgeR::topTags(
                    test,
                    n = Inf,
                    sort.by = "none"
                )$table
            )
        }, error = function(e) {
            warning(
                "edgeR failed for ", celltype, ": ",
                conditionMessage(e),
                call. = FALSE
            )
            NULL
        })
        if (is.null(fit_result)) {
            next
        }

        fitted_result <- fit_result$table
        if (!is.null(test_features)) {
            filtered_requested <- setdiff(
                test_features,
                rownames(fitted_result)
            )
            if (length(filtered_requested)) {
                stop(
                    "Requested test_features did not pass expression ",
                    "filtering in cell type '", celltype, "': ",
                    paste(filtered_requested, collapse = ", "),
                    ". Adjust the pre-specified support rule explicitly ",
                    "or report the target as not estimable."
                )
            }
        }
        tested <- if (is.null(test_features)) {
            rep(TRUE, nrow(fitted_result))
        } else {
            rownames(fitted_result) %in% test_features
        }
        result <- fitted_result[tested, , drop = FALSE]
        pvalue_underflow <- is.finite(result$PValue) & result$PValue == 0
        reported_pvalue <- result$PValue
        reported_pvalue[pvalue_underflow] <- .Machine$double.xmin
        signed_statistic <- sign(result$logFC) * sqrt(pmax(result$F, 0))
        ## glmQLFTest reports a quasi-likelihood F statistic, not a Wald
        ## statistic, and does not provide a coefficient standard error.
        ## Retain a signed square root only as a directional ranking statistic;
        ## it must never be inverted to manufacture an SE or confidence interval.
        standard_error <- rep(NA_real_, nrow(result))
        result_data <- data.frame(
            gene = rownames(result),
            celltype = celltype,
            log2FC = result$logFC,
            se = standard_error,
            statistic = signed_statistic,
            pvalue = reported_pvalue,
            pvalue_underflow_clamped = pvalue_underflow,
            padj = stats::p.adjust(reported_pvalue, method = "BH"),
            logCPM = result$logCPM,
            exposure = exposure,
            n_donors = length(donors),
            n_unexposed = if (length(group_counts)) unname(group_counts[[1]]) else NA_integer_,
            n_exposed = if (length(group_counts)) unname(group_counts[[2]]) else NA_integer_,
            min_cells = min(donor_cells),
            median_cells = stats::median(donor_cells),
            method = if (isTRUE(robust)) {
                "edgeR_robust_QL"
            } else {
                "edgeR_QL"
            },
            stringsAsFactors = FALSE
        )
        all_results[[celltype]] <- result_data
        diagnostics[[celltype]] <- list(
            n_donors = length(donors),
            donor_ids = donors,
            n_genes_input = nrow(pseudobulk),
            n_genes_fit = nrow(fitted_result),
            n_genes_tested = nrow(result_data),
            test_feature_ids = rownames(result),
            n_test_features_requested = if (is.null(test_features)) {
                nrow(fitted_result)
            } else {
                length(test_features)
            },
            filter_method = if (isTRUE(filter_genes)) {
                filter_method
            } else {
                "nonzero"
            },
            filter_min_cpm = filter_min_cpm,
            filter_min_donors = as.integer(filter_min_donors),
            filter_min_count = filter_min_count,
            filter_min_total_count = filter_min_total_count,
            full_library_size = full_library_size,
            normalisation_factors = dge$samples$norm.factors,
            design_columns = colnames(design),
            dropped_design_columns = dropped_design,
            design_rank = qr(design)$rank,
            residual_df = min(
                if (!is.null(fit_result$fit$df.residual.zeros)) {
                    fit_result$fit$df.residual.zeros
                } else {
                    fit_result$fit$df.residual
                },
                na.rm = TRUE
            ),
            min_cells = min(donor_cells),
            median_cells = stats::median(donor_cells)
        )
    }

    result_parameters <- list(
        exposure = exposure,
        covariates = covariates,
        min_cells = as.integer(min_cells),
        min_donors = as.integer(min_donors),
        min_group_donors = as.integer(min_group_donors),
        filter_genes = filter_genes,
        filter_method = filter_method,
        filter_min_cpm = filter_min_cpm,
        filter_min_donors = as.integer(filter_min_donors),
        filter_min_count = filter_min_count,
        filter_min_total_count = filter_min_total_count,
        test_features = test_features,
        n_test_features_requested = if (is.null(test_features)) {
            NULL
        } else {
            length(test_features)
        },
        robust = robust,
        schema_version = "exposomeSC_sc_exwas_edger_v2",
        se_method = "not_available_edgeR_QL",
        statistic_type = "signed_sqrt_qlf",
        multiplicity = "BH across all tested gene-by-cell-type hypotheses"
    )

    if (!length(all_results)) {
        result <- S4Vectors::DataFrame(
            gene = character(),
            celltype = character(),
            log2FC = numeric(),
            se = numeric(),
            statistic = numeric(),
            pvalue = numeric(),
            pvalue_underflow_clamped = logical(),
            padj = numeric(),
            padj_global = numeric(),
            logCPM = numeric(),
            exposure = character(),
            n_donors = integer(),
            n_unexposed = integer(),
            n_exposed = integer(),
            min_cells = integer(),
            median_cells = numeric(),
            method = character()
        )
        S4Vectors::metadata(result)$diagnostics <- diagnostics
        S4Vectors::metadata(result)$parameters <- result_parameters
        return(result)
    }

    output <- do.call(rbind, all_results)
    rownames(output) <- NULL
    output$padj_global <- stats::p.adjust(output$pvalue, method = "BH")
    output <- output[order(output$padj_global, output$pvalue), , drop = FALSE]
    result <- S4Vectors::DataFrame(output)
    S4Vectors::metadata(result)$diagnostics <- diagnostics
    S4Vectors::metadata(result)$parameters <- result_parameters
    result
}


# -------------------------------------------------------
# dreamlet backend (internal)
# Uses the variancePartition voom-dream functions:
#   1. Pseudobulk aggregation per donor x celltype
#   2. voomWithDreamWeights: precision weights from the mean-variance trend
#   3. dream(): ~ exposure + covariates. With one observation per donor there
#      are no random effects, so dream() delegates to limma and the moderated
#      t statistic uses the residual degrees of freedom.
# -------------------------------------------------------
.run_sc_exwas_dreamlet <- function(x, exposure, celltype_col,
    celltypes, sample_col, covariates, min_cells, min_donors,
    filter_genes) {

    for (pkg in c("variancePartition", "limma", "edgeR")) {
        if (!requireNamespace(pkg, quietly = TRUE))
            stop("Package '", pkg, "' required for method='dreamlet'. ",
                 "BiocManager::install('", pkg, "')")
    }

    exp_data <- exposureData(x)
    required_fields <- c(exposure, covariates)
    missing_fields <- setdiff(required_fields, colnames(exp_data))
    if (length(missing_fields)) {
        stop(
            "Exposure/covariate column(s) not found in exposureData: ",
            paste(missing_fields, collapse = ", ")
        )
    }

    cd <- SummarizedExperiment::colData(x)
    if (!all(c(sample_col, celltype_col) %in% colnames(cd))) {
        stop("sample_col and celltype_col must be present in colData(x).")
    }
    counts_mat <- SummarizedExperiment::assay(x, "counts")
    samples <- as.character(cd[[sample_col]])
    cell_types <- as.character(cd[[celltype_col]])

    if (is.null(celltypes))
        celltypes <- sort(unique(cell_types))

    all_results <- list()

    for (ct in celltypes) {
        ## Pseudobulk aggregation
        pb_result <- .pseudobulk_aggregate(
            counts_mat, samples, cell_types, ct,
            min_cells = min_cells)

        if (is.null(pb_result)) {
            warning("Skipping ", ct, ": no donors with >= ",
                    min_cells, " cells", call. = FALSE)
            next
        }

        valid_donors <- pb_result$valid_donors
        pb_mat <- pb_result$pb_mat
        n_cells_per_donor <- pb_result$n_cells

        if (length(valid_donors) < min_donors) {
            warning("Skipping ", ct, ": only ",
                    length(valid_donors), " donors",
                    call. = FALSE)
            next
        }

        ## Build metadata and retain only donors with a complete design.
        meta_df <- data.frame(
            exposure = as.numeric(exp_data[valid_donors, exposure]),
            n_cells = n_cells_per_donor,
            row.names = valid_donors,
            check.names = FALSE
        )
        for (cov in covariates) {
            meta_df[[cov]] <- as.numeric(exp_data[valid_donors, cov])
        }
        complete <- stats::complete.cases(
            meta_df[, c("exposure", covariates), drop = FALSE]
        )
        valid_donors <- valid_donors[complete]
        meta_df <- meta_df[complete, , drop = FALSE]
        pb_mat <- pb_mat[, complete, drop = FALSE]
        n_cells_per_donor <- n_cells_per_donor[complete]

        if (length(valid_donors) < min_donors) {
            warning(
                "Skipping ", ct, ": only ", length(valid_donors),
                " complete donors with >= ", min_cells, " cells (need >= ",
                min_donors, ").",
                call. = FALSE
            )
            next
        }
        exposure_levels <- sort(unique(meta_df$exposure))
        group_counts <- if (length(exposure_levels) == 2L) {
            table(factor(meta_df$exposure, levels = exposure_levels))
        } else {
            integer()
        }

        ## Gene filtering
        if (filter_genes) {
            min_samples <- max(2L, ceiling(0.2 * ncol(pb_mat)))
            keep <- rowSums(pb_mat > 0) >= min_samples
            pb_mat <- pb_mat[keep, , drop = FALSE]
        }
        if (nrow(pb_mat) == 0) {
            warning("Skipping ", ct, ": no gene passed expression filtering.",
                    call. = FALSE)
            next
        }

        ## Design formula
        form <- stats::reformulate(c("exposure", covariates))

        tryCatch({
            ## Create DGEList for voom
            dge <- edgeR::DGEList(counts = pb_mat)
            dge <- edgeR::normLibSizes(dge)

            ## voomWithDreamWeights: precision weights from the
            ## mean-variance trend
            vobj <- variancePartition::voomWithDreamWeights(
                dge, form, meta_df, plot = FALSE)

            ## dream(): no random effects in a one-row-per-donor design,
            ## so this is a weighted limma fit
            fit <- variancePartition::dream(
                vobj, form, meta_df)
            fit <- limma::eBayes(fit)

            ## Extract results for exposure coefficient
            tt <- limma::topTable(
                fit, coef = "exposure",
                number = Inf, sort.by = "none")

            standard_error <- ifelse(tt$t != 0,
                abs(tt$logFC / tt$t), NA_real_)

            res_df <- .sc_exwas_result_frame(
                gene = rownames(tt),
                celltype = ct,
                log2FC = tt$logFC,
                se = standard_error,
                statistic = tt$t,
                pvalue = tt$P.Value,
                exposure = exposure,
                n_donors = length(valid_donors),
                n_unexposed = if (length(group_counts)) {
                    unname(group_counts[[1L]])
                } else NA_integer_,
                n_exposed = if (length(group_counts)) {
                    unname(group_counts[[2L]])
                } else NA_integer_,
                min_cells = min(n_cells_per_donor),
                median_cells = stats::median(n_cells_per_donor),
                method = "voom_dream"
            )
            res_df$logCPM <- rowMeans(
                vobj$E[rownames(tt), , drop = FALSE]
            )
            all_results <- c(all_results, list(res_df))

        }, error = function(e) {
            warning("voom-dream failed for ", ct, ": ",
                    conditionMessage(e), call. = FALSE)
        })
    }

    if (length(all_results) == 0) {
        return(S4Vectors::DataFrame(
            gene = character(), celltype = character(),
            log2FC = numeric(), se = numeric(),
            statistic = numeric(), pvalue = numeric(),
            pvalue_underflow_clamped = logical(),
            padj = numeric(), padj_global = numeric(),
            logCPM = numeric(),
            exposure = character(), n_donors = integer(),
            n_unexposed = integer(), n_exposed = integer(),
            min_cells = integer(), median_cells = numeric(),
            method = character()))
    }

    out <- do.call(rbind, all_results)
    out$padj_global <- p.adjust(out$pvalue, method = "BH")
    S4Vectors::DataFrame(out)
}
