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
#' @param x A \code{\linkS4class{SingleCellExposomeExperiment}}. Each
#'   cell's donor is read from the \code{colData} column named by
#'   \code{\link{sampleColumn}(x)}.
#' @param exposure Character. Name of the exposure variable
#'   (column in \code{exposureData}).
#' @param celltype_col Character. Column in \code{colData}
#'   containing cell type annotations. Default \code{"cell_type"}.
#' @param celltypes Character vector (optional). Which cell
#'   types to test. Default: all.
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
#'   voom-dream backends retain genes detected in at least 20\% of donors.
#' @param method Character. Backend for differential expression.
#'   \code{"edgeR"} (default) uses the edgeR v4 quasi-likelihood pipeline
#'   (see Details). \code{"DESeq2"} uses DESeq2 on aggregated counts.
#'   \code{"voom_dream"} uses variancePartition's
#'   \code{voomWithDreamWeights()} and \code{dream()} with precision
#'   weights from the mean-variance trend (Hoffman and Roussos 2021).
#'   Pseudobulk data have one observation per donor and the model has no
#'   random effects, so \code{dream()} fits it with limma and the
#'   moderated t statistic uses the residual degrees of freedom. This
#'   backend does not use the \pkg{dreamlet} package.
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
#' @param robust Logical. For the edgeR backend, moderate the
#'   quasi-likelihood dispersions with robust empirical Bayes
#'   (\code{glmQLFit(robust = TRUE)}). Default \code{TRUE}.
#' @param ... Must be empty. Any further argument is an error, so that a
#'   misspelled argument name cannot silently change the analysis.
#'
#' @return A \code{DataFrame} with one row per tested gene and cell type and
#'   the same columns, in the same order, for every backend: \code{gene},
#'   \code{celltype}, \code{log2FC}, \code{se}, \code{statistic},
#'   \code{pvalue}, \code{pvalue_underflow_clamped} (\code{TRUE} where a
#'   p-value of exactly zero was replaced by the smallest positive double),
#'   \code{padj} (BH within the cell type), \code{padj_global} (BH across all
#'   tested gene-by-cell-type hypotheses), \code{logCPM}, \code{exposure},
#'   \code{n_donors}, \code{n_unexposed} and \code{n_exposed} (group sizes
#'   when the exposure takes two values, otherwise \code{NA}),
#'   \code{min_cells} and \code{median_cells} (cells per donor in the
#'   cell type), and \code{method}. Rows are ordered by \code{padj_global},
#'   then \code{pvalue}.
#'
#'   \code{logCPM} is the average log2 count per million reported by the
#'   backend: edgeR's \code{logCPM} for the edgeR backend, the mean of the
#'   voom log-CPM values for voom-dream, and \code{edgeR::aveLogCPM()} of the
#'   filtered pseudobulk counts for DESeq2.
#'
#'   Uncertainty fields are backend-specific: edgeR QL reports \code{se} as
#'   numeric \code{NA} and \code{statistic} as the directional signed square
#'   root of the QL F statistic; DESeq2 reports its coefficient standard error
#'   and Wald statistic; voom-dream reports the moderated standard error
#'   (\code{|log2FC / t|}) and the moderated t statistic.
#'   \code{metadata(result)$parameters} records the settings, the backend and
#'   its package version, \code{se_method} and \code{statistic_type};
#'   \code{metadata(result)$diagnostics} holds, for each fitted cell type,
#'   the donors used, the numbers of genes before and after filtering and the
#'   design columns.
#'
#' @details
#' The pseudobulk approach aggregates single-cell counts
#' within each (donor x cell type) stratum, creating one
#' observation per donor per cell type. This correctly
#' treats cells from the same donor as non-independent,
#' avoiding the inflated Type I error rates that plague
#' cell-level analyses (Squair et al. 2021, \emph{Nat Commun}).
#' The pseudobulk counts that enter the models are available as
#' \code{SummarizedExperiment} objects from \code{\link{donor_pseudobulk}}.
#'
#' The exposure variable and specified covariates enter one donor-level design
#' matrix. The default edgeR backend fits the edgeR v4 quasi-likelihood
#' model, \code{glmQLFit(legacy = FALSE)}, which sets the negative-binomial
#' dispersion internally from the most highly expressed genes and moderates
#' the quasi-likelihood dispersions by empirical Bayes; the exposure
#' coefficient is tested with a quasi-likelihood F-test (Chen et al. 2025).
#' No separate \code{estimateDisp()} step is run: from edgeR 4.10.0 the
#' method gives the same result whether or not the data hold such estimates,
#' whereas earlier edgeR 4.x releases did not. The edgeR version and
#' \code{legacy = FALSE} are recorded in \code{metadata(result)$parameters}.
#' Complete-transcriptome library sizes are retained when expression
#' filtering or \code{test_features} restricts the reported hypothesis
#' family.
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
#' \strong{What is delegated and what is kept.} The raw counts of each
#' donor and cell type are summed with
#' \code{scrapper::aggregateAcrossCells()}, and the models are fitted by
#' \pkg{edgeR}, \pkg{DESeq2} or \pkg{variancePartition} and \pkg{limma}.
#' Four choices are made by \pkg{exposomeSC} itself, on purpose, because
#' they define the hypotheses of an exposome-wide study:
#' \enumerate{
#'   \item The exposure is a donor-level variable from
#'     \code{exposureData(x)}, continuous or binary, and enters the design
#'     matrix together with the covariates; cells are never grouped by
#'     exposure.
#'   \item Gene filtering does not use the exposure. The default edgeR rule
#'     (\code{filter_method = "fixed_support"}) keeps genes with at least
#'     \code{filter_min_cpm} CPM in \code{filter_min_donors} donors, and the
#'     DESeq2 and voom-dream backends keep genes detected in at least 20\%
#'     of donors (and in at least two); which genes are tested therefore
#'     does not depend on the exposure values. By contrast,
#'     \code{scran::pseudoBulkDGE()}, for example, filters each label with
#'     \code{edgeR::filterByExpr()} using the design matrix.
#'   \item In the edgeR backend the library sizes are the column sums of the
#'     complete pseudobulk transcriptome, computed before filtering, and are
#'     kept when filtering or \code{test_features} reduces the tested genes.
#'     The voom-dream backend computes library sizes from the retained genes
#'     and DESeq2 uses its own size factors.
#'   \item \code{padj_global} applies the Benjamini-Hochberg procedure once
#'     across all tested gene-by-cell-type hypotheses, next to the
#'     within-cell-type \code{padj}.
#' }
#'
#' @references
#' Squair JW et al. (2021). Confronting false discoveries in
#' single-cell differential expression. \emph{Nat Commun}
#' 12:5692. \doi{10.1038/s41467-021-25960-2}
#'
#' Chen Y, Chen L, Lun ATL, Baldoni PL, Smyth GK (2025). edgeR v4: powerful
#' differential analysis of sequencing data with expanded functionality and
#' improved support for small counts and larger datasets. \emph{Nucleic Acids
#' Research} 53(2):gkaf018. \doi{10.1093/nar/gkaf018}
#'
#' Hoffman GE, Roussos P (2021). Dream: powerful differential expression
#' analysis for repeated measures designs. \emph{Bioinformatics}
#' 37(2):192--201. \doi{10.1093/bioinformatics/btaa687}
#'
#' Hoffman GE, Schadt EE (2016). variancePartition: interpreting drivers of
#' variation in complex gene expression studies. \emph{BMC Bioinformatics}
#' 17:483. \doi{10.1186/s12859-016-1323-z}
#'
#' @seealso \code{\link{donor_pseudobulk}} for the pseudobulk counts,
#'   \code{\link{run_multi_exwas}} for several exposures.
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
#' result <- run_sc_exwas(scee, exposure = "E1",
#'     celltype_col = "cell_type", min_cells = 5L)
#' head(result)
#' metadata(result)$parameters$statistic_type
setMethod("run_sc_exwas",
    "SingleCellExposomeExperiment",
    function(x, exposure, celltype_col = "cell_type", celltypes = NULL,
             covariates = NULL,
             min_cells = 10L,
             min_donors = 5L,
             filter_genes = TRUE,
             method = c("edgeR", "DESeq2", "voom_dream"),
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
    }
    switch(method,
        edgeR = .run_sc_exwas_edger(
            x = x,
            exposure = exposure,
            celltype_col = celltype_col,
            celltypes = celltypes,
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
        ),
        DESeq2 = .run_sc_exwas_deseq2(
            x, exposure, celltype_col, celltypes, covariates,
            min_cells, min_donors, filter_genes
        ),
        voom_dream = .run_sc_exwas_voom_dream(
            x, exposure, celltype_col, celltypes, covariates,
            min_cells, min_donors, filter_genes
        )
    )
})


# Settings recorded in metadata(result)$parameters by every backend.
.sc_exwas_common_parameters <- function(method, exposure, covariates,
                                         min_cells, min_donors,
                                         filter_genes, se_method,
                                         statistic_type, backend_version) {
    list(
        method = method,
        exposure = exposure,
        covariates = covariates,
        min_cells = as.integer(min_cells),
        min_donors = as.integer(min_donors),
        filter_genes = filter_genes,
        backend_version = backend_version,
        aggregation = "sum of raw counts per donor and cell type",
        schema_version = "exposomeSC_sc_exwas_v3",
        se_method = se_method,
        statistic_type = statistic_type,
        multiplicity = "BH across all tested gene-by-cell-type hypotheses"
    )
}

# Expression filter shared by the DESeq2 and voom-dream backends: genes
# detected (count > 0) in at least 20% of donors and in at least two.
.detected_in_donors <- function(pb_mat, filter_genes) {
    if (filter_genes) {
        min_samples <- max(2L, ceiling(0.2 * ncol(pb_mat)))
        rowSums(pb_mat > 0) >= min_samples
    } else {
        rowSums(pb_mat) > 0
    }
}

# Aggregate one cell type and keep donors with a complete design. Returns
# NULL, with the warnings the backends have always given, when the cell type
# has too few donors.
.sc_exwas_celltype_data <- function(inputs, celltype, exposure, covariates,
                                    min_cells, min_donors) {
    aggregated <- .pseudobulk_aggregate(
        inputs$counts, inputs$samples, inputs$cell_types, celltype,
        min_cells = min_cells
    )
    if (is.null(aggregated)) {
        warning("Skipping ", celltype, ": no donors with >= ",
                min_cells, " cells", call. = FALSE)
        return(NULL)
    }
    if (length(aggregated$valid_donors) < min_donors) {
        warning("Skipping ", celltype, ": only ",
                length(aggregated$valid_donors), " donors with >= ",
                min_cells, " cells (need >= ", min_donors, ")",
                call. = FALSE)
        return(NULL)
    }
    design <- .sc_exwas_design_data(
        inputs$exp_data, aggregated$valid_donors, exposure, covariates
    )
    donors <- aggregated$valid_donors[design$complete]
    if (length(donors) < min_donors) {
        warning(
            "Skipping ", celltype, ": only ", length(donors),
            " complete donors with >= ", min_cells, " cells (need >= ",
            min_donors, ").",
            call. = FALSE
        )
        return(NULL)
    }
    list(
        donors = donors,
        pb_mat = aggregated$pb_mat[, design$complete, drop = FALSE],
        n_cells = aggregated$n_cells[design$complete],
        design = design
    )
}

# Per-cell-type diagnostics of the DESeq2 and voom-dream backends.
.sc_exwas_basic_diagnostics <- function(stratum, n_genes_tested,
                                        filter_genes, design_columns) {
    list(
        n_donors = length(stratum$donors),
        donor_ids = stratum$donors,
        n_genes_input = nrow(stratum$pb_mat),
        n_genes_tested = n_genes_tested,
        filter_method = if (isTRUE(filter_genes)) {
            "detected_in_20pct_of_donors"
        } else {
            "nonzero"
        },
        design_columns = design_columns,
        min_cells = min(stratum$n_cells),
        median_cells = stats::median(stratum$n_cells)
    )
}


# DESeq2 backend (internal)
.run_sc_exwas_deseq2 <- function(x, exposure, celltype_col, celltypes,
                                 covariates, min_cells, min_donors,
                                 filter_genes) {
    if (!requireNamespace("DESeq2", quietly = TRUE)) {
        stop("Package 'DESeq2' is required for method = 'DESeq2'. ",
             "Install it with BiocManager::install('DESeq2').")
    }
    inputs <- .sc_exwas_inputs(x, exposure, covariates, celltype_col,
                               celltypes)
    celltypes <- inputs$celltypes
    design_formula <- stats::reformulate(c("exposure", covariates))
    tables <- stats::setNames(vector("list", length(celltypes)), celltypes)
    diagnostics <- list()

    for (ct in celltypes) {
        stratum <- .sc_exwas_celltype_data(
            inputs, ct, exposure, covariates, min_cells, min_donors
        )
        if (is.null(stratum)) next
        if (length(stratum$donors) < 8L) {
            warning("Cell type '", ct, "': only ",
                    length(stratum$donors),
                    " donors -- DESeq2 dispersion estimates ",
                    "may be unreliable. Consider ",
                    "interpreting results cautiously.",
                    call. = FALSE)
        }
        pb_mat <- stratum$pb_mat[
            .detected_in_donors(stratum$pb_mat, filter_genes), ,
            drop = FALSE
        ]
        if (nrow(pb_mat) == 0) {
            warning("Skipping ", ct, ": no gene passed expression filtering.",
                    call. = FALSE)
            next
        }

        table_ct <- tryCatch({
            dds <- DESeq2::DESeqDataSetFromMatrix(
                countData = pb_mat,
                colData = stratum$design$data,
                design = design_formula)
            dds <- DESeq2::DESeq(dds, quiet = TRUE)
            res <- DESeq2::results(dds, name = "exposure")
            .sc_exwas_result_frame(
                gene = rownames(res),
                celltype = ct,
                log2FC = res$log2FoldChange,
                se = res$lfcSE,
                statistic = res$stat,
                pvalue = res$pvalue,
                logCPM = edgeR::aveLogCPM(pb_mat[rownames(res), ,
                                                 drop = FALSE]),
                exposure = exposure,
                n_donors = length(stratum$donors),
                n_unexposed = stratum$design$n_unexposed,
                n_exposed = stratum$design$n_exposed,
                min_cells = min(stratum$n_cells),
                median_cells = stats::median(stratum$n_cells),
                method = "DESeq2_Wald"
            )
        }, error = function(e) {
            warning("DESeq2 failed for ", ct, ": ",
                    conditionMessage(e), call. = FALSE)
            NULL
        })
        if (is.null(table_ct)) next
        tables[[ct]] <- table_ct
        diagnostics[[ct]] <- .sc_exwas_basic_diagnostics(
            stratum, nrow(table_ct), filter_genes,
            c("(Intercept)", "exposure", covariates)
        )
    }

    .sc_exwas_finalise(
        tables, diagnostics,
        .sc_exwas_common_parameters(
            method = "DESeq2", exposure = exposure, covariates = covariates,
            min_cells = min_cells, min_donors = min_donors,
            filter_genes = filter_genes,
            se_method = "DESeq2_lfcSE",
            statistic_type = "wald_z",
            backend_version = c(
                DESeq2 = as.character(utils::packageVersion("DESeq2"))
            )
        )
    )
}


# Robust edgeR quasi-likelihood backend (internal)
#
# Library sizes are computed from the complete pseudobulk transcriptome before
# expression filtering. Filtering therefore changes only the tested hypothesis
# universe, never the offset for a retained gene.
.run_sc_exwas_edger <- function(x, exposure, celltype_col, celltypes,
                                covariates, min_cells,
                                min_donors, min_group_donors,
                                filter_genes, filter_min_count,
                                filter_min_total_count, filter_method,
                                filter_min_cpm, filter_min_donors,
                                test_features, robust) {
    .check_edger_arguments(
        min_cells = min_cells, min_donors = min_donors,
        min_group_donors = min_group_donors, filter_genes = filter_genes,
        filter_min_cpm = filter_min_cpm,
        filter_min_donors = filter_min_donors,
        filter_min_count = filter_min_count,
        filter_min_total_count = filter_min_total_count, robust = robust
    )
    filter_method <- match.arg(
        filter_method,
        choices = c("fixed_support", "filterByExpr")
    )
    inputs <- .sc_exwas_inputs(x, exposure, covariates, celltype_col,
                               celltypes)
    .check_test_features(test_features, rownames(inputs$counts))
    celltypes <- inputs$celltypes

    design_formula <- stats::reformulate(c("exposure", covariates))
    tables <- stats::setNames(vector("list", length(celltypes)), celltypes)
    diagnostics <- list()

    for (celltype in celltypes) {
        aggregated <- .pseudobulk_aggregate(
            inputs$counts,
            inputs$samples,
            inputs$cell_types,
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

        design_info <- .sc_exwas_design_data(
            inputs$exp_data, aggregated$valid_donors, exposure, covariates
        )
        donors <- aggregated$valid_donors[design_info$complete]
        design_data <- design_info$data
        pseudobulk <- aggregated$pb_mat[, design_info$complete, drop = FALSE]
        donor_cells <- aggregated$n_cells[design_info$complete]

        if (length(donors) < min_donors) {
            warning(
                "Skipping ", celltype, ": only ", length(donors),
                " complete donors with >= ", min_cells, " cells.",
                call. = FALSE
            )
            next
        }
        group_counts <- design_info$group_counts
        if (length(group_counts) && any(group_counts < min_group_donors)) {
            warning(
                "Skipping ", celltype,
                ": binary exposure group has fewer than ",
                min_group_donors, " complete donors.",
                call. = FALSE
            )
            next
        }

        design_full <- stats::model.matrix(design_formula, data = design_data)
        design <- .identifiable_design(design_full)
        if (is.null(design)) {
            warning(
                "Skipping ", celltype,
                ": the exposure is not identifiable.",
                call. = FALSE
            )
            next
        }
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
        keep <- .edger_keep(
            dge, pseudobulk, design, filter_genes, filter_method,
            filter_min_cpm, filter_min_donors, filter_min_count,
            filter_min_total_count
        )
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
            ## edgeR v4 quasi-likelihood: the NB dispersion is set inside
            ## glmQLFit(), so no estimateDisp() step is needed.
            fit <- edgeR::glmQLFit(
                dge,
                design,
                robust = robust,
                legacy = FALSE
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
        ## glmQLFTest reports a quasi-likelihood F statistic, not a Wald
        ## statistic, and does not provide a coefficient standard error.
        ## Retain a signed square root only as a directional ranking statistic;
        ## it must never be inverted to manufacture an SE or confidence
        ## interval.
        tables[[celltype]] <- .sc_exwas_result_frame(
            gene = rownames(result),
            celltype = celltype,
            log2FC = result$logFC,
            se = rep(NA_real_, nrow(result)),
            statistic = sign(result$logFC) * sqrt(pmax(result$F, 0)),
            pvalue = result$PValue,
            logCPM = result$logCPM,
            exposure = exposure,
            n_donors = length(donors),
            n_unexposed = design_info$n_unexposed,
            n_exposed = design_info$n_exposed,
            min_cells = min(donor_cells),
            median_cells = stats::median(donor_cells),
            method = if (isTRUE(robust)) "edgeR_robust_QL" else "edgeR_QL"
        )
        diagnostics[[celltype]] <- list(
            n_donors = length(donors),
            donor_ids = donors,
            n_genes_input = nrow(pseudobulk),
            n_genes_fit = nrow(fitted_result),
            n_genes_tested = nrow(result),
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

    parameters <- c(
        .sc_exwas_common_parameters(
            method = "edgeR", exposure = exposure, covariates = covariates,
            min_cells = min_cells, min_donors = min_donors,
            filter_genes = filter_genes,
            se_method = "not_available_edgeR_QL",
            statistic_type = "signed_sqrt_qlf",
            backend_version = c(
                edgeR = as.character(utils::packageVersion("edgeR"))
            )
        ),
        list(
            min_group_donors = as.integer(min_group_donors),
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
            edger_version = as.character(utils::packageVersion("edgeR")),
            ql_legacy = FALSE
        )
    )
    .sc_exwas_finalise(tables, diagnostics, parameters)
}

# Argument checks of the edgeR backend.
.check_edger_arguments <- function(min_cells, min_donors, min_group_donors,
                                   filter_genes, filter_min_cpm,
                                   filter_min_donors, filter_min_count,
                                   filter_min_total_count, robust) {
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
    invisible(NULL)
}

# test_features must name unique rows of x.
.check_test_features <- function(test_features, feature_ids) {
    if (is.null(test_features)) return(invisible(NULL))
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
    invisible(NULL)
}

# Keep the intercept and exposure columns and add each nuisance column only
# when it raises the rank. Returns NULL when the exposure is not
# identifiable.
.identifiable_design <- function(design_full) {
    protected <- match(c("(Intercept)", "exposure"), colnames(design_full))
    if (anyNA(protected) ||
            qr(design_full[, protected, drop = FALSE])$rank < 2L) {
        return(NULL)
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
    design_full[, keep_design, drop = FALSE]
}

# Expression filter of the edgeR backend.
.edger_keep <- function(dge, pseudobulk, design, filter_genes, filter_method,
                        filter_min_cpm, filter_min_donors, filter_min_count,
                        filter_min_total_count) {
    if (!isTRUE(filter_genes)) {
        return(rowSums(pseudobulk) > 0)
    }
    if (identical(filter_method, "fixed_support")) {
        donor_cpm <- edgeR::cpm(
            dge,
            normalized.lib.sizes = FALSE,
            log = FALSE
        )
        return(rowSums(donor_cpm >= filter_min_cpm) >=
            as.integer(filter_min_donors) &
            rowSums(pseudobulk) >= filter_min_total_count)
    }
    edgeR::filterByExpr(
        dge,
        design = design,
        min.count = filter_min_count,
        min.total.count = filter_min_total_count
    )
}


# -------------------------------------------------------
# voom-dream backend (internal)
# Uses the variancePartition voom-dream functions:
#   1. Pseudobulk aggregation per donor x celltype
#   2. voomWithDreamWeights: precision weights from the mean-variance trend
#   3. dream(): ~ exposure + covariates. With one observation per donor there
#      are no random effects, so dream() delegates to limma and the moderated
#      t statistic uses the residual degrees of freedom.
# -------------------------------------------------------
.run_sc_exwas_voom_dream <- function(x, exposure, celltype_col,
    celltypes, covariates, min_cells, min_donors, filter_genes) {

    if (!requireNamespace("variancePartition", quietly = TRUE)) {
        stop("Package 'variancePartition' is required for ",
             "method = 'voom_dream'. Install it with ",
             "BiocManager::install('variancePartition').")
    }
    inputs <- .sc_exwas_inputs(x, exposure, covariates, celltype_col,
                               celltypes)
    celltypes <- inputs$celltypes
    form <- stats::reformulate(c("exposure", covariates))
    tables <- stats::setNames(vector("list", length(celltypes)), celltypes)
    diagnostics <- list()

    for (ct in celltypes) {
        stratum <- .sc_exwas_celltype_data(
            inputs, ct, exposure, covariates, min_cells, min_donors
        )
        if (is.null(stratum)) next
        meta_df <- stratum$design$data
        pb_mat <- stratum$pb_mat
        if (filter_genes) {
            pb_mat <- pb_mat[.detected_in_donors(pb_mat, TRUE), ,
                             drop = FALSE]
        }
        if (nrow(pb_mat) == 0) {
            warning("Skipping ", ct, ": no gene passed expression filtering.",
                    call. = FALSE)
            next
        }

        table_ct <- tryCatch({
            dge <- edgeR::normLibSizes(edgeR::DGEList(counts = pb_mat))
            ## voomWithDreamWeights: precision weights from the
            ## mean-variance trend
            vobj <- variancePartition::voomWithDreamWeights(
                dge, form, meta_df, plot = FALSE)
            ## dream(): no random effects in a one-row-per-donor design,
            ## so this is a weighted limma fit
            fit <- limma::eBayes(variancePartition::dream(vobj, form, meta_df))
            tt <- limma::topTable(
                fit, coef = "exposure",
                number = Inf, sort.by = "none")
            .sc_exwas_result_frame(
                gene = rownames(tt),
                celltype = ct,
                log2FC = tt$logFC,
                se = ifelse(tt$t != 0, abs(tt$logFC / tt$t), NA_real_),
                statistic = tt$t,
                pvalue = tt$P.Value,
                logCPM = rowMeans(vobj$E[rownames(tt), , drop = FALSE]),
                exposure = exposure,
                n_donors = length(stratum$donors),
                n_unexposed = stratum$design$n_unexposed,
                n_exposed = stratum$design$n_exposed,
                min_cells = min(stratum$n_cells),
                median_cells = stats::median(stratum$n_cells),
                method = "voom_dream"
            )
        }, error = function(e) {
            warning("voom-dream failed for ", ct, ": ",
                    conditionMessage(e), call. = FALSE)
            NULL
        })
        if (is.null(table_ct)) next
        tables[[ct]] <- table_ct
        diagnostics[[ct]] <- .sc_exwas_basic_diagnostics(
            stratum, nrow(table_ct), filter_genes,
            c("(Intercept)", "exposure", covariates)
        )
    }

    .sc_exwas_finalise(
        tables, diagnostics,
        .sc_exwas_common_parameters(
            method = "voom_dream", exposure = exposure,
            covariates = covariates, min_cells = min_cells,
            min_donors = min_donors, filter_genes = filter_genes,
            se_method = "limma_moderated_from_t",
            statistic_type = "moderated_t",
            backend_version = c(
                variancePartition = as.character(
                    utils::packageVersion("variancePartition")
                ),
                limma = as.character(utils::packageVersion("limma"))
            )
        )
    )
}
