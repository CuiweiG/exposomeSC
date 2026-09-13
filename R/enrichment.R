# R/enrichment.R
# Pathway enrichment on ExWAS results

#' @include AllClasses.R
#' @importFrom stats setNames
NULL

#' Gene set enrichment on sc-ExWAS results
#'
#' Runs fast gene set enrichment analysis (fgsea) on the
#' ranked gene list from \code{\link{run_sc_exwas}}, separately
#' for each cell type. This bridges the gap between
#' gene-level associations and biological interpretation --
#' a step that every scRNA-seq exposure paper performs
#' manually.
#'
#' @param exwas_result A \code{DataFrame} returned by
#'   \code{\link{run_sc_exwas}} or \code{\link{run_multi_exwas}}.
#' @param gene_sets A named list of character vectors (gene
#'   sets). Each element is a pathway name mapping to gene
#'   symbols. Compatible with MSigDB collections from
#'   \code{msigdbr}.
#' @param rank_by Character. Column used to rank genes.
#'   Default \code{"statistic"}, whose scale depends on the sc-ExWAS backend:
#'   signed square-root QL F for edgeR, Wald statistic for DESeq2, and
#'   moderated t statistic for voom-dream. The current \code{"log2FC"} field
#'   and legacy aliases \code{"stat"} and \code{"log2FoldChange"} are also
#'   supported.
#' @param min_size Integer. Minimum gene set size. Default 10.
#' @param max_size Integer. Maximum gene set size. Default 500.
#'
#' @return A \code{DataFrame} with columns: pathway, celltype,
#'   pval, padj, NES, size, exposure, leadingEdge (semicolon-
#'   separated).
#'
#' @details
#' This function requires the \code{fgsea} package
#' (Bioconductor). Gene sets can be obtained from MSigDB via
#' \code{msigdbr::msigdbr()} or from any named list.
#'
#' Genes are ranked by the backend-specific association statistic (default)
#' within each cell type, then fgsea is run per cell type. For edgeR QL this is
#' the signed square root of the QL F statistic and is not a Wald z statistic.
#' Results are combined with a final cross-cell-type FDR adjustment.
#'
#' This directly addresses a common workflow gap: researchers
#' routinely extract ExWAS hit genes, then manually run GO/
#' KEGG enrichment in a separate tool. By integrating
#' enrichment into the SCEE pipeline, the ranking preserves direction and
#' within-backend evidence ordering while FDR is handled across cell types.
#' Statistic magnitudes are not calibrated for comparison across different
#' backends; use \code{rank_by = "log2FC"} when an effect-size ranking is the
#' scientific target.
#'
#' @references
#' Korotkevich G et al. (2021). Fast gene set enrichment
#' analysis. \emph{bioRxiv}. \doi{10.1101/060012}
#'
#' Subramanian A et al. (2005). Gene set enrichment analysis.
#' \emph{PNAS} 102:15545-15550.
#' \doi{10.1073/pnas.0506580102}
#'
#' @export
#' @examples
#' # Requires fgsea package and gene sets
#' # See vignette for full example with MSigDB
#' showClass("SingleCellExposomeExperiment")
run_gsea <- function(exwas_result, gene_sets,
                      rank_by = "statistic",
                      min_size = 10L,
                      max_size = 500L) {

    if (!requireNamespace("fgsea", quietly = TRUE))
        stop("Package 'fgsea' required. ",
             "BiocManager::install('fgsea')")

    stopifnot(is.data.frame(exwas_result) ||
              is(exwas_result, "DataFrame"))
    stopifnot(is.list(gene_sets))

    res <- as.data.frame(exwas_result)

    rank_candidates <- switch(
        rank_by,
        statistic = c("statistic", "stat"),
        stat = c("stat", "statistic"),
        log2FC = c("log2FC", "log2FoldChange"),
        log2FoldChange = c("log2FoldChange", "log2FC"),
        rank_by
    )
    rank_column <- .resolve_exwas_result_column(
        res,
        rank_candidates,
        "gene-ranking"
    )

    required <- c("gene", "celltype")
    missing <- setdiff(required, colnames(res))
    if (length(missing)) {
        stop(
            "exwas_result is missing required column(s): ",
            paste(missing, collapse = ", "),
            "."
        )
    }
    if (!"exposure" %in% colnames(res)) {
        res$exposure <- NA_character_
    }
    strata <- unique(res[, c("celltype", "exposure"), drop = FALSE])
    all_gsea <- list()

    for (stratum_index in seq_len(nrow(strata))) {
        ct <- strata$celltype[[stratum_index]]
        exposure <- strata$exposure[[stratum_index]]
        same_exposure <- if (is.na(exposure)) {
            is.na(res$exposure)
        } else {
            !is.na(res$exposure) & res$exposure == exposure
        }
        ct_res <- res[
            !is.na(res$celltype) & res$celltype == ct & same_exposure,
            ,
            drop = FALSE
        ]
        ct_res <- ct_res[!is.na(ct_res[[rank_column]]), ]

        if (nrow(ct_res) < 10) next

        ## Build ranked gene vector
        ranks <- setNames(ct_res[[rank_column]], ct_res$gene)
        ranks <- sort(ranks, decreasing = TRUE)
        ## Remove duplicates (keep first = highest rank)
        ranks <- ranks[!duplicated(names(ranks))]

        gsea_out <- tryCatch(
            fgsea::fgsea(pathways = gene_sets,
                         stats = ranks,
                         minSize = min_size,
                         maxSize = max_size),
            error = function(e) {
                warning("fgsea failed for ", ct, ": ",
                        conditionMessage(e), call. = FALSE)
                NULL
            })

        if (is.null(gsea_out) || nrow(gsea_out) == 0) next

        gsea_df <- as.data.frame(gsea_out)
        gsea_df$celltype <- ct
        gsea_df$exposure <- exposure

        ## Convert leadingEdge list to semicolon-separated
        gsea_df$leadingEdge <- vapply(
            gsea_out$leadingEdge,
            function(x) paste(x, collapse = ";"),
            character(1))

        all_gsea <- c(all_gsea, list(gsea_df))
    }

    if (length(all_gsea) == 0) {
        return(S4Vectors::DataFrame(
            pathway = character(),
            celltype = character(),
            pval = numeric(),
            padj = numeric(),
            NES = numeric(),
            size = integer(),
            exposure = character(),
            leadingEdge = character()))
    }

    out <- do.call(rbind, all_gsea)
    ## Global FDR across cell types
    out$padj_global <- stats::p.adjust(out$pval, "BH")
    S4Vectors::DataFrame(out)
}
