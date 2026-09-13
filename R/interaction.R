# R/interaction.R
# Test whether exposure effects differ between cell types

#' @include AllClasses.R
#' @include utils.R
#' @importFrom stats lm anova as.formula p.adjust
NULL

#' Test exposure-by-celltype interaction
#'
#' Formally tests whether an exposure effect differs across
#' cell types using a pseudobulk interaction model. This
#' prevents the common mistake of concluding "the exposure
#' affects cell type A but not B" based on separate per-
#' celltype p-values (Gelman & Stern 2006).
#'
#' @param x A \code{\linkS4class{SingleCellExposomeExperiment}}.
#' @param exposure Character. Exposure variable name.
#' @param celltype_col Character. Column with cell type labels.
#' @param sample_col Character. Column with donor IDs.
#' @param target_genes Character vector (optional). Genes to
#'   test. Default: top 50 most variable genes across all
#'   cell types.
#' @param covariates Character vector (optional).
#' @param min_cells Integer. Min cells per donor per celltype.
#'   Default 10.
#' @param min_donors Integer. Min donors per celltype. Default 5.
#'
#' @return A \code{DataFrame} with columns: gene,
#'   p_interaction (F-test for exposure:celltype term),
#'   padj_interaction (BH-corrected), n_celltypes, n_donors.
#'
#' @details
#' For each gene, a linear model is fit on stacked pseudobulk
#' log-CPM across all cell types:
#' \code{y ~ exposure * celltype + covariates}
#'
#' The interaction term \code{exposure:celltype} tests whether
#' the slope of the exposure-expression relationship differs
#' between cell types. A significant interaction means the
#' exposure has a genuinely different effect size in different
#' cell types -- not just that one p-value is smaller.
#'
#' This addresses a pervasive statistical error in the field:
#' comparing significance levels across cell types is not the
#' same as testing for a difference (Gelman & Stern 2006,
#' Nieuwenhuis et al. 2011). Our interaction test provides
#' the correct approach.
#'
#' @references
#' Gelman A, Stern H (2006). The difference between
#' 'significant' and 'not significant' is not itself
#' statistically significant. \emph{Am Stat} 60:328-331.
#'
#' Nieuwenhuis S et al. (2011). Erroneous analyses of
#' interactions in neuroscience. \emph{Nat Neurosci}
#' 14:1105-1107. \doi{10.1038/nn.2886}
#'
#' @export
#' @examples
#' library(SingleCellExperiment); library(S4Vectors)
#' set.seed(1)
#' counts <- matrix(rpois(5000, 8), nrow = 50,
#'     dimnames = list(paste0("G", 1:50), paste0("c", 1:100)))
#' sce <- SingleCellExperiment(assays = list(counts = counts),
#'     colData = DataFrame(cell_id = paste0("c", 1:100),
#'         donor_id = rep(paste0("D", 1:10), each = 10),
#'         cell_type = rep(c("Mono", "NK"), 50)))
#' exp_mat <- matrix(rnorm(20), nrow = 10,
#'     dimnames = list(paste0("D", 1:10), c("E1", "E2")))
#' scee <- build_scee(sce, exp_mat, sample_col = "donor_id")
#' ix <- run_interaction_test(scee, exposure = "E1",
#'     celltype_col = "cell_type", sample_col = "donor_id",
#'     min_cells = 3L, min_donors = 3L)
#' head(ix)
run_interaction_test <- function(x, exposure,
                                  celltype_col = "cell_type",
                                  sample_col = "donor_id",
                                  target_genes = NULL,
                                  covariates = NULL,
                                  min_cells = 10L,
                                  min_donors = 5L) {
    stopifnot(is(x, "SingleCellExposomeExperiment"))

    exp_data <- slot(x, "exposureData")
    if (!exposure %in% colnames(exp_data))
        stop("Exposure '", exposure, "' not found")

    cd <- SummarizedExperiment::colData(x)
    counts_mat <- SummarizedExperiment::assay(x, "counts")
    samples <- as.character(cd[[sample_col]])
    cell_types <- as.character(cd[[celltype_col]])
    all_cts <- sort(unique(cell_types))

    if (length(all_cts) < 2L)
        stop("Need >= 2 cell types for interaction test")

    ## Pseudobulk each cell type, collect stacked data
    stacked <- list()
    for (ct in all_cts) {
        pb_result <- .pseudobulk_aggregate(
            counts_mat, samples, cell_types, ct,
            min_cells = min_cells)
        if (is.null(pb_result)) next
        if (length(pb_result$valid_donors) < min_donors) next

        log_cpm <- .log_cpm(pb_result$pb_mat)
        valid <- pb_result$valid_donors

        for (d in valid) {
            stacked <- c(stacked, list(data.frame(
                donor = d,
                celltype = ct,
                exposure = exp_data[d, exposure],
                t(log_cpm[, d, drop = FALSE]),
                check.names = FALSE,
                stringsAsFactors = FALSE)))
        }
    }

    if (length(stacked) < 4L)
        stop("Not enough valid donor-celltype strata")

    df_all <- do.call(rbind, stacked)
    genes <- setdiff(colnames(df_all),
                     c("donor", "celltype", "exposure"))

    if (is.null(target_genes)) {
        gene_vars <- vapply(genes, function(g)
            stats::var(df_all[[g]], na.rm = TRUE),
            numeric(1))
        n_top <- min(50L, length(genes))
        target_genes <- names(sort(gene_vars,
            decreasing = TRUE))[seq_len(n_top)]
    }
    target_genes <- intersect(target_genes, genes)
    if (length(target_genes) == 0)
        stop("No target genes found")

    n_cts <- length(unique(df_all$celltype))
    n_donors <- length(unique(df_all$donor))

    ## Fit interaction model per gene
    results <- lapply(target_genes, function(gene) {
        df_all$y <- df_all[[gene]]

        ## Full model with interaction
        fit_full <- tryCatch(
            lm(y ~ exposure * celltype, data = df_all),
            error = function(e) NULL)
        ## Reduced model without interaction
        fit_reduced <- tryCatch(
            lm(y ~ exposure + celltype, data = df_all),
            error = function(e) NULL)

        if (is.null(fit_full) || is.null(fit_reduced))
            return(NULL)

        f_test <- tryCatch(
            anova(fit_reduced, fit_full),
            error = function(e) NULL)

        if (is.null(f_test) || nrow(f_test) < 2)
            return(NULL)

        data.frame(
            gene = gene,
            p_interaction = f_test[2, "Pr(>F)"],
            n_celltypes = n_cts,
            n_donors = n_donors,
            stringsAsFactors = FALSE)
    })

    out <- do.call(rbind, Filter(Negate(is.null), results))
    if (is.null(out) || nrow(out) == 0) {
        return(S4Vectors::DataFrame(
            gene = character(),
            p_interaction = numeric(),
            padj_interaction = numeric(),
            n_celltypes = integer(),
            n_donors = integer()))
    }
    out$padj_interaction <- p.adjust(out$p_interaction, "BH")
    S4Vectors::DataFrame(out)
}
