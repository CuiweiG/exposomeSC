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
#'   p_interaction (F-test for the exposure-by-cell-type terms),
#'   df_interaction (numerator degrees of freedom),
#'   padj_interaction (BH-corrected), n_celltypes, n_donors.
#'
#' @details
#' For each gene, pseudobulk log-CPM values of all cell types are stacked,
#' one row per donor and cell type, and two linear models are compared with
#' an F-test:
#' \preformatted{
#' reduced: y ~ donor + celltype + covariates:celltype
#' full:    y ~ donor + celltype + covariates:celltype + exposure:celltype
#' }
#' Donor fixed effects absorb every donor-level quantity, including the
#' exposure and covariate main effects, and the shared donor component of
#' expression in different cell types, so the test compares exposure slopes
#' between cell types through within-donor contrasts rather than treating
#' cell types from the same donor as independent. Covariate-by-cell-type
#' terms allow covariate effects to differ between cell types. Donors
#' observed in only one cell type do not inform the interaction.
#'
#' Comparing significance levels across cell types is not the same as
#' testing for a difference between them (Gelman & Stern 2006,
#' Nieuwenhuis et al. 2011); this function tests the difference directly.
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
    covariates <- as.character(covariates)
    missing_covariates <- setdiff(covariates, colnames(exp_data))
    if (length(missing_covariates))
        stop("Covariate(s) not found in exposureData: ",
             paste(missing_covariates, collapse = ", "))

    cd <- SummarizedExperiment::colData(x)
    counts_mat <- SummarizedExperiment::assay(x, "counts")
    samples <- as.character(cd[[sample_col]])
    cell_types <- as.character(cd[[celltype_col]])
    all_cts <- sort(unique(cell_types))

    if (length(all_cts) < 2L)
        stop("Need >= 2 cell types for interaction test")

    ## Pseudobulk each cell type; keep the design separate from the
    ## expression values so that gene names cannot collide with it
    design_rows <- list()
    expression_blocks <- list()
    for (ct in all_cts) {
        pb_result <- .pseudobulk_aggregate(
            counts_mat, samples, cell_types, ct,
            min_cells = min_cells)
        if (is.null(pb_result)) next
        if (length(pb_result$valid_donors) < min_donors) next

        valid <- pb_result$valid_donors
        design_rows[[ct]] <- data.frame(
            .donor = valid, .celltype = ct, stringsAsFactors = FALSE)
        expression_blocks[[ct]] <- .log_cpm(pb_result$pb_mat)[, valid,
                                                              drop = FALSE]
    }
    if (length(design_rows) < 2L)
        stop("Need >= 2 cell types with enough donors for the ",
             "interaction test")

    design <- do.call(rbind, design_rows)
    genes <- Reduce(intersect, lapply(expression_blocks, rownames))
    expression <- do.call(cbind, lapply(expression_blocks, function(block)
        block[genes, , drop = FALSE]))
    design$.exposure <- as.numeric(exp_data[design$.donor, exposure])
    covariate_columns <- if (length(covariates)) {
        paste0(".covariate", seq_along(covariates))
    } else {
        character()
    }
    for (k in seq_along(covariates)) {
        design[[covariate_columns[k]]] <-
            as.numeric(exp_data[design$.donor, covariates[k]])
    }
    complete <- stats::complete.cases(design)
    design <- design[complete, , drop = FALSE]
    expression <- expression[, complete, drop = FALSE]
    design$.donor <- factor(design$.donor)
    design$.celltype <- factor(design$.celltype)
    if (nrow(design) < 4L)
        stop("Not enough valid donor-celltype strata")

    if (is.null(target_genes)) {
        gene_vars <- apply(expression, 1, stats::var, na.rm = TRUE)
        n_top <- min(50L, length(genes))
        target_genes <- names(sort(gene_vars,
            decreasing = TRUE))[seq_len(n_top)]
    }
    target_genes <- intersect(target_genes, genes)
    if (length(target_genes) == 0)
        stop("No target genes found")

    n_cts <- length(unique(design$.celltype))
    n_donors <- length(unique(design$.donor))
    reduced_terms <- c(".donor", ".celltype",
                       if (length(covariate_columns))
                           paste0(covariate_columns, ":.celltype"))
    reduced_formula <- stats::reformulate(reduced_terms, response = ".y")
    full_formula <- stats::reformulate(
        c(reduced_terms, ".exposure:.celltype"), response = ".y")

    ## Fit reduced and full models per gene
    results <- lapply(target_genes, function(gene) {
        model_data <- design
        model_data$.y <- as.numeric(expression[gene, ])

        fit_full <- tryCatch(lm(full_formula, data = model_data),
                             error = function(e) NULL)
        fit_reduced <- tryCatch(lm(reduced_formula, data = model_data),
                                error = function(e) NULL)
        if (is.null(fit_full) || is.null(fit_reduced))
            return(NULL)

        f_test <- tryCatch(
            anova(fit_reduced, fit_full),
            error = function(e) NULL)
        if (is.null(f_test) || nrow(f_test) < 2 ||
                !is.finite(f_test[2, "Df"]) || f_test[2, "Df"] < 1)
            return(NULL)

        data.frame(
            gene = gene,
            p_interaction = f_test[2, "Pr(>F)"],
            df_interaction = as.integer(f_test[2, "Df"]),
            n_celltypes = n_cts,
            n_donors = n_donors,
            stringsAsFactors = FALSE)
    })

    out <- do.call(rbind, Filter(Negate(is.null), results))
    if (is.null(out) || nrow(out) == 0) {
        return(S4Vectors::DataFrame(
            gene = character(),
            p_interaction = numeric(),
            df_interaction = integer(),
            padj_interaction = numeric(),
            n_celltypes = integer(),
            n_donors = integer()))
    }
    out$padj_interaction <- p.adjust(out$p_interaction, "BH")
    S4Vectors::DataFrame(out)
}
