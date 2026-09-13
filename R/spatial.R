# R/spatial.R
# Spatial exposomics: integration with SpatialExperiment

#' @include AllClasses.R
#' @include utils.R
#' @importFrom stats lm coef p.adjust cor.test
NULL

#' Spatial exposure-expression association
#'
#' Tests whether exposure effects on gene expression vary
#' across spatial regions within a tissue section. Integrates
#' \code{SpatialExperiment} with donor-level exposures to
#' enable spatially-resolved exposome analysis.
#'
#' @param spe A \code{SpatialExperiment} object with spatial
#'   coordinates. Must have \code{spatialCoords()}.
#' @param exposure_matrix Numeric matrix. Rows = samples,
#'   columns = exposure variables.
#' @param exposure Character. Exposure variable name.
#' @param sample_col Character. Column in \code{colData} with
#'   donor/sample IDs.
#' @param region_col Character (optional). Column with spatial
#'   region annotations (e.g., tissue compartment). If
#'   \code{NULL}, spatial regions are defined by k-means on
#'   coordinates.
#' @param n_regions Integer. Number of spatial regions for
#'   k-means (only if \code{region_col = NULL}). Default 4.
#' @param target_genes Character vector (optional). Genes to
#'   test. Default: top 50 most spatially variable.
#' @param min_spots Integer. Minimum spots per region per
#'   donor. Default 10.
#'
#' @return A \code{DataFrame} with columns: gene, region,
#'   coefficient (exposure effect), se, pvalue, padj,
#'   n_donors, method (DESeq2 or lm).
#'
#' @details
#' The analysis proceeds per spatial region:
#' \enumerate{
#'   \item Spots are grouped by region (annotated or k-means)
#'   \item Within each region, spots are pseudobulked per donor
#'   \item Exposure-expression association is tested via DESeq2
#'     on pseudobulk counts (consistent with core sc-ExWAS);
#'     falls back to linear regression on log-CPM when DESeq2
#'     is unavailable or fewer than 5 donors
#'   \item Results across regions reveal whether exposure
#'     effects are spatially heterogeneous
#' }
#'
#' This enables questions like: "Does PM2.5 affect gene
#' expression differently in the airway epithelium vs
#' submucosal tissue?" -- a spatial dimension that
#' dissociated scRNA-seq cannot address.
#'
#' @references
#' Righelli D et al. (2022). SpatialExperiment: infrastructure
#' for spatially-resolved transcriptomics data in R using
#' Bioconductor. \emph{Bioinformatics} 38:5397-5401.
#' \doi{10.1093/bioinformatics/btac299}
#'
#' @export
#' @examples
#' # Requires SpatialExperiment package
#' # See vignette for usage with Visium data
#' cat("See vignette for spatial examples\n")
run_spatial_exwas <- function(spe, exposure_matrix, exposure,
                               sample_col = "sample_id",
                               region_col = NULL,
                               n_regions = 4L,
                               target_genes = NULL,
                               min_spots = 10L) {

    if (!requireNamespace("SpatialExperiment", quietly = TRUE))
        stop("Package 'SpatialExperiment' required. ",
             "BiocManager::install('SpatialExperiment')")

    stopifnot(is(spe, "SpatialExperiment"))
    stopifnot(is.matrix(exposure_matrix))

    if (!exposure %in% colnames(exposure_matrix))
        stop("Exposure '", exposure, "' not found in matrix")

    cd <- SummarizedExperiment::colData(spe)
    if (!sample_col %in% colnames(cd))
        stop("'", sample_col, "' not in colData")

    samples <- as.character(cd[[sample_col]])
    donors <- unique(samples)
    exp_donors <- intersect(donors, rownames(exposure_matrix))

    if (length(exp_donors) < 3)
        stop("Need >= 3 donors with exposure data")

    ## Assign spatial regions
    if (is.null(region_col)) {
        coords <- SpatialExperiment::spatialCoords(spe)
        if (is.null(coords) || nrow(coords) == 0)
            stop("No spatial coordinates found")
        ## K-means clustering on spatial coordinates
        km <- stats::kmeans(coords, centers = n_regions,
            nstart = 10)
        regions <- paste0("Region_", km$cluster)
    } else {
        if (!region_col %in% colnames(cd))
            stop("'", region_col, "' not in colData")
        regions <- as.character(cd[[region_col]])
    }

    counts_mat <- SummarizedExperiment::assay(spe, "counts")
    region_levels <- sort(unique(regions))

    ## Gene selection: top spatially variable
    if (is.null(target_genes)) {
        ## Use variance across regions as proxy
        gene_vars <- apply(counts_mat, 1, stats::var)
        n_top <- min(50L, nrow(counts_mat))
        target_genes <- names(sort(gene_vars,
            decreasing = TRUE))[seq_len(n_top)]
    }
    target_genes <- intersect(target_genes, rownames(counts_mat))
    if (length(target_genes) == 0)
        stop("No target genes found")

    all_results <- list()

    for (reg in region_levels) {
        reg_mask <- regions == reg
        reg_samples <- samples[reg_mask]
        reg_counts <- counts_mat[target_genes, reg_mask,
            drop = FALSE]

        ## Pseudobulk per donor within region
        donor_counts <- lapply(exp_donors, function(d) {
            d_mask <- reg_samples == d
            if (sum(d_mask) < min_spots) return(NULL)
            sub <- reg_counts[, d_mask, drop = FALSE]
            if (methods::is(sub, "sparseMatrix")) {
                Matrix::rowSums(sub)
            } else if (sum(d_mask) == 1L) {
                as.numeric(sub)
            } else {
                rowSums(sub)
            }
        })
        names(donor_counts) <- exp_donors

        ## Remove NULL entries (donors with < min_spots)
        keep <- !vapply(donor_counts, is.null, logical(1))
        if (sum(keep) < 3) next

        valid <- exp_donors[keep]
        pb_mat <- do.call(cbind, donor_counts[keep])
        colnames(pb_mat) <- valid

        ## Log-CPM
        lib <- colSums(pb_mat)
        if (any(lib == 0)) next
        log_cpm <- log2(t(t(pb_mat) / lib * 1e6) + 1)

        dose <- exposure_matrix[valid, exposure]

        ## Use DESeq2 if available (consistent with core
        ## sc-ExWAS), else fall back to lm on log-CPM
        has_deseq2 <- requireNamespace("DESeq2",
            quietly = TRUE)

        if (has_deseq2 && length(valid) >= 5L) {
            ## DESeq2 on pseudobulk counts (consistent
            ## with run_sc_exwas methodology)
            design_df <- data.frame(
                sample = valid,
                exposure = dose,
                stringsAsFactors = FALSE)
            rownames(design_df) <- valid

            ## Filter genes: expressed in >= 20% of donors
            min_samp <- max(2L,
                ceiling(0.2 * ncol(pb_mat)))
            keep_g <- rowSums(pb_mat > 0) >= min_samp
            pb_filt <- pb_mat[keep_g, , drop = FALSE]

            if (nrow(pb_filt) > 0) {
                reg_df <- tryCatch({
                    dds <- DESeq2::DESeqDataSetFromMatrix(
                        countData = pb_filt,
                        colData = design_df,
                        design = ~ exposure)
                    dds <- DESeq2::DESeq(dds, quiet = TRUE)
                    res <- DESeq2::results(dds,
                        name = "exposure")
                    data.frame(
                        gene = rownames(res),
                        region = reg,
                        coefficient = res$log2FoldChange,
                        se = res$lfcSE,
                        pvalue = res$pvalue,
                        n_donors = length(valid),
                        method = "DESeq2",
                        stringsAsFactors = FALSE)
                }, error = function(e) NULL)
            } else {
                reg_df <- NULL
            }
        } else {
            ## Fallback: lm on log-CPM
            reg_res <- lapply(target_genes, function(gene) {
                y <- log_cpm[gene, ]
                fit <- tryCatch(
                    summary(lm(y ~ dose)),
                    error = function(e) NULL)
                if (is.null(fit)) return(NULL)
                cf <- fit$coefficients
                if (!"dose" %in% rownames(cf)) return(NULL)
                data.frame(
                    gene = gene,
                    region = reg,
                    coefficient = cf["dose", "Estimate"],
                    se = cf["dose", "Std. Error"],
                    pvalue = cf["dose", "Pr(>|t|)"],
                    n_donors = length(valid),
                    method = "lm",
                    stringsAsFactors = FALSE)
            })
            reg_df <- do.call(rbind,
                Filter(Negate(is.null), reg_res))
        }

        if (!is.null(reg_df) && nrow(reg_df) > 0) {
            reg_df$padj <- p.adjust(reg_df$pvalue, "BH")
            all_results <- c(all_results, list(reg_df))
        }
    }

    if (length(all_results) == 0) {
        return(S4Vectors::DataFrame(
            gene = character(), region = character(),
            coefficient = numeric(), se = numeric(),
            pvalue = numeric(), padj = numeric(),
            n_donors = integer(), method = character()))
    }

    out <- do.call(rbind, all_results)
    out$padj_global <- p.adjust(out$pvalue, "BH")
    S4Vectors::DataFrame(out)
}
