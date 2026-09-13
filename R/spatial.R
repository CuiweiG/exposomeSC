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
#'   coordinates, which requires a single tissue section or
#'   \code{coordinates_registered = TRUE}.
#' @param n_regions Integer. Number of spatial regions for
#'   k-means (only if \code{region_col = NULL}). Default 4.
#' @param target_genes Character vector (optional). Genes to
#'   report. Default: the 50 genes with the largest variance of
#'   spot counts.
#' @param min_spots Integer. Minimum spots per region per
#'   donor. Default 10.
#' @param coordinates_registered Logical. Set to \code{TRUE} only
#'   when the spatial coordinates of all sections share one
#'   registered frame, so that k-means regions correspond between
#'   sections. Default \code{FALSE}.
#' @param seed Integer seed for k-means. The caller's random
#'   number generator state is restored on return. Default 1.
#'
#' @return A \code{DataFrame} with columns: gene, region,
#'   coefficient (exposure effect), se, pvalue, padj (within
#'   region), n_donors, method (DESeq2 or lm) and padj_global
#'   (Benjamini-Hochberg across regions).
#'
#' @details
#' The analysis proceeds per spatial region:
#' \enumerate{
#'   \item Spots are grouped by region (annotated or k-means)
#'   \item Within each region, all genes are pseudobulked per
#'     donor, so that normalisation uses the whole transcriptome
#'   \item Exposure-expression association is tested via DESeq2
#'     on pseudobulk counts and reported for the target genes;
#'     with fewer than 5 donors, or without DESeq2, linear
#'     regression on log-CPM is used instead
#'   \item Results across regions describe whether exposure
#'     effects differ between regions; comparing regions
#'     formally requires an interaction model
#' }
#'
#' Coordinates of separate tissue sections are generally not in
#' a common frame, so k-means regions computed from pooled
#' coordinates would not correspond between donors. Without
#' \code{region_col}, more than one section is therefore an
#' error unless \code{coordinates_registered = TRUE}.
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
#' if (requireNamespace("SpatialExperiment", quietly = TRUE)) {
#'     set.seed(1)
#'     donors <- paste0("D", 1:6)
#'     spot_donor <- rep(donors, each = 40)
#'     region <- rep(rep(c("epithelium", "stroma"), each = 20), times = 6)
#'     counts <- matrix(stats::rpois(200 * length(spot_donor), 5),
#'         nrow = 200, dimnames = list(paste0("G", 1:200),
#'             paste0("s", seq_along(spot_donor))))
#'     spe <- SpatialExperiment::SpatialExperiment(
#'         assays = list(counts = counts),
#'         colData = S4Vectors::DataFrame(sample_id = spot_donor,
#'             region = region),
#'         spatialCoords = cbind(x = stats::runif(240),
#'             y = stats::runif(240)))
#'     exposure <- matrix(stats::rnorm(6), ncol = 1,
#'         dimnames = list(donors, "PM2.5"))
#'     run_spatial_exwas(spe, exposure, "PM2.5", region_col = "region",
#'         target_genes = paste0("G", 1:10))
#' }
run_spatial_exwas <- function(spe, exposure_matrix, exposure,
                               sample_col = "sample_id",
                               region_col = NULL,
                               n_regions = 4L,
                               target_genes = NULL,
                               min_spots = 10L,
                               coordinates_registered = FALSE,
                               seed = 1L) {

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
    exp_donors <- exp_donors[!is.na(exposure_matrix[exp_donors, exposure])]

    if (length(exp_donors) < 3)
        stop("Need >= 3 donors with exposure data")

    ## Assign spatial regions
    if (is.null(region_col)) {
        if (length(donors) > 1L && !isTRUE(coordinates_registered))
            stop("k-means regions from coordinates pooled across tissue ",
                 "sections correspond between sections only when the ",
                 "coordinates share a registered frame. Supply region_col, ",
                 "or set coordinates_registered = TRUE.")
        coords <- SpatialExperiment::spatialCoords(spe)
        if (is.null(coords) || nrow(coords) == 0)
            stop("No spatial coordinates found")
        .local_rng_scope(seed)
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

    ## Default genes: largest variance of spot counts, computed without
    ## densifying a sparse matrix
    if (is.null(target_genes)) {
        n_spots <- ncol(counts_mat)
        row_mean <- Matrix::rowMeans(counts_mat)
        gene_vars <- (Matrix::rowSums(counts_mat^2) - n_spots * row_mean^2) /
            (n_spots - 1)
        names(gene_vars) <- rownames(counts_mat)
        n_top <- min(50L, nrow(counts_mat))
        target_genes <- names(sort(gene_vars,
            decreasing = TRUE))[seq_len(n_top)]
    }
    target_genes <- intersect(target_genes, rownames(counts_mat))
    if (length(target_genes) == 0)
        stop("No target genes found")
    has_deseq2 <- requireNamespace("DESeq2", quietly = TRUE)

    all_results <- list()

    for (reg in region_levels) {
        reg_mask <- regions == reg
        reg_samples <- samples[reg_mask]
        keep_donor <- vapply(exp_donors, function(d)
            sum(reg_samples == d) >= min_spots, logical(1))
        if (sum(keep_donor) < 3) next
        valid <- exp_donors[keep_donor]

        ## Pseudobulk every gene per donor within the region
        pb_mat <- vapply(valid, function(d) {
            as.numeric(Matrix::rowSums(
                counts_mat[, reg_mask & samples == d, drop = FALSE]))
        }, numeric(nrow(counts_mat)))
        pb_mat <- matrix(pb_mat, nrow = nrow(counts_mat),
            dimnames = list(rownames(counts_mat), valid))

        lib <- colSums(pb_mat)
        if (any(lib == 0)) next
        dose <- exposure_matrix[valid, exposure]

        reg_df <- NULL
        if (has_deseq2 && length(valid) >= 5L) {
            ## DESeq2 on all genes expressed in >= 20% of donors;
            ## results are reported for the target genes
            min_samp <- max(2L, ceiling(0.2 * ncol(pb_mat)))
            pb_filt <- pb_mat[rowSums(pb_mat > 0) >= min_samp, ,
                              drop = FALSE]
            tested <- intersect(target_genes, rownames(pb_filt))
            if (length(tested)) {
                reg_df <- tryCatch({
                    design_df <- data.frame(exposure = dose,
                                            row.names = valid)
                    dds <- DESeq2::DESeqDataSetFromMatrix(
                        countData = round(pb_filt),
                        colData = design_df,
                        design = ~ exposure)
                    dds <- DESeq2::DESeq(dds, quiet = TRUE)
                    res <- DESeq2::results(dds, name = "exposure")
                    res <- res[tested, , drop = FALSE]
                    data.frame(
                        gene = tested,
                        region = reg,
                        coefficient = res$log2FoldChange,
                        se = res$lfcSE,
                        pvalue = res$pvalue,
                        n_donors = length(valid),
                        method = "DESeq2",
                        stringsAsFactors = FALSE)
                }, error = function(e) {
                    warning("DESeq2 failed for region ", reg, ": ",
                            conditionMessage(e), call. = FALSE)
                    NULL
                })
            }
        } else {
            ## Fallback: lm on log-CPM from whole-transcriptome library
            ## sizes
            log_cpm <- log2(t(t(pb_mat) / lib * 1e6) + 1)
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
            n_donors = integer(), method = character(),
            padj_global = numeric()))
    }

    out <- do.call(rbind, all_results)
    out$padj_global <- p.adjust(out$pvalue, "BH")
    S4Vectors::DataFrame(out)
}
