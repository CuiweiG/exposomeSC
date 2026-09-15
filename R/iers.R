# R/iers.R
# Integrated Exposure Response Score

#' @importFrom stats quantile
NULL

#' Integrated Exposure Response Score (IERS)
#'
#' Combines sc-ExWAS results with an optional exposure-response
#' decomposition and optional cell coupling results into a single
#' per-gene score that captures the multi-scale exposure response.
#'
#' @param exwas A \code{DataFrame} from \code{run_sc_exwas}.
#' @param erd Optional \code{data.frame} with columns \code{gene} and
#'   \code{pct_compositional}, the percentage of the gene's exposure
#'   effect attributed to cell-type composition. No exposomeSC function
#'   produces this table: \code{run_decomposed_exwas()} has been removed
#'   and the developmental interventional decomposition reports a
#'   different quantity, so it has to be supplied by the caller. If NULL,
#'   the directness component is skipped.
#' @param coupling A \code{data.frame} from
#'   \code{run_cell_coupling}. Optional; if NULL, coupling
#'   component is skipped.
#' @param celltype Character. Cell type to compute IERS for; requires
#'   a \code{celltype} column in \code{exwas}.
#' @param weights Non-negative numeric vector of length 3 with a
#'   positive sum: weights for sc-ExWAS effect, directness, and
#'   coupling rewiring. Default \code{c(1, 1, 1)} (equal weights).
#'
#' @return A \code{data.frame} with columns: gene,
#'   score_exwas, score_direct, score_coupling, IERS,
#'   IERS_rank, IERS_percentile.
#'
#' @details
#' The IERS for gene \eqn{g} is:
#' \deqn{\text{IERS}_g = w_1 \cdot s_1(g) + w_2 \cdot s_2(g)
#'   + w_3 \cdot s_3(g)}
#' where:
#' \describe{
#'   \item{\eqn{s_1}}{Absolute backend-specific association statistic from
#'     sc-ExWAS (signed square-root QL F for edgeR, Wald statistic for DESeq2,
#'     or moderated t statistic for voom-dream)}
#'   \item{\eqn{s_2}}{Directness: \eqn{1 - pct\_comp/100},
#'     measuring how much of the effect is direct
#'     transcriptional rather than compositional}
#'   \item{\eqn{s_3}}{Maximum coupling |beta1| across all
#'     proteins (rewiring magnitude)}
#' }
#' All components are rank-normalised to [0, 1] before
#' combining, ensuring equal scale contribution. The sc-ExWAS component is a
#' within-input evidence rank rather than a calibrated z-score; IERS values
#' should not be compared across analyses fitted with different backends.
#'
#' @examples
#' exwas_df <- S4Vectors::DataFrame(
#'     gene = paste0("Gene", 1:5),
#'     log2FC = c(0.8, -0.4, 0.2, 0, 1.1),
#'     statistic = c(4, -2, 1, 0, 5),
#'     padj = c(0.001, 0.04, 0.2, 0.8, 0.0001))
#' iers <- compute_iers(exwas_df)
#' head(iers)
#' @export
compute_iers <- function(exwas, erd = NULL, coupling = NULL,
                          celltype = NULL,
                          weights = c(1, 1, 1)) {

    exwas_df <- as.data.frame(exwas)
    if (!"gene" %in% colnames(exwas_df)) {
        stop("exwas must contain a gene column.")
    }
    statistic_column <- .resolve_exwas_result_column(
        exwas_df,
        c("statistic", "stat"),
        "association-statistic"
    )

    if (!is.numeric(weights) || length(weights) != 3L ||
            any(!is.finite(weights)) || any(weights < 0) ||
            sum(weights) <= 0) {
        stop("weights must be three non-negative numbers with a ",
             "positive sum.")
    }

    ## Filter to celltype if specified
    if (!is.null(celltype)) {
        if (!"celltype" %in% colnames(exwas_df))
            stop("celltype was supplied but exwas has no celltype column.")
        exwas_df <- exwas_df[exwas_df$celltype == celltype, ]
    }

    if (!nrow(exwas_df)) {
        return(data.frame(
            gene = character(),
            score_exwas = numeric(),
            score_direct = numeric(),
            score_coupling = numeric(),
            IERS = numeric(),
            IERS_rank = numeric(),
            IERS_percentile = numeric(),
            stringsAsFactors = FALSE
        ))
    }

    genes <- exwas_df$gene

    ## Component 1: absolute backend-specific sc-ExWAS association statistic
    s1 <- abs(exwas_df[[statistic_column]])
    s1[is.na(s1)] <- 0
    s1_rank <- rank(s1) / length(s1)  # normalize to [0,1]

    ## Component 2: Directness from ERD
    s2_rank <- rep(0.5, length(genes))  # neutral if no ERD
    if (!is.null(erd)) {
        erd_df <- as.data.frame(erd)
        pct_comp <- erd_df$pct_compositional[
            match(genes, erd_df$gene)]
        ## Directness = 1 - pct_comp/100 (higher = more direct)
        ## Handle >100% (Simpson's paradox): cap at 0
        directness <- pmax(1 - pct_comp / 100, 0)
        directness[is.na(directness)] <- 0.5  # neutral
        s2_rank <- rank(directness) / length(directness)
    }

    ## Component 3: Max coupling rewiring
    s3_rank <- rep(0, length(genes))  # zero if no coupling
    if (!is.null(coupling)) {
        coupling_df <- as.data.frame(coupling)
        ## For each gene, find max |beta1| across all proteins
        max_rewiring <- vapply(genes, function(g) {
            rows <- coupling_df[coupling_df$gene == g, ]
            if (nrow(rows) == 0) return(0)
            max(abs(rows$beta1), na.rm = TRUE)
        }, numeric(1))
        max_rewiring[!is.finite(max_rewiring)] <- 0
        s3_rank <- rank(max_rewiring) / length(max_rewiring)
    }

    w <- weights / sum(weights)

    ## Compute IERS
    iers <- w[1] * s1_rank + w[2] * s2_rank + w[3] * s3_rank

    out <- data.frame(
        gene = genes,
        score_exwas = s1_rank,
        score_direct = s2_rank,
        score_coupling = s3_rank,
        IERS = iers,
        stringsAsFactors = FALSE)
    out$IERS_rank <- rank(-out$IERS)
    out$IERS_percentile <- 100 * rank(out$IERS) / nrow(out)
    out <- out[order(out$IERS_rank), ]
    rownames(out) <- NULL
    out
}
