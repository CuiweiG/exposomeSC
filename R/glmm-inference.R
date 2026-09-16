# R/glmm-inference.R
# Donor-level inference for sc-ExWAS, replacing the anticonservative GLMM
# path.
#
# The retired path fitted the mixed model with ML (lme4 REML = FALSE) and
# referred a Wald statistic to the normal (p = 2 * pnorm(|z|)). Both choices
# inflate type-I error for the small donor samples of an sc-ExWAS: ML biases
# the variance components downward, and the Wald-normal reference ignores the
# finite denominator degrees of freedom. This file replaces that path with two
# estimators whose calibration is defensible:
#
#   (a) run_sc_exwas_pb_offset() -- donor offset-pseudobulk negative-binomial
#       GLM. For a donor-level exposure (all of a donor's cells share A_i), the
#       donor NB GLM with offset log s_i = log(donor total counts) has point
#       estimate AND standard error equal to the cell-level NB-GLMM up to
#       O(1/sqrt(n)), and is an order of magnitude faster (Lee & Han 2024,
#       Bioinformatics 40(8):btae498). This is the recommended default for the
#       direct/marginal channel. It does NOT extend to within-donor-varying
#       terms.
#
#   (b) run_sc_exwas_glmm() -- cell-level linear mixed model
#       y ~ exposure * celltype + covariates + (1 | donor), fitted with
#       lmerTest::lmer(REML = TRUE). p-values come from a finite denominator df
#       -- Satterthwaite (default) or Kenward-Roger -- NOT the Wald normal. REML
#       is MANDATORY: the Kenward-Roger covariance and the Satterthwaite
#       Hessian are defined only on the REML fit.
#
# Calibration: do not quote a nominal type-I rate for either estimator without
# a donor-permutation check (permute the exposure ACROSS donors, refit the
# frozen pipeline, and report a Wilson interval for the empirical FPR).

#' @include AllClasses.R
#' @include AllGenerics.R
#' @include utils.R
#' @importFrom stats p.adjust as.formula coef vcov anova setNames
NULL

# ---- Parallel backend ---------------------------------------------------
# Both entry points run serially unless the caller asks for otherwise: the
# default BPPARAM is SerialParam(), so nothing here starts a worker process on
# its own. Examples, vignettes and tests must use at most two cores, and a
# default scaled to the host's core count would breach that on a large
# machine. A caller who wants parallel gene fits passes an explicit param,
# for example BiocParallel::SnowParam(workers = 2L).

# bplapply when a BiocParallel param is supplied, else a serial lapply.
.bpapply <- function(X, FUN, BPPARAM) {
    if (!is.null(BPPARAM) && requireNamespace("BiocParallel", quietly = TRUE)) {
        BiocParallel::bplapply(X, FUN, BPPARAM = BPPARAM)
    } else {
        lapply(X, FUN)
    }
}

# Fit the offset-NB GLM per gene for ONE cell type, parallelised here so the
# worker closure captures only this cell type's small pseudobulk (pb_mat) and
# NEVER the full single-cell counts matrix in the caller's frame.
.pb_fit_celltype <- function(pb_mat, A, logoff, covariates, X, genes, BPPARAM) {
    rhs <- "A + offset(logoff)"
    if (length(covariates))
        rhs <- paste(rhs, "+", paste(covariates, collapse = " + "))
    form <- stats::as.formula(paste("y ~", rhs))
    one_gene <- function(g) {
        d <- data.frame(y = as.numeric(pb_mat[g, ]), A = A, logoff = logoff)
        if (!is.null(X)) d <- cbind(d, X)
        fit <- tryCatch(MASS::glm.nb(form, data = d), error = function(e) NULL)
        if (is.null(fit))
            return(c(log2FC = NA_real_, se = NA_real_,
                     statistic = NA_real_, pvalue = NA_real_))
        s <- tryCatch(stats::coef(summary(fit)), error = function(e) NULL)
        if (is.null(s) || !"A" %in% rownames(s))
            return(c(log2FC = NA_real_, se = NA_real_,
                     statistic = NA_real_, pvalue = NA_real_))
        beta <- s["A", "Estimate"]; se_nat <- s["A", "Std. Error"]
        c(log2FC = beta / log(2), se = se_nat / log(2),
          statistic = s["A", "z value"], pvalue = s["A", "Pr(>|z|)"])
    }
    do.call(rbind, .bpapply(genes, one_gene, BPPARAM = BPPARAM))
}

# cell-level log2 counts-per-10k (+1); the LMM response for family = "gaussian".
# Library size uses the FULL count matrix (all genes), not the target subset.
.log_cp10k_cell <- function(counts_full, genes) {
    lib <- pmax(Matrix::colSums(counts_full), 1)
    sub <- as.matrix(counts_full[genes, , drop = FALSE])
    log2(t(t(sub) / lib * 1e4) + 1)
}

#' Donor offset-pseudobulk negative-binomial sc-ExWAS
#'
#' A donor-level negative-binomial alternative to \code{\link{run_sc_exwas}}.
#' Per cell type and per gene, aggregates cells to a donor pseudobulk count and
#' fits a negative-binomial GLM with a library-size offset
#' \eqn{\log s_i = \log(\text{donor total counts})}. Lee and Han (2024) report
#' that pseudobulk models with proper offsets have the same statistical
#' properties as generalised linear mixed models in single-cell case-control
#' studies, without the pseudoreplication that inflates cell-level tests.
#' P-values are Wald tests with an estimated dispersion and can be
#' anticonservative when there are few donors.
#'
#' @param scee A \linkS4class{SingleCellExposomeExperiment}.
#' @param exposure Character; a donor-level column of \code{exposureData(scee)}.
#' @param celltype_col Character; \code{colData} column with cell-type labels.
#' @param sample_col Character; donor id column (default \code{"donor_id"}).
#' @param covariates Character vector; donor covariate columns in
#'   \code{exposureData} (numeric), added as GLM main effects.
#' @param target_genes Character; genes to test (default: all rows of
#'   \code{scee}).
#' @param BPPARAM A \code{BiocParallel} param used to parallelise over genes.
#'   Defaults to \code{BiocParallel::SerialParam()}, so no worker process is
#'   started unless one is asked for; pass, for example,
#'   \code{BiocParallel::SnowParam(workers = 2L)} to fit genes in parallel.
#'   Pass \code{NULL} to fall back to a plain \code{lapply}.
#'
#' @return A \code{data.frame}, one row per gene x cell type:
#'   \code{gene}, \code{celltype}, \code{log2FC} (exposure effect on the
#'   \eqn{\log_2} rate scale), \code{se} (on the same scale), \code{statistic}
#'   (donor-level Wald z from the NB GLM), \code{pvalue}, and BH-adjusted
#'   \code{padj} across all strata. Rows ordered by \code{padj}.
#'
#' @note Report empirical calibration from a donor-permutation null, with the
#'   exposure permuted across donors, before quoting any type-I rate.
#'
#' @references
#' Lee H & Han B (2024). Pseudobulk with proper offsets has the same statistical
#' properties as generalized linear mixed models in single-cell case-control
#' studies. \emph{Bioinformatics} 40(8):btae498.
#'
#' Squair JW et al. (2021). Confronting false discoveries in single-cell
#' differential expression. \emph{Nat Commun} 12:5692.
#'
#' @seealso \code{\link{run_sc_exwas_glmm}} for the cell-level LMM reference.
#' @examples
#' donor_ids <- paste0("D", seq_len(10L))
#' exposure <- rep(c(0, 1), each = 5L)
#' cell_donor <- rep(donor_ids, each = 10L)
#' cell_type <- rep("TypeA", length(cell_donor))
#' cell_ids <- paste0("cell", seq_along(cell_donor))
#' pb_targets <- rbind(
#'     G1 = c(30, 140, 50, 180, 40, 110, 260, 70, 220, 130),
#'     G2 = c(200, 40, 240, 60, 300, 50, 270, 80, 330, 70),
#'     G3 = c(100, 120, 90, 150, 110, 130, 170, 100, 160, 140),
#'     G4 = c(250, 200, 300, 180, 270, 220, 320, 210, 290, 230)
#' )
#' counts <- matrix(0L, nrow = 4L, ncol = length(cell_ids))
#' for (donor in seq_along(donor_ids)) {
#'     counts[, cell_donor == donor_ids[donor]] <- pb_targets[, donor] / 10L
#' }
#' storage.mode(counts) <- "integer"
#' dimnames(counts) <- list(paste0("G", seq_len(4L)), cell_ids)
#' sce <- SingleCellExperiment::SingleCellExperiment(
#'     assays = list(counts = counts),
#'     colData = S4Vectors::DataFrame(
#'         cell_id = cell_ids,
#'         donor_id = cell_donor,
#'         cell_type = cell_type
#'     )
#' )
#' donor_design <- matrix(exposure, ncol = 1L,
#'     dimnames = list(donor_ids, "exposure"))
#' scee <- build_scee(sce, donor_design, sample_col = "donor_id")
#' offset_result <- run_sc_exwas_pb_offset(
#'     scee,
#'     exposure = "exposure",
#'     celltype_col = "cell_type",
#'     target_genes = c("G1", "G2"),
#'     BPPARAM = BiocParallel::SerialParam()
#' )
#' head(offset_result)
#' @export
run_sc_exwas_pb_offset <- function(scee, exposure, celltype_col,
                                    sample_col = "donor_id",
                                    covariates = NULL, target_genes = NULL,
                                    BPPARAM = BiocParallel::SerialParam()) {
    stopifnot(methods::is(scee, "SingleCellExposomeExperiment"))
    if (!requireNamespace("MASS", quietly = TRUE))
        stop("Package 'MASS' is required for the negative-binomial GLM")

    exp_data <- exposureData(scee)
    if (!exposure %in% colnames(exp_data))
        stop("Exposure '", exposure, "' not found in exposureData(scee)")
    if (length(covariates) && !all(covariates %in% colnames(exp_data)))
        stop("Covariate(s) not in exposureData: ",
             paste(setdiff(covariates, colnames(exp_data)), collapse = ", "))

    cd <- SummarizedExperiment::colData(scee)
    counts_mat <- SummarizedExperiment::assay(scee, "counts")
    samples <- as.character(cd[[sample_col]])
    cell_types <- as.character(cd[[celltype_col]])
    celltypes <- sort(unique(cell_types))

    ## The donor offset absorbs library-size differences, so we keep the
    ## per-donor cell threshold permissive (1) rather than silently dropping
    ## donors; the NB likelihood down-weights sparse pseudobulk columns.
    min_cells <- 1L

    per_ct <- lapply(celltypes, function(ct) {
        agg <- .pseudobulk_aggregate(counts_mat, samples, cell_types, ct,
                                     min_cells = min_cells)
        if (is.null(agg)) return(NULL)
        donors <- agg$valid_donors
        if (length(donors) < 3L) {
            message("Skipping cell type '", ct, "': < 3 donors with cells")
            return(NULL)
        }
        pb_mat <- agg$pb_mat
        logoff <- log(pmax(colSums(pb_mat), 1))  # log donor total counts
        A <- as.numeric(exp_data[donors, exposure])
        X <- if (length(covariates))
            as.data.frame(exp_data[donors, covariates, drop = FALSE]) else NULL
        genes <- if (is.null(target_genes)) rownames(pb_mat)
                 else intersect(target_genes, rownames(pb_mat))
        if (!length(genes)) return(NULL)

        est <- .pb_fit_celltype(pb_mat, A, logoff, covariates, X, genes,
                                BPPARAM)
        data.frame(gene = genes, celltype = ct,
                   log2FC = est[, "log2FC"], se = est[, "se"],
                   statistic = est[, "statistic"], pvalue = est[, "pvalue"],
                   row.names = NULL, stringsAsFactors = FALSE)
    })

    out <- do.call(rbind, Filter(Negate(is.null), per_ct))
    if (is.null(out) || !nrow(out))
        stop("No gene/cell-type stratum yielded an estimable NB GLM")
    out$padj <- stats::p.adjust(out$pvalue, "BH")
    out[order(out$padj, out$pvalue), , drop = FALSE]
}

#' Cell-level linear mixed model sc-ExWAS with finite-df inference
#'
#' Fits, per gene, the cell-level model
#' \code{y ~ exposure * celltype + covariates + (1 | donor) +
#' (1 | donor:celltype)} and reports
#' cell-type-specific exposure effects. For \code{family = "gaussian"} the
#' model is a linear mixed model fitted by \code{lmerTest::lmer} with
#' \strong{REML = TRUE} (mandatory), and p-values are referred to a finite
#' denominator degrees of freedom -- \strong{Satterthwaite} (default) or
#' \strong{Kenward-Roger} -- via \code{lmerTest::contest1D}. This replaces the
#' anticonservative Wald-normal test of the retired \code{glmm.R} path.
#'
#' @details
#' \strong{Why REML.} Maximum likelihood biases the variance components
#' downward, which shrinks the residual scale and inflates the exposure
#' \eqn{t}-statistic; the Kenward-Roger adjusted covariance and the
#' Satterthwaite degrees-of-freedom approximation are both defined only on the
#' REML fit. The function therefore forces \code{REML = TRUE} and does not
#' expose a switch.
#'
#' \strong{Random effects.} The exposure varies between donors, so the
#' replicates for a cell-type-specific exposure effect are donor-by-cell-type
#' units, not cells. The model therefore includes a donor-by-cell-type random
#' intercept as well as a donor random intercept; without it, variation shared
#' by the cells of one donor and cell type is treated as independent and the
#' cell-type-specific tests are anticonservative. With a single cell type the
#' model reduces to \code{y ~ exposure + covariates + (1 | donor)}. A variance
#' component estimated at zero (a singular fit) is expected when that source
#' of variation is absent and is not reported.
#'
#' \strong{Cell-type-specific effects.} Under treatment contrasts the
#' \code{exposure} coefficient is the effect at the reference cell type and each
#' \code{exposure:celltype} coefficient is a difference from it. Each cell
#' type's effect is therefore a linear combination of fixed effects, tested with
#' \code{contest1D} so that the degrees of freedom and standard error account
#' for the covariance between the main effect and the interaction.
#'
#' \strong{family = "nbinom".} Satterthwaite and Kenward-Roger are linear
#' -mixed-model constructions and do not extend to a generalised LMM, so the
#' negative-binomial path fits \code{lme4::glmer.nb} on an equivalent
#' cell-type-nested reparametrisation of \code{exposure * celltype} (one
#' exposure slope column per cell type) and tests each cell type's effect by a
#' likelihood-ratio test (\code{anova} of the full versus the column-dropped
#' model) -- again not a Wald-normal test. This path is markedly more expensive
#' (one \code{glmer.nb} refit per cell type per gene); for count data prefer the
#' donor-level \code{\link{run_sc_exwas_pb_offset}}, which is calibration-
#' equivalent and far faster (Lee & Han 2024). The \code{ddf} argument is
#' ignored for \code{family = "nbinom"}.
#'
#' @param scee A \linkS4class{SingleCellExposomeExperiment}.
#' @param exposure Character; a donor-level column of \code{exposureData(scee)}.
#' @param celltype_col Character; \code{colData} column with cell-type labels.
#' @param sample_col Character; donor id column (default \code{"donor_id"}).
#' @param covariates Character vector; donor covariate columns in
#'   \code{exposureData} (numeric), added as fixed-effect main effects.
#' @param target_genes Character; genes to test (default: all rows of
#'   \code{scee}). A cell-level LMM per gene is costly -- passing a focused gene
#'   set is strongly recommended.
#' @param ddf Denominator degrees-of-freedom method for
#'   \code{family = "gaussian"}: \code{"Satterthwaite"} (default) or
#'   \code{"Kenward-Roger"} (needs \pkg{pbkrtest}; use for \eqn{n \lesssim 30}
#'   donors).
#' @param family \code{"gaussian"} (LMM on cell-level \eqn{\log_2} CP10k,
#'   default) or \code{"nbinom"} (negative-binomial GLMM on cell counts).
#' @param BPPARAM A \code{BiocParallel} param used to parallelise over genes.
#'   Defaults to \code{BiocParallel::SerialParam()}, so no worker process is
#'   started unless one is asked for; pass, for example,
#'   \code{BiocParallel::SnowParam(workers = 2L)} to fit genes in parallel.
#'
#' @return A \code{data.frame}, one row per gene x cell type: \code{gene},
#'   \code{celltype}, \code{log2FC}, \code{se}, \code{statistic} (Satterthwaite
#'   /KR \eqn{t} for gaussian, LR \eqn{\chi^2} for nbinom), \code{df} (finite
#'   denominator df; 1 for the nbinom LRT), \code{pvalue}, BH-adjusted
#'   \code{padj}, and \code{method}. Rows ordered by \code{padj}.
#'
#' @note Calibrate against a donor-permutation null, with the exposure permuted
#'   across donors, before quoting any type-I rate;
#'   Satterthwaite/KR control the finite-sample error only under the working
#'   model, so the empirical check is the headline.
#'
#' @references
#' Kuznetsova A, Brockhoff PB & Christensen RHB (2017). lmerTest package: tests
#' in linear mixed effects models. \emph{J Stat Softw} 82(13).
#'
#' Kenward MG & Roger JH (1997). Small sample inference for fixed effects from
#' restricted maximum likelihood. \emph{Biometrics} 53:983-997.
#'
#' Halekoh U & Hojsgaard S (2014). A Kenward-Roger approximation and parametric
#' bootstrap methods for tests in linear mixed models -- the R package pbkrtest.
#' \emph{J Stat Softw} 59(9).
#'
#' @seealso \code{\link{run_sc_exwas_pb_offset}} for the donor-level default.
#' @examples
#' donor_ids <- paste0("D", seq_len(10L))
#' exposure <- rep(c(0, 1), each = 5L)
#' cell_donor <- rep(donor_ids, each = 12L)
#' cell_type <- rep(
#'     rep(c("TypeA", "TypeB"), each = 6L),
#'     times = length(donor_ids)
#' )
#' cell_ids <- paste0("cell", seq_along(cell_donor))
#' counts <- outer(
#'     seq_len(2L),
#'     seq_along(cell_ids),
#'     function(gene, cell) 2L + ((7L * gene + 5L * cell) %% 11L)
#' )
#' exposed_type_a <- exposure[match(cell_donor, donor_ids)] == 1 &
#'     cell_type == "TypeA"
#' counts[1L, exposed_type_a] <- counts[1L, exposed_type_a] + 3L
#' storage.mode(counts) <- "integer"
#' dimnames(counts) <- list(c("G1", "G2"), cell_ids)
#' sce <- SingleCellExperiment::SingleCellExperiment(
#'     assays = list(counts = counts),
#'     colData = S4Vectors::DataFrame(
#'         cell_id = cell_ids,
#'         donor_id = cell_donor,
#'         cell_type = cell_type
#'     )
#' )
#' donor_design <- matrix(exposure, ncol = 1L,
#'     dimnames = list(donor_ids, "exposure"))
#' scee <- build_scee(sce, donor_design, sample_col = "donor_id")
#' if (requireNamespace("lme4", quietly = TRUE) &&
#'         requireNamespace("lmerTest", quietly = TRUE)) {
#'     mixed_result <- run_sc_exwas_glmm(
#'         scee,
#'         exposure = "exposure",
#'         celltype_col = "cell_type",
#'         target_genes = "G1",
#'         family = "gaussian",
#'         BPPARAM = BiocParallel::SerialParam()
#'     )
#'     mixed_result
#' }
#' @export
run_sc_exwas_glmm <- function(scee, exposure, celltype_col,
                              sample_col = "donor_id",
                              covariates = NULL, target_genes = NULL,
                              ddf = c("Satterthwaite", "Kenward-Roger"),
                              family = c("gaussian", "nbinom"),
                              BPPARAM = BiocParallel::SerialParam()) {
    stopifnot(methods::is(scee, "SingleCellExposomeExperiment"))
    ddf <- match.arg(ddf); family <- match.arg(family)
    if (!requireNamespace("lmerTest", quietly = TRUE))
        stop("Package 'lmerTest' is required (REML LMM + Satterthwaite/KR df)")
    if (!requireNamespace("lme4", quietly = TRUE))
        stop("Package 'lme4' is required")
    if (family == "gaussian" && ddf == "Kenward-Roger" &&
        !requireNamespace("pbkrtest", quietly = TRUE))
        stop("Kenward-Roger df require the 'pbkrtest' package")
    if (family == "nbinom" && !requireNamespace("MASS", quietly = TRUE))
        stop("Package 'MASS' is required for lme4::glmer.nb (theta estimation)")

    exp_data <- exposureData(scee)
    if (!exposure %in% colnames(exp_data))
        stop("Exposure '", exposure, "' not found in exposureData(scee)")
    if (length(covariates) && !all(covariates %in% colnames(exp_data)))
        stop("Covariate(s) not in exposureData: ",
             paste(setdiff(covariates, colnames(exp_data)), collapse = ", "))

    cd <- SummarizedExperiment::colData(scee)
    counts_mat <- SummarizedExperiment::assay(scee, "counts")
    samples <- as.character(cd[[sample_col]])
    cell_types <- as.character(cd[[celltype_col]])
    genes <- if (is.null(target_genes)) rownames(scee)
             else intersect(target_genes, rownames(scee))
    if (!length(genes)) stop("No target genes present in scee")

    ## shared cell-level design (donor, cell type, exposure, covariates) ----
    ct_f <- factor(cell_types)
    lvls <- levels(ct_f)
    meta <- data.frame(donor = factor(samples), celltype = ct_f,
                       exposure = as.numeric(exp_data[samples, exposure]),
                       stringsAsFactors = FALSE)
    if (length(covariates))
        for (cv in covariates) meta[[cv]] <- as.numeric(exp_data[samples, cv])
    cov_rhs <- if (length(covariates))
        paste("+", paste(covariates, collapse = " + ")) else ""
    multi_celltype <- length(lvls) > 1L
    re_rhs <- if (multi_celltype) {
        "+ (1 | donor) + (1 | donor:celltype)"
    } else {
        "+ (1 | donor)"
    }
    singular_ok <- lme4::.makeCC(action = "ignore", tol = 1e-4)

    if (family == "gaussian") {
        Y <- .log_cp10k_cell(counts_mat, genes)            # genes x cells
        fixed_rhs <- if (multi_celltype) "exposure * celltype" else "exposure"
        form <- stats::as.formula(paste("y ~", fixed_rhs, cov_rhs, re_rhs))
        control <- lme4::lmerControl(check.conv.singular = singular_ok)

        one_gene <- function(gi) {
            df <- meta; df$y <- as.numeric(Y[gi, ])
            fit <- tryCatch(lmerTest::lmer(form, data = df, REML = TRUE,
                                           control = control),
                            error = function(e) NULL)
            if (is.null(fit)) return(NULL)
            nm <- names(lme4::fixef(fit))
            if (!"exposure" %in% nm) return(NULL)
            rows <- lapply(lvls, function(ct) {
                L <- stats::setNames(numeric(length(nm)), nm)
                L["exposure"] <- 1
                if (ct != lvls[1]) {
                    ix <- paste0("exposure:celltype", ct)
                    # interaction aliased/dropped
                    if (!ix %in% nm) return(NULL)
                    L[ix] <- 1
                }
                tst <- tryCatch(
                    lmerTest::contest1D(fit, L, ddf = ddf, confint = FALSE),
                    error = function(e) NULL)
                if (is.null(tst)) return(NULL)
                ## y is already log2 CP10k, so the estimate IS a log2 effect.
                data.frame(gene = genes[gi], celltype = ct,
                           log2FC = tst[["Estimate"]], se = tst[["Std. Error"]],
                           statistic = tst[["t value"]], df = tst[["df"]],
                           pvalue = tst[["Pr(>|t|)"]],
                           stringsAsFactors = FALSE)
            })
            do.call(rbind, Filter(Negate(is.null), rows))
        }
        res <- .bpapply(seq_along(genes), one_gene, BPPARAM = BPPARAM)
        method_lbl <- paste0("lmerTest/", ddf)

    } else {
        ## nbinom: cell-type-nested reparametrisation so each cell type's
        ## exposure effect is a single coefficient; test by LRT (non-Wald).
        message("family='nbinom': Satterthwaite/KR are LMM-only; using ",
                "likelihood-ratio tests (ddf ignored). This refits glmer.nb ",
                "per cell type per gene -- prefer run_sc_exwas_pb_offset for ",
                "count data.")
        Yc <- as.matrix(counts_mat[genes, , drop = FALSE])
        meta$logoff <- log(pmax(Matrix::colSums(counts_mat), 1))  # cell offset
        expo_cols <- paste0("expo_", make.names(lvls))
        for (k in seq_along(lvls))
            meta[[expo_cols[k]]] <- meta$exposure * (meta$celltype == lvls[k])
        full_rhs <- paste("0 + celltype +",
                          paste(expo_cols, collapse = " + "), cov_rhs,
                          "+ offset(logoff)", re_rhs)
        full_form <- stats::as.formula(paste("y ~", full_rhs))
        control <- lme4::glmerControl(check.conv.singular = singular_ok)

        one_gene <- function(gi) {
            df <- meta
            df$y <- as.integer(round(pmax(Yc[gi, ], 0)))
            full <- tryCatch(lme4::glmer.nb(full_form, data = df,
                                            control = control),
                             error = function(e) NULL)
            if (is.null(full)) return(NULL)
            fe <- lme4::fixef(full)
            V <- tryCatch(as.matrix(stats::vcov(full)),
                          error = function(e) NULL)
            rows <- lapply(seq_along(lvls), function(k) {
                cn <- expo_cols[k]
                if (!cn %in% names(fe)) return(NULL)
                beta <- fe[[cn]]
                se_nat <- if (!is.null(V) && cn %in% rownames(V))
                    sqrt(max(V[cn, cn], 0)) else NA_real_
                other_slopes <- setdiff(expo_cols, cn)
                red_rhs <- paste("0 + celltype",
                                 if (length(other_slopes))
                                     paste("+", paste(other_slopes,
                                                      collapse = " + ")),
                                 cov_rhs, "+ offset(logoff)", re_rhs)
                red <- tryCatch(lme4::glmer.nb(
                    stats::as.formula(paste("y ~", red_rhs)), data = df,
                    control = control),
                    error = function(e) NULL)
                if (is.null(red)) {
                    chi <- NA_real_
                    pv <- NA_real_
                } else {
                    an <- tryCatch(stats::anova(red, full),
                                   error = function(e) NULL)
                    if (is.null(an)) {
                        chi <- NA_real_
                        pv <- NA_real_
                    } else {
                        chi <- an[["Chisq"]][2]
                        pv <- an[["Pr(>Chisq)"]][2]
                    }
                }
                ## glmer.nb uses a natural-log link -> rescale to log2.
                data.frame(gene = genes[gi], celltype = lvls[k],
                           log2FC = beta / log(2), se = se_nat / log(2),
                           statistic = chi, df = 1, pvalue = pv,
                           stringsAsFactors = FALSE)
            })
            do.call(rbind, Filter(Negate(is.null), rows))
        }
        res <- .bpapply(seq_along(genes), one_gene, BPPARAM = BPPARAM)
        method_lbl <- "lme4::glmer.nb/LRT"
    }

    out <- do.call(rbind, Filter(Negate(is.null), res))
    if (is.null(out) || !nrow(out))
        stop("No gene/cell-type stratum yielded an estimable model")
    out$padj <- stats::p.adjust(out$pvalue, "BH")
    out$method <- method_lbl
    out <- out[, c("gene", "celltype", "log2FC", "se", "statistic",
                   "df", "pvalue", "padj", "method")]
    out[order(out$padj, out$pvalue), , drop = FALSE]
}
