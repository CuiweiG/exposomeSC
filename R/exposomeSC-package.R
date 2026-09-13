#' exposomeSC: Single-Cell Exposome-Wide Association Studies
#'
#' A framework for cell-type-specific exposome-wide association
#' studies (sc-ExWAS) bridging donor-level environmental
#' exposures with single-cell resolution gene expression.
#'
#' @section Core analysis:
#' \describe{
#'   \item{\code{\link{build_scee}}}{Build a
#'     SingleCellExposomeExperiment from SCE + exposure matrix}
#'   \item{\code{\link{run_sc_exwas}}}{Cell-type-specific ExWAS
#'     via pseudobulk (edgeR, DESeq2 or voom-dream backend)}
#'   \item{\code{\link{run_sc_exwas_glmm}}}{GLMM-based ExWAS
#'     with donor random effects}
#'   \item{\code{\link{run_multi_exwas}}}{Batch ExWAS across
#'     multiple exposures with global FDR}
#' }
#'
#' @section Experimental and legacy association utilities:
#' \describe{
#'   \item{\code{\link{run_causal_mediation}}}{Model-based exploratory
#'     mediation calculations requiring user-supplied identification assumptions}
#' }
#'
#' @section Dose-response and mixtures:
#' \describe{
#'   \item{\code{\link{run_dose_response}}}{Polynomial
#'     dose-response (AIC selection)}
#'   \item{\code{\link{run_dose_response_gam}}}{GAM
#'     dose-response via mgcv}
#'   \item{\code{\link{run_sc_mixture}}}{Quantile-scored
#'     mixture screening}
#'   \item{\code{\link{run_mixture_qgcomp}}}{Formal qgcomp
#'     with bootstrap CI}
#'   \item{\code{\link{run_interaction_test}}}{Exposure x
#'     celltype interaction F-test}
#' }
#'
#' @section Cell-level coupling:
#' \describe{
#'   \item{\code{\link{run_cell_coupling}}}{Gene-protein
#'     coupling via meta-regression with confounders}
#'   \item{\code{\link{run_state_coupling}}}{Cell-state
#'     coupling: exposure x trajectory interaction}
#' }
#'
#' @section Enrichment and meta-analysis:
#' \describe{
#'   \item{\code{\link{run_gsea}}}{Pathway enrichment via
#'     fgsea}
#'   \item{\code{\link{run_meta_exwas}}}{Cross-cohort
#'     meta-analysis (fixed/random effects)}
#'   \item{\code{\link{run_spatial_exwas}}}{Region-resolved
#'     ExWAS via SpatialExperiment}
#'   \item{\code{\link{estimate_power}}}{Power/sample size
#'     calculator}
#' }
#'
#' @section Preprocessing:
#' \describe{
#'   \item{\code{\link{exposure_impute_lod}}}{Below-LOD
#'     imputation}
#'   \item{\code{\link{as_scee}}}{Convert from rexposome}
#'   \item{\code{\link{seurat_to_exposure}}}{Extract from
#'     Seurat}
#' }
#'
#' @references
#' Squair JW et al. (2021). Confronting false discoveries in
#' single-cell differential expression. \emph{Nat Commun}
#' 12:5692. \doi{10.1038/s41467-021-25960-2}
#'
#' VanderWeele TJ (2015). \emph{Explanation in Causal
#' Inference}. Oxford University Press.
#'
#' VanderWeele TJ, Ding P (2017). Sensitivity analysis in
#' observational research: introducing the E-value.
#' \emph{Ann Intern Med} 167:268-274.
#' \doi{10.7326/M16-2607}
#'
#' Aitchison J (1986). \emph{The Statistical Analysis of
#' Compositional Data}. Chapman and Hall.
#'
#' Imai K, Keele L, Tingley D (2010). A general approach to
#' causal mediation analysis. \emph{Psychol Methods}
#' 15:309-334. \doi{10.1037/a0020761}
#'
#' @docType package
#' @name exposomeSC-package
#' @aliases exposomeSC
#' @keywords package
"_PACKAGE"
