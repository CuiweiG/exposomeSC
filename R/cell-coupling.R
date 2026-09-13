# R/cell-coupling.R
# Cell-level gene-protein coupling and exposure-driven rewiring

# Suppress R CMD check NOTEs for metafor formula variables
utils::globalVariables(c("yi", "vi"))

#' Cell-Level Gene-Protein Coupling Analysis
#'
#' Computes within-donor Spearman partial correlations between
#' gene expression and protein abundance at the single-cell level,
#' then uses random-effects meta-regression to test whether
#' coupling strength varies with a donor-level exposure.
#'
#' @param scee A \code{SingleCellExposomeExperiment} object.
#'   Must contain an altExp named \code{altexp_name} with
#'   cell-level protein counts.
#' @param genes Character vector of gene names to test.
#' @param proteins Character vector of protein names to test
#'   (must be present in \code{altExp}).
#' @param exposure Character; name of the exposure variable
#'   in \code{exposureData(scee)}.
#' @param celltype Character; cell type to analyse.
#' @param celltype_col Character; column in \code{colData}
#'   containing cell type labels. Default \code{"celltype"}.
#' @param sample_col Character; column in \code{colData}
#'   containing donor/sample IDs. Default \code{"donor"}.
#' @param altexp_name Character; name of the altExp containing
#'   protein data. Default \code{"CITE"}.
#' @param confounders Character vector of column names in
#'   \code{colData(scee)} to partial out (e.g.,
#'   \code{c("S.Score", "G2M.Score", "percent.mt")}).
#'   Library sizes (nCount_RNA, nCount_protein) are always
#'   controlled regardless of this parameter.
#'   Default \code{NULL} (library sizes only).
#' @param min_cells Integer; minimum cells per donor to include.
#'   Default 30.
#' @param min_donors Integer; minimum donors with sufficient
#'   cells. Default 10.
#' @param meta_moderators Formula or character vector of
#'   additional donor-level moderators for the meta-regression
#'   (e.g., \code{"age"} or \code{c("age", "sex")}).
#'   These must be columns in \code{exposureData(scee)}.
#'   Default \code{NULL} (exposure only).
#' @param var_method Character. Method for estimating the
#'   sampling variance of Fisher z-transformed correlations.
#'   \code{"asymptotic"} (default) uses the formula
#'   \eqn{1.06 / (n - q - 3)}. \code{"bootstrap"} resamples
#'   cells within each donor to compute empirical variance,
#'   which is more accurate for small cell counts or
#'   non-normal data but slower.
#' @param n_boot Integer. Number of bootstrap replicates
#'   per donor when \code{var_method = "bootstrap"}.
#'   Default 200.
#'
#' @return A \code{data.frame} with one row per gene-protein pair:
#'   \describe{
#'     \item{gene, protein}{Feature names.}
#'     \item{n_donors}{Number of donors included.}
#'     \item{beta0, se0, pval0}{Intercept (baseline coupling).}
#'     \item{beta1, se1, pval1}{Exposure coefficient (rewiring).}
#'     \item{tau2}{Between-donor heterogeneity variance (REML).}
#'     \item{I2}{Heterogeneity percentage.}
#'     \item{mean_r}{Mean partial correlation across donors.}
#'     \item{mean_n_cells}{Mean cells per donor.}
#'     \item{n_confounders}{Number of variables partialed out.}
#'   }
#'
#' @details
#' For each gene-protein pair, the method proceeds in two stages:
#'
#' \strong{Stage 1 (within-donor):} For each donor \eqn{d} with
#' at least \code{min_cells} cells of the specified cell type,
#' compute the Spearman partial correlation between gene
#' expression and protein abundance, controlling for library
#' sizes and any user-specified confounders:
#' \deqn{r_d = \rho_S(\text{gene}_i, \text{protein}_j \mid
#'   \text{nCount\_RNA}, \text{nCount\_protein},
#'   \text{confounders})}
#' Then apply Fisher's z-transform:
#' \deqn{z_d = \text{arctanh}(r_d), \quad
#'   \text{Var}(z_d) \approx 1.06 / (n_d - q - 3)}
#' where \eqn{q} is the number of conditioning variables
#' (2 + length of \code{confounders}). The factor 1.06
#' corrects for the efficiency loss of rank-based correlation
#' (Fieller et al. 1957).
#'
#' \strong{Stage 2 (across-donor meta-regression):} Fit a
#' random-effects model via REML with Knapp-Hartung adjustment:
#' \deqn{z_d = \beta_0 + \beta_1 \cdot \text{exposure}_d
#'   + \sum_k \gamma_k \cdot \text{moderator}_{dk}
#'   + u_d + \varepsilon_d}
#' where \eqn{\varepsilon_d \sim N(0, v_d)} (known sampling
#' variance) and \eqn{u_d \sim N(0, \tau^2)} (between-donor
#' heterogeneity). The Knapp-Hartung adjustment provides
#' more accurate confidence intervals for small numbers of
#' studies (donors), using a t-distribution with
#' \eqn{D - p - 1} degrees of freedom rather than a
#' z-distribution.
#'
#' \strong{Confounder control:} Cell-level confounders (cell
#' cycle, mitochondrial content, ambient RNA scores) are
#' partialed out in Stage 1, ensuring that the correlation
#' reflects genuine gene-protein coupling rather than
#' shared technical or biological confounds. Donor-level
#' moderators (age, sex) are included in Stage 2.
#'
#' @references
#' Viechtbauer W (2010). Conducting meta-analyses in R with
#' the metafor package. \emph{J Stat Softw}, 36(3), 1-48.
#' \doi{10.18637/jss.v036.i03}
#'
#' Fieller EC, Hartley HO, Pearson ES (1957). Tests for rank
#' correlation coefficients. I. \emph{Biometrika}, 44, 470-481.
#'
#' Kim S (2015). ppcor: An R Package for a Fast Calculation
#' to Semi-partial Correlation Coefficients. \emph{Commun Stat
#' Appl Methods}, 22(6), 665-674.
#'
#' @importFrom stats sd complete.cases
#' @importFrom SingleCellExperiment altExpNames altExp counts
#' @importFrom SummarizedExperiment assay colData
#' @importFrom Matrix colSums
#' @export
#' @examples
#' set.seed(1)
#' donors <- sprintf("D%02d", 1:12)
#' donor <- rep(donors, each = 60)
#' gene <- matrix(stats::rpois(2 * length(donor), 10), nrow = 2,
#'     dimnames = list(c("Gene1", "Gene2"), paste0("c", seq_along(donor))))
#' protein <- matrix(stats::rpois(2 * length(donor), 5), nrow = 2,
#'     dimnames = list(c("Prot1", "Prot2"), colnames(gene)))
#' sce <- SingleCellExperiment::SingleCellExperiment(
#'     assays = list(counts = gene),
#'     colData = S4Vectors::DataFrame(donor = donor, celltype = "T",
#'         pseudotime = stats::runif(length(donor))))
#' SingleCellExperiment::altExp(sce, "CITE") <-
#'     SummarizedExperiment::SummarizedExperiment(
#'         assays = list(counts = protein))
#' exp_mat <- matrix(seq(0, 2, length.out = 12), ncol = 1,
#'     dimnames = list(donors, "exposure"))
#' scee <- build_scee(sce, exp_mat, sample_col = "donor")
#' if (requireNamespace("metafor", quietly = TRUE) &&
#'     requireNamespace("ppcor", quietly = TRUE)) {
#'     run_cell_coupling(scee, genes = "Gene1", proteins = "Prot1",
#'         exposure = "exposure", celltype = "T",
#'         min_cells = 20L, min_donors = 10L)
#' }
run_cell_coupling <- function(scee,
                              genes,
                              proteins,
                              exposure,
                              celltype,
                              celltype_col = "celltype",
                              sample_col = "donor",
                              altexp_name = "CITE",
                              confounders = NULL,
                              min_cells = 30L,
                              min_donors = 10L,
                              meta_moderators = NULL,
                              var_method = c("asymptotic",
                                  "bootstrap"),
                              n_boot = 200L) {

    stopifnot(is(scee, "SingleCellExposomeExperiment"))
    stopifnot(altexp_name %in% altExpNames(scee))
    var_method <- match.arg(var_method)

    if (!requireNamespace("metafor", quietly = TRUE))
        stop("Package 'metafor' required. ",
             "Install with: install.packages('metafor')")
    if (!requireNamespace("ppcor", quietly = TRUE))
        stop("Package 'ppcor' required. ",
             "Install with: install.packages('ppcor')")

    ## Validate exposure
    exp_data <- exposureData(scee)
    if (!exposure %in% colnames(exp_data))
        stop("Exposure '", exposure, "' not in exposureData")
    exp_vec <- setNames(exp_data[, exposure], rownames(exp_data))

    ## Validate meta_moderators
    mod_mat <- NULL
    if (!is.null(meta_moderators)) {
        missing_mods <- setdiff(meta_moderators, colnames(exp_data))
        if (length(missing_mods) > 0)
            stop("Moderators not in exposureData: ",
                 paste(missing_mods, collapse = ", "))
        mod_mat <- exp_data[, meta_moderators, drop = FALSE]
    }

    ## Validate altExp
    prot_se <- altExp(scee, altexp_name)
    missing_prot <- setdiff(proteins, rownames(prot_se))
    if (length(missing_prot) > 0)
        warning("Proteins not found in altExp: ",
                paste(head(missing_prot, 5), collapse = ", "))
    proteins <- intersect(proteins, rownames(prot_se))
    if (length(proteins) == 0)
        stop("No matching proteins in altExp")

    ## Validate genes
    missing_genes <- setdiff(genes, rownames(scee))
    if (length(missing_genes) > 0)
        warning("Genes not found: ",
                paste(head(missing_genes, 5), collapse = ", "))
    genes <- intersect(genes, rownames(scee))
    if (length(genes) == 0)
        stop("No matching genes in SCE")

    ## Validate confounders
    cd <- colData(scee)
    if (!is.null(confounders)) {
        missing_conf <- setdiff(confounders, colnames(cd))
        if (length(missing_conf) > 0)
            stop("Confounders not in colData: ",
                 paste(missing_conf, collapse = ", "))
    }

    ## Cell type indices
    ct_idx <- which(cd[[celltype_col]] == celltype)
    if (length(ct_idx) == 0)
        stop("No cells found for celltype '", celltype, "'")

    ct_donors <- as.character(cd[[sample_col]][ct_idx])

    ## Pre-compute library sizes (always controlled)
    ncount_rna  <- Matrix::colSums(counts(scee)[, ct_idx, drop = FALSE])
    prot_counts <- assay(prot_se, "counts")
    ncount_prot <- Matrix::colSums(prot_counts[, ct_idx, drop = FALSE])

    ## Pre-extract confounder data for cell type
    conf_data <- NULL
    n_conf <- 2L  # library sizes are always included
    if (!is.null(confounders)) {
        conf_data <- as.data.frame(cd[ct_idx, confounders, drop = FALSE])
        n_conf <- n_conf + length(confounders)
    }

    ## Identify valid donors
    donor_tab <- table(ct_donors)
    ## Need enough cells for partial correlation:
    ## n_cells > n_conf + 3 + safety margin
    min_for_pcor <- max(min_cells, n_conf + 5L)
    valid_donors <- names(donor_tab[donor_tab >= min_for_pcor])
    if (length(valid_donors) < min_donors) {
        warning("Only ", length(valid_donors), " donors with >= ",
                min_for_pcor, " cells (need ", min_donors, ")")
        return(data.frame())
    }

    n_pairs <- length(genes) * length(proteins)
    message(sprintf(
        "[exposomeSC] Cell coupling: %s, %d genes x %d proteins = %d pairs, %d confounders",
        celltype, length(genes), length(proteins), n_pairs, n_conf))

    ## Build meta-regression formula
    meta_formula <- as.formula("~ donor_exp")
    if (!is.null(meta_moderators)) {
        mod_terms <- paste(meta_moderators, collapse = " + ")
        meta_formula <- as.formula(paste("~ donor_exp +", mod_terms))
    }

    ## Results accumulator
    results <- vector("list", n_pairs)
    idx <- 0L

    for (gene in genes) {
        gene_row <- match(gene, rownames(scee))
        rna_all <- counts(scee)[gene_row, ct_idx]
        if (is(rna_all, "sparseVector"))
            rna_all <- as.numeric(rna_all)

        for (protein in proteins) {
            idx <- idx + 1L
            prot_row <- match(protein, rownames(prot_se))
            prot_all <- prot_counts[prot_row, ct_idx]
            if (is(prot_all, "sparseVector"))
                prot_all <- as.numeric(prot_all)

            ## Per-donor partial correlation
            donor_z <- numeric()
            donor_v <- numeric()
            donor_r <- numeric()
            donor_n <- integer()
            donor_exp <- numeric()
            donor_ids <- character()

            for (d in valid_donors) {
                d_idx <- which(ct_donors == d)
                n_d <- length(d_idx)
                if (n_d < min_for_pcor) next

                rna_d  <- rna_all[d_idx]
                prot_d <- prot_all[d_idx]

                if (sd(rna_d) < 1e-10 || sd(prot_d) < 1e-10) next

                ## Build data frame for ppcor
                ## Columns: gene, protein, lib_rna, lib_prot, [confounders]
                dat <- data.frame(
                    gene = rank(rna_d),
                    protein = rank(prot_d),
                    lib_rna = rank(ncount_rna[d_idx]),
                    lib_prot = rank(ncount_prot[d_idx]))

                ## Add confounders if specified
                if (!is.null(conf_data)) {
                    for (cv in confounders) {
                        cv_vals <- conf_data[[cv]][d_idx]
                        if (sd(cv_vals, na.rm = TRUE) < 1e-10) {
                            ## No variance in confounder for this donor;
                            ## use a small jitter to avoid singular matrix
                            cv_vals <- cv_vals + rnorm(length(cv_vals), 0, 1e-8)
                        }
                        dat[[cv]] <- rank(cv_vals, na.last = "keep")
                    }
                    ## Remove rows with NA confounders
                    complete <- complete.cases(dat)
                    if (sum(complete) < min_for_pcor) next
                    dat <- dat[complete, , drop = FALSE]
                    n_d <- nrow(dat)
                }

                pc <- tryCatch(
                    ppcor::pcor(dat, method = "spearman"),
                    error = function(e) NULL)
                if (is.null(pc)) next

                ## r_partial is [1,2] (gene-protein)
                r_p <- pc$estimate[1, 2]
                r_clamped <- max(min(r_p, 0.999), -0.999)
                z_val <- atanh(r_clamped)

                ## Variance estimation
                if (var_method == "bootstrap") {
                    ## Bootstrap within-donor variance
                    boot_z <- vapply(seq_len(n_boot), function(b) {
                        bi <- sample(nrow(dat), replace = TRUE)
                        dat_b <- dat[bi, , drop = FALSE]
                        ## Check variance
                        if (sd(dat_b$gene) < 1e-10 ||
                            sd(dat_b$protein) < 1e-10)
                            return(NA_real_)
                        pc_b <- tryCatch(
                            ppcor::pcor(dat_b, method = "spearman"),
                            error = function(e) NULL)
                        if (is.null(pc_b)) return(NA_real_)
                        r_b <- pc_b$estimate[1, 2]
                        r_b <- max(min(r_b, 0.999), -0.999)
                        atanh(r_b)
                    }, numeric(1))
                    boot_z <- boot_z[is.finite(boot_z)]
                    if (length(boot_z) < n_boot * 0.5) next
                    v_val <- var(boot_z)
                } else {
                    ## Asymptotic: 1.06 / (n - q - 3)
                    v_val <- 1.06 / (n_d - n_conf - 3)
                }
                if (v_val <= 0 || !is.finite(v_val)) next

                e_d <- exp_vec[d]
                if (is.na(e_d)) next

                donor_z <- c(donor_z, z_val)
                donor_v <- c(donor_v, v_val)
                donor_r <- c(donor_r, r_p)
                donor_n <- c(donor_n, n_d)
                donor_exp <- c(donor_exp, e_d)
                donor_ids <- c(donor_ids, d)
            }

            if (length(donor_z) < min_donors) next

            ## Build meta-regression data
            meta_df <- data.frame(
                yi = donor_z, vi = donor_v,
                donor_exp = donor_exp)

            if (!is.null(mod_mat)) {
                for (mm in meta_moderators) {
                    meta_df[[mm]] <- mod_mat[donor_ids, mm]
                }
            }

            ## Meta-regression: REML + Knapp-Hartung
            fit <- tryCatch(
                metafor::rma(
                    yi = yi, vi = vi,
                    mods = meta_formula,
                    data = meta_df,
                    method = "REML", test = "knha"),
                error = function(e) NULL)
            if (is.null(fit)) next

            results[[idx]] <- data.frame(
                gene = gene,
                protein = protein,
                n_donors = length(donor_z),
                beta0 = fit$beta[1, 1],
                se0 = fit$se[1],
                pval0 = fit$pval[1],
                beta1 = fit$beta[2, 1],
                se1 = fit$se[2],
                pval1 = fit$pval[2],
                tau2 = fit$tau2,
                I2 = fit$I2,
                mean_r = mean(donor_r),
                mean_n_cells = mean(donor_n),
                n_confounders = n_conf,
                stringsAsFactors = FALSE)
        }
    }

    out <- do.call(rbind, results[!vapply(results, is.null, logical(1))])
    if (is.null(out)) out <- data.frame()

    message(sprintf("[exposomeSC] %d pairs tested for %s",
                    nrow(out), celltype))
    out
}
