# R/cell-coupling.R
# Cell-level gene-protein coupling and its association with an exposure

# Suppress R CMD check NOTEs for metafor formula variables
utils::globalVariables(c("yi", "vi"))

#' Cell-Level Gene-Protein Coupling Analysis
#'
#' Computes within-donor Spearman partial correlations between
#' gene expression and protein abundance at the single-cell level,
#' then uses random-effects meta-regression to test whether
#' coupling strength varies with a donor-level exposure.
#'
#' @param x A \code{SingleCellExposomeExperiment} object.
#'   Each cell's donor is read from the \code{colData} column named by
#'   \code{\link{sampleColumn}(x)}. Must contain an altExp named
#'   \code{altexp_name} with cell-level protein counts.
#' @param genes Character vector of gene names to test.
#' @param proteins Character vector of protein names to test
#'   (must be present in \code{altExp}).
#' @param exposure Character; name of the exposure variable
#'   in \code{exposureData(x)}.
#' @param celltype Character; cell type to analyse.
#' @param celltype_col Character; column in \code{colData}
#'   containing cell type labels. Default \code{"cell_type"}.
#' @param altexp_name Character; name of the altExp containing
#'   protein data. Default \code{"CITE"}.
#' @param confounders Character vector of column names in
#'   \code{colData(x)} to partial out (e.g.,
#'   \code{c("S.Score", "G2M.Score", "percent.mt")}).
#'   Library sizes (nCount_RNA, nCount_protein) are always
#'   controlled regardless of this parameter. Each confounder is ranked
#'   within a donor, so it must be numeric, logical, an ordered factor or
#'   a label with at most two values; an unordered label with more
#'   categories, such as a cell-cycle phase, is an error and must be
#'   given as 0/1 indicator columns for all its categories but one,
#'   leaving out a category that every donor has cells of (for G1, S and
#'   G2M, for example S and G2M when every donor has G1 cells). The names
#'   \code{gene}, \code{protein}, \code{lib_rna} and \code{lib_prot} are
#'   used for the variables of the partial correlation and cannot be
#'   confounder names. Nor can a confounder be a column that is not passed
#'   on as one column under its own name: a matrix or nested table with
#'   several columns or with a column name other than that of the
#'   confounder, or a column with the name of an argument of
#'   \code{\link[base]{data.frame}} (\code{row.names}, \code{check.rows},
#'   \code{check.names}, \code{fix.empty.names} or
#'   \code{stringsAsFactors}).
#'   Default \code{NULL} (library sizes only).
#' @param min_cells Integer; minimum cells per donor to include.
#'   Default 30.
#' @param min_donors Integer; minimum number of donors with enough cells
#'   (at least \code{min_cells}, and at least 7 plus the number of
#'   confounders, because the two library sizes are always adjusted for)
#'   and values of the exposure and of every moderator. A gene-protein pair
#'   also needs a partial correlation in that many donors. Default 10.
#' @param meta_moderators Character vector or one-sided formula of
#'   additional donor-level moderators for the meta-regression
#'   (e.g., \code{"age"}, \code{c("age", "sex")} or
#'   \code{~ age + sex}). These must be columns in
#'   \code{exposureData(x)} with syntactic R names (see
#'   \code{\link[base]{make.names}}) other than \code{.}, \code{...} and
#'   names such as \code{..1}, which R reads specially in a formula; a
#'   formula may only list such columns, without transformations or
#'   interactions. The exposure itself
#'   cannot be a moderator, and the names \code{yi}, \code{vi},
#'   \code{donor_exp}, \code{intrcpt} and \code{meta_formula}, which the
#'   meta-regression uses for its own variables, cannot be moderator names.
#'   Default \code{NULL} (exposure only).
#' @param var_method Character. Method for estimating the
#'   sampling variance of Fisher z-transformed correlations.
#'   \code{"asymptotic"} (default) uses the formula
#'   \eqn{1.06 / (n - q - 3)}. \code{"bootstrap"} resamples
#'   cells within each donor to compute empirical variance,
#'   which is more accurate for small cell counts or
#'   non-normal data but slower. The bootstrap uses the random
#'   number generator, so set a seed for reproducible variances.
#' @param n_boot Integer. Number of bootstrap replicates
#'   per donor when \code{var_method = "bootstrap"}.
#'   Default 200.
#'
#' @return A \code{DataFrame} with one row per gene-protein pair
#'   (zero rows when too few donors qualify). A pair whose exposure does
#'   not vary among its donors, or is collinear with the moderators, is
#'   not tested, with a warning that counts such pairs:
#'   \describe{
#'     \item{gene, protein}{Feature names.}
#'     \item{n_donors}{Number of donors in the meta-regression: donors
#'       with enough cells, a partial correlation and values of the
#'       exposure and of every moderator.}
#'     \item{beta0, se0, pval0}{Intercept (baseline coupling).}
#'     \item{beta1, se1, pval1}{Exposure coefficient: the slope, across
#'       donors, of the Fisher z of the coupling on the exposure.}
#'     \item{tau2}{Between-donor heterogeneity variance (REML).}
#'     \item{I2}{Heterogeneity percentage.}
#'     \item{mean_r}{Mean partial correlation across donors.}
#'     \item{mean_n_cells}{Mean cells per donor.}
#'     \item{n_confounders}{Number of conditioning variables
#'       requested (library sizes plus confounders); variables
#'       constant within a donor are dropped for that donor.}
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
#' (2 + length of \code{confounders}, less any that are constant
#' within the donor and therefore dropped). The factor 1.06
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
#'     colData = S4Vectors::DataFrame(donor_id = donor, cell_type = "T"))
#' SingleCellExperiment::altExp(sce, "CITE") <-
#'     SummarizedExperiment::SummarizedExperiment(
#'         assays = list(counts = protein))
#' exp_mat <- matrix(seq(0, 2, length.out = 12), ncol = 1,
#'     dimnames = list(donors, "exposure"))
#' scee <- build_scee(sce, exp_mat, sample_col = "donor_id")
#' if (requireNamespace("metafor", quietly = TRUE) &&
#'     requireNamespace("ppcor", quietly = TRUE)) {
#'     run_cell_coupling(scee, genes = "Gene1", proteins = "Prot1",
#'         exposure = "exposure", celltype = "T",
#'         min_cells = 20L, min_donors = 10L)
#' }
run_cell_coupling <- function(x,
                              genes,
                              proteins,
                              exposure,
                              celltype,
                              celltype_col = "cell_type",
                              altexp_name = "CITE",
                              confounders = NULL,
                              min_cells = 30L,
                              min_donors = 10L,
                              meta_moderators = NULL,
                              var_method = c("asymptotic",
                                  "bootstrap"),
                              n_boot = 200L) {

    stopifnot(is(x, "SingleCellExposomeExperiment"))
    exposure <- .single_name(exposure, "exposure")
    celltype <- .single_name(celltype, "celltype")
    stopifnot(altexp_name %in% altExpNames(x))
    sample_col <- sampleColumn(x)
    var_method <- match.arg(var_method)

    if (!requireNamespace("metafor", quietly = TRUE))
        stop("Package 'metafor' required. ",
             "Install with: install.packages('metafor')")
    if (!requireNamespace("ppcor", quietly = TRUE))
        stop("Package 'ppcor' required. ",
             "Install with: install.packages('ppcor')")

    ## Validate exposure
    exp_data <- exposureData(x)
    if (!exposure %in% colnames(exp_data))
        stop("Exposure '", exposure, "' not in exposureData")
    exp_vec <- setNames(exp_data[, exposure], rownames(exp_data))

    ## Validate meta_moderators: column names, or a formula that lists them
    if (inherits(meta_moderators, "formula"))
        meta_moderators <- .moderator_names(meta_moderators)
    mod_mat <- NULL
    if (!is.null(meta_moderators)) {
        meta_moderators <- as.character(meta_moderators)
        missing_mods <- setdiff(meta_moderators, colnames(exp_data))
        if (length(missing_mods) > 0)
            stop("Moderators not in exposureData: ",
                 paste(missing_mods, collapse = ", "))
        .check_formula_names(meta_moderators, "moderator")
        if (exposure %in% meta_moderators)
            stop("The exposure must not also be listed in ",
                "'meta_moderators'.", call. = FALSE)
        ## the meta-regression data hold yi, vi and donor_exp, metafor calls
        ## the intercept intrcpt, and it looks meta_formula, the name its
        ## 'mods' argument is given below, up in those data first
        .check_reserved_names(meta_moderators, "moderator",
            c("yi", "vi", "donor_exp", "intrcpt", "meta_formula"))
        mod_mat <- exp_data[, meta_moderators, drop = FALSE]
    }

    ## Validate altExp
    prot_se <- altExp(x, altexp_name)
    missing_prot <- setdiff(proteins, rownames(prot_se))
    if (length(missing_prot) > 0)
        warning("Proteins not found in altExp: ",
                paste(head(missing_prot, 5), collapse = ", "))
    proteins <- intersect(proteins, rownames(prot_se))
    if (length(proteins) == 0)
        stop("No matching proteins in altExp")

    ## Validate genes
    missing_genes <- setdiff(genes, rownames(x))
    if (length(missing_genes) > 0)
        warning("Genes not found: ",
                paste(head(missing_genes, 5), collapse = ", "))
    genes <- intersect(genes, rownames(x))
    if (length(genes) == 0)
        stop("No matching genes in SCE")

    ## Validate confounders
    cd <- colData(x)
    if (!is.null(confounders)) {
        missing_conf <- setdiff(confounders, colnames(cd))
        if (length(missing_conf) > 0)
            stop("Confounders not in colData: ",
                 paste(missing_conf, collapse = ", "))
        ## the columns of the partial correlation in .coupling_donor_stat()
        .check_reserved_names(confounders, "confounder",
            c("gene", "protein", "lib_rna", "lib_prot"), "colData(x)")
    }

    ## Cell type indices
    ct_idx <- which(cd[[celltype_col]] == celltype)
    if (length(ct_idx) == 0)
        stop("No cells found for celltype '", celltype, "'")

    ct_donors <- as.character(cd[[sample_col]][ct_idx])

    ## Pre-compute library sizes (always controlled)
    ncount_rna  <- Matrix::colSums(counts(x)[, ct_idx, drop = FALSE])
    prot_counts <- assay(prot_se, "counts")
    ncount_prot <- Matrix::colSums(prot_counts[, ct_idx, drop = FALSE])

    ## Pre-extract confounder data for cell type
    conf_data <- NULL
    n_conf <- 2L  # library sizes are always included
    if (!is.null(confounders)) {
        conf_data <- .coupling_confounders(cd, ct_idx, confounders)
        n_conf <- n_conf + length(confounders)
        ## each confounder is ranked within a donor, so an unordered label
        ## with more than two values would be adjusted for in an arbitrary
        ## order (G1 < G2M < S)
        categorical <- vapply(conf_data, function(value) {
            (is.character(value) || (is.factor(value) && !is.ordered(value))) &&
                length(unique(value[!is.na(value)])) > 2L
        }, logical(1))
        if (any(categorical)) {
            categorical_names <- paste(confounders[categorical],
                collapse = ", ")
            stop("Confounder(s) with more than two unordered categories: ",
                categorical_names, ". Each confounder is ranked, so give ",
                "such a variable as 0/1 indicator columns of colData(x), ",
                "one for each of its categories except one that every ",
                "donor has.", call. = FALSE)
        }
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
        return(S4Vectors::DataFrame())
    }

    ## A donor enters the meta-regression only with an exposure value and a
    ## value of every moderator
    analysable <- !is.na(exp_vec)
    if (!is.null(mod_mat))
        analysable <- analysable & stats::complete.cases(mod_mat)
    n_incomplete <- sum(!analysable[valid_donors])
    if (n_incomplete > 0L)
        message("[exposomeSC] Leaving ", n_incomplete, " donor(s) with a ",
                "missing exposure or moderator value out of the ",
                "meta-regression.")
    n_analysable <- length(valid_donors) - n_incomplete
    if (n_analysable < min_donors) {
        warning("Only ", n_analysable, " donors with >= ", min_for_pcor,
                " cells and complete exposure and moderator values (need ",
                min_donors, ")")
        return(S4Vectors::DataFrame())
    }

    n_pairs <- length(genes) * length(proteins)
    fmt <- paste0("[exposomeSC] Cell coupling: %s, %d genes x %d proteins ",
                  "= %d pairs, %d confounders")
    message(sprintf(fmt, celltype, length(genes), length(proteins),
                    n_pairs, n_conf))

    ## Build meta-regression formula
    meta_formula <- as.formula("~ donor_exp")
    if (!is.null(meta_moderators)) {
        mod_terms <- paste(meta_moderators, collapse = " + ")
        ## the exposure comes last: metafor drops the later of two collinear
        ## columns, so an exposure collinear with the moderators is dropped
        ## and the pair is counted as not tested below
        meta_formula <- as.formula(paste("~", mod_terms, "+ donor_exp"))
    }

    ## Results accumulator
    results <- vector("list", n_pairs)
    idx <- 0L
    n_untestable <- 0L

    for (gene in genes) {
        gene_row <- match(gene, rownames(x))
        rna_all <- counts(x)[gene_row, ct_idx]
        if (is(rna_all, "sparseVector"))
            rna_all <- as.numeric(rna_all)

        for (protein in proteins) {
            idx <- idx + 1L
            prot_row <- match(protein, rownames(prot_se))
            prot_all <- prot_counts[prot_row, ct_idx]
            if (is(prot_all, "sparseVector"))
                prot_all <- as.numeric(prot_all)

            ## Per-donor partial correlation, Fisher z and its variance
            n_valid <- length(valid_donors)
            donor_z <- rep(NA_real_, n_valid)
            donor_v <- rep(NA_real_, n_valid)
            donor_r <- rep(NA_real_, n_valid)
            donor_n <- rep(NA_integer_, n_valid)
            for (k in seq_len(n_valid)) {
                d <- valid_donors[k]
                d_idx <- which(ct_donors == d)
                stat <- .coupling_donor_stat(
                    rna_d = rna_all[d_idx],
                    prot_d = prot_all[d_idx],
                    lib_rna_d = ncount_rna[d_idx],
                    lib_prot_d = ncount_prot[d_idx],
                    conf_d = if (is.null(conf_data)) NULL else
                        conf_data[d_idx, , drop = FALSE],
                    min_for_pcor = min_for_pcor,
                    var_method = var_method,
                    n_boot = n_boot
                )
                ## the exposure and moderators are checked after the donor
                ## statistic so that bootstrap draws do not depend on
                ## missing values
                if (is.null(stat) || !analysable[[d]]) next
                donor_z[k] <- stat[["z"]]
                donor_v[k] <- stat[["v"]]
                donor_r[k] <- stat[["r"]]
                donor_n[k] <- as.integer(stat[["n"]])
            }
            used <- !is.na(donor_z)
            donor_ids <- valid_donors[used]
            donor_z <- donor_z[used]
            donor_v <- donor_v[used]
            donor_r <- donor_r[used]
            donor_n <- donor_n[used]
            donor_exp <- unname(exp_vec[donor_ids])

            if (length(donor_z) < min_donors) next
            ## the meta-regression needs an exposure that varies
            if (length(unique(donor_exp)) < 2L) {
                n_untestable <- n_untestable + 1L
                next
            }

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
            ## coefficients are taken by name: metafor drops a predictor
            ## that is collinear with the others
            k0 <- match("intrcpt", rownames(fit$beta))
            k1 <- match("donor_exp", rownames(fit$beta))
            if (is.na(k0) || is.na(k1)) {
                n_untestable <- n_untestable + 1L
                next
            }

            results[[idx]] <- data.frame(
                gene = gene,
                protein = protein,
                n_donors = length(donor_z),
                beta0 = fit$beta[k0, 1],
                se0 = fit$se[k0],
                pval0 = fit$pval[k0],
                beta1 = fit$beta[k1, 1],
                se1 = fit$se[k1],
                pval1 = fit$pval[k1],
                tau2 = fit$tau2,
                I2 = fit$I2,
                mean_r = mean(donor_r),
                mean_n_cells = mean(donor_n),
                n_confounders = n_conf,
                ## row.names = NULL: the coefficient taken from fit$beta
                ## would otherwise name the row "intrcpt"
                row.names = NULL,
                stringsAsFactors = FALSE)
        }
    }

    if (n_untestable > 0L)
        warning(n_untestable, " gene-protein pair(s) not tested: the ",
                "exposure does not vary among their donors, or is collinear ",
                "with the moderators.", call. = FALSE)
    out <- do.call(rbind, results[!vapply(results, is.null, logical(1))])
    if (is.null(out)) out <- data.frame()

    message(sprintf("[exposomeSC] %d pairs tested for %s",
                    nrow(out), celltype))
    S4Vectors::DataFrame(out, check.names = FALSE)
}

# The moderators named by a one-sided formula such as ~ age + sex. A term
# that transforms or combines variables, such as I(age^2) or age:sex, has no
# column to take, so such a formula is refused rather than reduced to its
# variables.
.moderator_names <- function(moderators) {
    model_terms <- tryCatch(stats::terms(moderators),
                            error = function(e) NULL)
    labels <- if (is.null(model_terms)) NULL else
        gsub("`", "", attr(model_terms, "term.labels"), fixed = TRUE)
    if (length(moderators) != 2L || !length(labels) ||
        attr(model_terms, "intercept") != 1L ||
        !identical(labels, all.vars(moderators))) {
        stop("A 'meta_moderators' formula must be one-sided and only list ",
            "columns of exposureData(x), for example ~ age + sex.",
            call. = FALSE)
    }
    labels
}

# The confounder columns of the chosen cells as a data frame. The conversion
# hands the columns of colData(x) to data.frame() by name, so a column
# called like one of its arguments (stringsAsFactors, say) is taken for that
# argument, and a matrix or a nested table arrives under the names of its
# own columns. A confounder that does not arrive as one column under its own
# name would be left out of the partial correlation or put in place of
# another variable, so it is refused.
.coupling_confounders <- function(cd, cells, confounders) {
    selected <- cd[cells, confounders, drop = FALSE]
    wanted <- colnames(selected)
    converted <- function(k) {
        tryCatch(as.data.frame(selected[, k, drop = FALSE]),
            error = function(e) NULL)
    }
    arrives <- function(columns, names) {
        identical(names(columns), names) &&
            all(vapply(columns, NCOL, integer(1)) == 1L)
    }
    conf_data <- converted(seq_along(wanted))
    if (!arrives(conf_data, wanted)) {
        usable <- vapply(seq_along(wanted), function(k) {
            arrives(converted(k), wanted[k])
        }, logical(1))
        unusable <- unique(wanted[!usable])
        if (!length(unusable)) unusable <- unique(wanted)
        stop("Confounder(s) that cannot be used: ",
            paste0("'", unusable, "'", collapse = ", "), ". A confounder ",
            "cannot be a matrix or nested table with several columns or ",
            "with a column name other than that of the confounder, and ",
            "cannot have the name of an argument of data.frame() ",
            "(row.names, check.rows, check.names, fix.empty.names, ",
            "stringsAsFactors).", call. = FALSE)
    }
    conf_data
}

# Spearman partial correlation of gene and protein given the other columns
# of dat (all already ranked); NA when it cannot be computed.
.coupling_partial_r <- function(dat, q_d) {
    if (q_d == 0L) {
        return(stats::cor(dat$gene, dat$protein, method = "spearman"))
    }
    pc <- tryCatch(
        ppcor::pcor(dat, method = "spearman"),
        error = function(e) NULL)
    if (is.null(pc)) NA_real_ else pc$estimate[1, 2]
}

# Bootstrap variance of the Fisher z of one donor, resampling its cells; NA
# when fewer than half of the replicates give a finite value.
.coupling_boot_var <- function(dat, q_d, n_boot) {
    boot_z <- vapply(seq_len(n_boot), function(b) {
        dat_b <- dat[sample(nrow(dat), replace = TRUE), , drop = FALSE]
        if (sd(dat_b$gene) < 1e-10 || sd(dat_b$protein) < 1e-10)
            return(NA_real_)
        r_b <- .coupling_partial_r(dat_b, q_d)
        if (!is.finite(r_b)) return(NA_real_)
        atanh(max(min(r_b, 0.999), -0.999))
    }, numeric(1))
    boot_z <- boot_z[is.finite(boot_z)]
    if (length(boot_z) < n_boot * 0.5) return(NA_real_)
    var(boot_z)
}

# Stage 1 for one donor: Spearman partial correlation of gene and protein
# given library sizes and confounders, its Fisher z and the variance of z.
# Returns NULL when the donor cannot contribute.
.coupling_donor_stat <- function(rna_d, prot_d, lib_rna_d, lib_prot_d,
                                 conf_d, min_for_pcor, var_method, n_boot) {
    n_d <- length(rna_d)
    if (n_d < min_for_pcor) return(NULL)
    if (sd(rna_d) < 1e-10 || sd(prot_d) < 1e-10) return(NULL)

    ## Columns: gene, protein, lib_rna, lib_prot, [confounders];
    ## run_cell_coupling() refuses confounders with the first four names
    dat <- data.frame(
        gene = rank(rna_d),
        protein = rank(prot_d),
        lib_rna = rank(lib_rna_d),
        lib_prot = rank(lib_prot_d))
    if (!is.null(conf_d)) {
        for (cv in colnames(conf_d)) {
            dat[[cv]] <- rank(conf_d[[cv]], na.last = "keep")
        }
        ## Remove rows with NA confounders
        complete <- complete.cases(dat)
        if (sum(complete) < min_for_pcor) return(NULL)
        dat <- dat[complete, , drop = FALSE]
        n_d <- nrow(dat)
    }

    ## A conditioning variable that is constant within this donor carries
    ## no information and would make the partial correlation singular, so it
    ## is dropped for this donor
    conditioning <- setdiff(names(dat), c("gene", "protein"))
    constant <- conditioning[vapply(dat[conditioning],
        function(v) isTRUE(sd(v) < 1e-10), logical(1))]
    if (length(constant))
        dat <- dat[, setdiff(names(dat), constant), drop = FALSE]
    q_d <- ncol(dat) - 2L

    r_p <- .coupling_partial_r(dat, q_d)
    if (!is.finite(r_p)) return(NULL)
    z_val <- atanh(max(min(r_p, 0.999), -0.999))

    v_val <- if (var_method == "bootstrap") {
        .coupling_boot_var(dat, q_d, n_boot)
    } else {
        ## Asymptotic: 1.06 / (n - q - 3)
        1.06 / (n_d - q_d - 3)
    }
    if (!is.finite(v_val) || v_val <= 0) return(NULL)
    c(z = z_val, v = v_val, r = r_p, n = n_d)
}
