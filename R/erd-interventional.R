# R/erd-interventional.R
# Rebuilt Exposure-Response Decomposition (ERD) on the interventional
# rate-ratio scale, estimated by counterfactual g-computation with a
# compositional (ILR) cell-type-composition mediator.
#
# This replaces the earlier difference-of-coefficients estimator, which does
# not identify a causal indirect effect (non-collapsibility of the log link +
# a multivariate simplex mediator whose natural effects are unidentified).
#
# Estimands (VanderWeele, Vansteelandt & Robins 2014; Vansteelandt & Daniel
# 2017), with G_{a'} a random draw from the population ILR-mediator law under
# exposure a' and psi(a,a') = E[Y_{a, G_{a'}}]:
#   IDE_RR = psi(a, a*) / psi(a*, a*)      within-composition DIRECT effect
#   IIE_RR = psi(a, a)  / psi(a, a*)       composition-DRIVEN effect
#   OE_RR  = IDE_RR * IIE_RR               exact multiplicative decomposition
# Closed form under an additive working model (no A x M interaction, exposure-
# independent Sigma): log IDE_RR = theta_A ; log IIE_RR = theta_M^T eta_A ;
# coordinate-wise log CIIE_k = theta_{M,k} eta_{A,k} (summing to log IIE_RR).
# The Monte-Carlo g-computation below is the default (interaction-agnostic);
# the closed form is returned as a cross-check.

#' @include AllClasses.R
#' @importFrom stats lm residuals cov coef model.matrix rnorm quantile
#' @importFrom stats p.adjust pnorm
NULL

# ---- ILR basis (canonical Egozcue orthonormal contrast) -----------------
# Returns a D x (D-1) matrix V with V^T V = I_{D-1} and V^T 1 = 0, so that
# M = t(V) %*% log(p) are ILR coordinates (Egozcue et al. 2003).
.ilr_V <- function(D) {
    if (D < 2L) stop("Need >= 2 cell types for an ILR basis")
    V <- matrix(0, nrow = D, ncol = D - 1L)
    for (i in seq_len(D - 1L)) {
        r <- D - i
        V[i, i] <- sqrt(r / (r + 1))
        V[(i + 1L):D, i] <- -sqrt(1 / (r * (r + 1)))
    }
    V
}

# Zero-replace a proportion matrix (rows sum to 1) by a small multiplicative
# replacement, then re-close. delta defaults below the smallest nonzero mass.
.zero_replace <- function(P, delta = NULL) {
    if (is.null(delta)) {
        pos <- P[P > 0]
        delta <- if (length(pos)) min(pos) / 2 else 1e-6
    }
    P[P <= 0] <- delta
    P / rowSums(P)
}

# donors x (D-1) ILR coordinate matrix from a donors x D proportion matrix.
.ilr_coords <- function(P, V) {
    P <- .zero_replace(P)
    log(P) %*% V            # (n x D) %*% (D x (D-1)) = n x (D-1)
}

#' Mediational E-value (Smith & VanderWeele 2019) on the rate-ratio scale
#'
#' @param rr Numeric rate ratio (point estimate or CI limit).
#' @return The E-value \code{rr + sqrt(rr*(rr-1))}, with \code{rr} inverted
#'   when below 1; \code{1} when the limit already crosses the null.
#' @examples
#' mediational_evalue(c(1, 1.5, 0.5))
#' @export
mediational_evalue <- function(rr) {
    rr <- ifelse(is.finite(rr) & rr > 0 & rr < 1, 1 / rr, rr)
    out <- ifelse(is.finite(rr) & rr >= 1, rr + sqrt(rr * (rr - 1)), NA_real_)
    out[is.finite(rr) & rr <= 1] <- 1
    out
}

#' E-value for the confidence limit nearest the null
#'
#' Returns one when the confidence interval includes the null rate ratio. For a
#' non-null interval, the endpoint closest to one is used.
#'
#' @param rr Point-estimate rate ratio.
#' @param lower,upper Confidence-interval limits.
#' @return Numeric E-value for the confidence interval.
#' @examples
#' mediational_evalue_ci(rr = 1.8, lower = 1.2, upper = 2.7)
#' mediational_evalue_ci(rr = 1.2, lower = 0.9, upper = 1.6)
#' @export
mediational_evalue_ci <- function(rr, lower, upper) {
    args <- cbind(rr = rr, lower = lower, upper = upper)
    out <- rep(NA_real_, nrow(args))
    valid <- apply(args, 1, function(value) {
        all(is.finite(value)) && all(value > 0) && value[[2]] <= value[[3]]
    })
    crosses_null <- valid & args[, "lower"] <= 1 & args[, "upper"] >= 1
    out[crosses_null] <- 1
    non_null <- valid & !crosses_null
    if (any(non_null)) {
        lower_distance <- abs(log(args[non_null, "lower"]))
        upper_distance <- abs(log(args[non_null, "upper"]))
        nearest <- ifelse(
            lower_distance <= upper_distance,
            args[non_null, "lower"],
            args[non_null, "upper"]
        )
        out[non_null] <- mediational_evalue(nearest)
    }
    out
}

# ---- donor-level composition + pseudobulk -------------------------------
# Cell-type proportions per donor (n_donor x D) and, for a target cell type,
# donor-level pseudobulk counts (genes x n_donor) with library sizes.
.donor_composition <- function(cell_type, donor) {
    tab <- table(donor, cell_type)
    P <- as.matrix(tab) / rowSums(tab)
    P
}

# ---- core: g-computation point estimate for one gene --------------------
# Given the outcome design pieces and pre-drawn mediator matrices, return
# the three interventional means and the RR effects. Mdraw_a / Mdraw_as are
# (n*R) x (D-1) matrices of joint ILR draws under exposure a and a*; Xrep is
# the (n*R) x p covariate matrix (donor covariates repeated over draws).
.gcomp_effects <- function(cf, a, as, Mdraw_a, Mdraw_as, Xrep, ncoord) {
    # linear predictor (offset dropped -> rates); cf is the named coef vector
    lp <- function(aY, Mmat) {
        eta <- cf[["(Intercept)"]] + cf[["A"]] * aY
        for (k in seq_len(ncoord)) {
            eta <- eta + cf[[paste0("M", k)]] * Mmat[, k] +
                cf[[paste0("A:M", k)]] * aY * Mmat[, k]
        }
        if (ncol(Xrep)) {
            for (nm in colnames(Xrep)) eta <- eta + cf[[nm]] * Xrep[, nm]
        }
        eta
    }
    psi_aa   <- mean(exp(lp(a,  Mdraw_a)))    # E[Y_{a,  G_a}]
    psi_aas  <- mean(exp(lp(a,  Mdraw_as)))   # E[Y_{a,  G_a*}]
    psi_asas <- mean(exp(lp(as, Mdraw_as)))   # E[Y_{a*, G_a*}]
    IDE <- psi_aas / psi_asas
    IIE <- psi_aa  / psi_aas
    c(IDE_RR = IDE, IIE_RR = IIE, OE_RR = IDE * IIE)
}

#' Interventional Exposure-Response Decomposition (ERD)
#'
#' Decomposes a donor-level exposure's effect on within-cell-type pseudobulk
#' expression into an interventional direct effect (within-composition) and an
#' interventional indirect effect (cell-composition-driven), on the rate-ratio
#' scale, by counterfactual g-computation with a compositional ILR mediator.
#'
#' This is a developmental function and is not exported. Its outcome is a
#' within-cell-type rate, which cannot identify the contribution of cell-type
#' abundance to a tissue mixture, so its direct and indirect effects and
#' E-values must not be read causally or used for confirmatory claims. It is
#' kept, and tested, so the estimand it was built around can be revisited; call
#' it as \code{exposomeSC:::run_erd_interventional()}.
#'
#' @param scee A \linkS4class{SingleCellExposomeExperiment}.
#' @param exposure Character; column of \code{exposureData(scee)} (donor-level).
#' @param celltype Character; the target cell type whose expression is the
#'   outcome (pseudobulk within this type).
#' @param celltype_col Character; \code{colData} column with cell-type labels.
#' @param sample_col Character; donor id column (default "donor_id").
#' @param covariates Character vector; donor covariate columns in
#'   \code{exposureData} (numeric).
#' @param target_genes Character; genes to test (default: all).
#' @param contrast Numeric length-2, \code{c(a, a*)} exposure levels
#'   (default \code{c(1, 0)}; for continuous exposure use two chosen values).
#' @param n_mc Monte-Carlo mediator draws per donor (default 200).
#' @param n_boot Donor bootstrap replicates for CIs (default 2000).
#' @param interaction Logical; include A x M in the outcome model. The default
#'   is \code{FALSE}; interaction models require substantially more donors.
#' @param min_cells Integer; minimum target-cell count per donor.
#' @param min_donors Integer; minimum eligible donors.
#' @param min_group_donors Integer; for a binary exposure, minimum eligible
#'   donors in each group.
#' @param min_boot_valid Numeric in (0, 1]; minimum fraction of finite bootstrap
#'   estimates required for an interval and p-value.
#' @param adjust Multiplicity adjustment scope. The default \code{"none"}
#'   prevents a per-cell-type call from being mistaken for study-wide FDR.
#'   Use \code{"within_call"} only when this call is the complete hypothesis
#'   family; otherwise combine calls and adjust globally downstream.
#' @param seed Integer RNG seed.
#' @param BPPARAM A \code{BiocParallel} param; use \code{SnowParam} on Windows.
#'
#' @return A data.frame, one row per gene: IDE_RR, IIE_RR, OE_RR with bootstrap
#'   95\% CIs, proportion mediated, per-effect mediational E-values (point and
#'   CI-limit), the closed-form cross-check, finite-resampling p-values, and
#'   optional within-call BH values.
#' @keywords internal
run_erd_interventional <- function(scee, exposure, celltype, celltype_col,
                                    sample_col = "donor_id", covariates = NULL,
                                     target_genes = NULL, contrast = c(1, 0),
                                     n_mc = 200L, n_boot = 2000L,
                                     interaction = FALSE,
                                     min_cells = 10L,
                                     min_donors = 30L,
                                     min_group_donors = 10L,
                                     min_boot_valid = 0.90,
                                     adjust = c("none", "within_call"),
                                     seed = 20260703L,
                                     BPPARAM = NULL) {
    for (p in c("MASS", "SummarizedExperiment", "SingleCellExperiment")) {
        if (!requireNamespace(p, quietly = TRUE)) stop("Package ", p, " required")
    }
    adjust <- match.arg(adjust)
    if (!is.numeric(contrast) || length(contrast) != 2L ||
            any(!is.finite(contrast)) || contrast[[1]] == contrast[[2]]) {
        stop("contrast must contain two distinct finite numeric levels.")
    }
    if (!is.numeric(n_mc) || length(n_mc) != 1L || n_mc < 1 ||
            !is.numeric(n_boot) || length(n_boot) != 1L || n_boot < 2) {
        stop("n_mc must be >= 1 and n_boot must be >= 2.")
    }
    if (!is.numeric(min_boot_valid) || length(min_boot_valid) != 1L ||
            min_boot_valid <= 0 || min_boot_valid > 1) {
        stop("min_boot_valid must lie in (0, 1].")
    }
    a <- contrast[1]; as <- contrast[2]
    cd <- SummarizedExperiment::colData(scee)
    donor <- as.character(cd[[sample_col]])
    ct <- as.character(cd[[celltype_col]])
    exp_mat <- exposureData(scee)

    # ---- donor-level composition mediator (ILR) --------------------------
    P <- .donor_composition(ct, donor)             # n_donor x D
    D <- ncol(P); V <- .ilr_V(D)
    M_all <- .ilr_coords(P, V)                     # n_donor x (D-1)
    ncoord <- D - 1L
    colnames(M_all) <- paste0("M", seq_len(ncoord))
    donors <- rownames(P)

    # ---- donor-level pseudobulk for the target cell type -----------------
    keep <- ct == celltype
    counts <- SummarizedExperiment::assay(scee, "counts")
    if (is.null(target_genes)) target_genes <- rownames(scee)
    missing_genes <- setdiff(target_genes, rownames(scee))
    if (length(missing_genes)) {
        stop("target_genes absent from scee: ",
             paste(utils::head(missing_genes, 5L), collapse = ", "))
    }
    target_genes <- unique(target_genes)
    cmat <- as.matrix(counts[target_genes, keep, drop = FALSE])
    dvec <- donor[keep]
    pb <- vapply(donors, function(d) rowSums(cmat[, dvec == d, drop = FALSE]),
                 numeric(length(target_genes)))            # genes x n_donor
    pb <- matrix(
        pb,
        nrow = length(target_genes),
        dimnames = list(target_genes, donors)
    )
    # The offset is defined from the complete target-cell transcriptome before
    # target_genes is applied. Candidate-set changes must never alter a retained
    # gene's offset or effect estimate.
    cell_library <- Matrix::colSums(counts[, keep, drop = FALSE])
    libsize <- vapply(donors, function(d) {
        sum(cell_library[dvec == d])
    }, numeric(1))
    logoff <- log(pmax(libsize, 1))

    target_cell_counts <- table(factor(dvec, levels = donors))
    eligible <- as.integer(target_cell_counts) >= min_cells
    if (sum(eligible) < min_donors) {
        stop(
            "run_erd_interventional('", celltype, "'): only ",
            sum(eligible), " donor(s) have >= ", min_cells,
            " target cells; need >= ", min_donors, "."
        )
    }
    donors <- donors[eligible]
    M_all <- M_all[eligible, , drop = FALSE]
    pb <- pb[, eligible, drop = FALSE]
    logoff <- logoff[eligible]

    # ---- donor-level exposure + covariates -------------------------------
    A <- exp_mat[donors, exposure]
    X <- if (length(covariates)) {
        as.matrix(exp_mat[donors, covariates, drop = FALSE])
    } else matrix(0, nrow = length(donors), ncol = 0)

    # ---- complete-case on exposure + covariates --------------------------
    # A donor with a missing covariate makes the mediator model's predict()
    # return NA for that donor, which poisons the shared Monte-Carlo mediator
    # draws and turns EVERY gene's IDE/IIE into NaN (the closed form survives
    # because it never touches per-donor draws). Restrict the whole analysis to
    # donors complete on (exposure, covariates) so the mediator draws, the
    # outcome GLM and the donor bootstrap all use one consistent donor set.
    ok <- stats::complete.cases(A, X)
    if (!all(ok)) {
        dropping_format <- paste0("run_erd_interventional('%s'): dropping ",
                                  "%d/%d donor(s) with NA in exposure/",
                                  "covariates (complete-case).")
        message(sprintf(dropping_format, celltype, sum(!ok), length(ok)))
        donors  <- donors[ok]
        A       <- A[ok]
        X       <- X[ok, , drop = FALSE]
        M_all   <- M_all[ok, , drop = FALSE]
        pb      <- pb[, ok, drop = FALSE]
        logoff  <- logoff[ok]
    }
    if (length(donors) < 5L)
        stop("run_erd_interventional('", celltype, "'): fewer than 5 donors ",
             "with complete exposure/covariates; cannot decompose.")
    exposure_levels <- sort(unique(A))
    if (length(exposure_levels) == 2L) {
        group_counts <- table(factor(A, levels = exposure_levels))
        if (any(group_counts < min_group_donors)) {
            stop(
                "run_erd_interventional('", celltype, "'): binary exposure ",
                "group has fewer than ", min_group_donors,
                " eligible donors."
            )
        }
    }

    # one-off closure that computes all-gene effects on a (possibly bootstrap)
    # index set of donors; returns genes x 3 matrix of RRs (+ closed form).
    est_fun <- function(idx, rng_off) {
        Ab <- A[idx]; Xb <- X[idx, , drop = FALSE]; Mb <- M_all[idx, , drop = FALSE]
        logoffb <- logoff[idx]
        # mediator model: M ~ A + X (multivariate Gaussian), one shared fit
        Dm <- data.frame(Ab = Ab, Xb)
        mm_fit <- stats::lm(Mb ~ ., data = Dm)
        etaA <- stats::coef(mm_fit)["Ab", ]                # (D-1) exposure slope
        Sig <- stats::cov(stats::residuals(mm_fit))
        Lchol <- chol(Sig + diag(1e-8, ncoord))
        # pre-draw joint ILR mediator values under a and a* (shared over genes)
        nb <- length(idx)
        mu_pred <- function(level) stats::predict(
            mm_fit, newdata = data.frame(Ab = rep(level, nb), Xb))
        draw <- function(level) {
            mu <- mu_pred(level)                            # nb x (D-1)
            mu <- mu[rep(seq_len(nb), each = n_mc), , drop = FALSE]
            z <- matrix(stats::rnorm(nrow(mu) * ncoord), ncol = ncoord)
            out <- mu + z %*% Lchol
            colnames(out) <- colnames(M_all); out
        }
        Md_a <- draw(a); Md_as <- draw(as)
        Xrep <- Xb[rep(seq_len(nb), each = n_mc), , drop = FALSE]
        # per-gene outcome NB GLM + g-computation
        res <- matrix(NA_real_, nrow = length(target_genes), ncol = 5,
                      dimnames = list(target_genes,
                                      c("IDE_RR", "IIE_RR", "OE_RR",
                                        "cf_IDE", "cf_IIE")))
        odat <- data.frame(A = Ab, Mb, Xb, logoff = logoffb)
        form <- paste("y ~ A +", paste(colnames(Mb), collapse = " + "))
        if (interaction) form <- paste(form, "+",
            paste0("A:", colnames(Mb), collapse = " + "))
        if (ncol(Xb)) form <- paste(form, "+", paste(colnames(Xb), collapse = " + "))
        form <- stats::as.formula(paste(form, "+ offset(logoff)"))
        for (g in seq_along(target_genes)) {
            odat$y <- pb[g, idx]
            fit <- tryCatch(MASS::glm.nb(form, data = odat),
                            error = function(e) NULL)
            if (is.null(fit)) next
            cf <- stats::coef(fit)
            # ensure interaction names exist even if dropped
            for (k in seq_len(ncoord)) {
                nmk <- paste0("A:M", k)
                if (is.na(cf[nmk])) cf[nmk] <- 0
            }
            eff <- .gcomp_effects(cf, a, as, Md_a, Md_as, Xrep, ncoord)
            # closed-form cross-check (additive working model)
            thetaA <- cf[["A"]]
            thetaM <- vapply(seq_len(ncoord), function(k) cf[[paste0("M", k)]], 0)
            res[g, ] <- c(eff, exp(thetaA), exp(sum(thetaM * etaA)))
        }
        res
    }

    .local_rng_scope(seed)
    point <- est_fun(seq_along(donors), 0L)

    ## Release the large single-cell objects BEFORE the bootstrap so the PSOCK
    ## workers receive only the small donor-level data (pb, M_all, A, X, logoff)
    ## captured by est_fun -- never the full counts matrix.
    rm(list = intersect(c("counts", "cmat", "cd", "ct", "donor", "exp_mat",
                          "P", "keep", "dvec", "scee"), ls())); gc()

    # ---- donor bootstrap (PSOCK-parallel) --------------------------------
    if (is.null(BPPARAM)) {
        if (requireNamespace("BiocParallel", quietly = TRUE)) {
            BPPARAM <- BiocParallel::SerialParam(RNGseed = seed)
        } else BPPARAM <- NULL
    }
    boot_one <- function(b) {
        idx <- sample.int(length(donors), replace = TRUE)
        est_fun(idx, b)[, c("IDE_RR", "IIE_RR", "OE_RR"), drop = FALSE]
    }
    boots <- if (!is.null(BPPARAM)) {
        BiocParallel::bplapply(seq_len(n_boot), boot_one, BPPARAM = BPPARAM)
    } else lapply(seq_len(n_boot), boot_one)

    # ---- assemble output with percentile CIs + E-values ------------------
    # Never winsorise the inferential distribution: truncating bootstrap
    # replicates changes interval coverage. A result is non-estimable when too
    # few finite positive replicates remain.
    min_valid <- ceiling(min_boot_valid * n_boot)
    boot_matrix <- function(effname, log_scale = FALSE) {
        arr <- vapply(boots, function(m) {
            value <- m[, effname]
            if (log_scale) log(value) else value
        }, numeric(length(target_genes)))
        matrix(arr, nrow = length(target_genes))
    }
    ci <- function(effname) {
        arr <- boot_matrix(effname)
        t(vapply(seq_len(nrow(arr)), function(g) {
            value <- arr[g, ]
            value <- value[is.finite(value) & value > 0]
            if (length(value) < min_valid) return(c(NA_real_, NA_real_))
            stats::quantile(value, probs = c(0.025, 0.975), names = FALSE)
        }, numeric(2)))
    }
    ide_ci <- ci("IDE_RR"); iie_ci <- ci("IIE_RR")
    # Centred bootstrap test with a finite-resampling +1 correction. This is an
    # approximate test; product-path indirect effects remain non-regular near a
    # double-zero null and require study-level calibration.
    pval <- function(effname, point_value) {
        arr <- boot_matrix(effname, log_scale = TRUE)
        vapply(seq_len(nrow(arr)), function(g) {
            v <- arr[g, ]; v <- v[is.finite(v)]
            observed <- log(point_value[[g]])
            if (length(v) < min_valid || !is.finite(observed)) return(NA_real_)
            null_draw <- v - observed
            (1 + sum(abs(null_draw) >= abs(observed))) / (length(v) + 1)
        }, 0)
    }
    boot_valid <- function(effname) {
        arr <- boot_matrix(effname)
        rowSums(is.finite(arr) & arr > 0)
    }
    out <- data.frame(
        gene = target_genes, celltype = celltype,
        IDE_RR = point[, "IDE_RR"], IDE_lcl = ide_ci[, 1], IDE_ucl = ide_ci[, 2],
        IIE_RR = point[, "IIE_RR"], IIE_lcl = iie_ci[, 1], IIE_ucl = iie_ci[, 2],
        OE_RR = point[, "OE_RR"],
        prop_mediated = (point[, "IDE_RR"] * (point[, "IIE_RR"] - 1)) /
            (point[, "IDE_RR"] * point[, "IIE_RR"] - 1),
        cf_IDE_RR = point[, "cf_IDE"], cf_IIE_RR = point[, "cf_IIE"],
        evalue_IDE = mediational_evalue(point[, "IDE_RR"]),
        evalue_IIE = mediational_evalue(point[, "IIE_RR"]),
        p_IDE = pval("IDE_RR", point[, "IDE_RR"]),
        p_IIE = pval("IIE_RR", point[, "IIE_RR"]),
        boot_valid_IDE = boot_valid("IDE_RR"),
        boot_valid_IIE = boot_valid("IIE_RR"),
        row.names = NULL, stringsAsFactors = FALSE)
    out$evalue_IDE_ci <- mediational_evalue_ci(
        out$IDE_RR, out$IDE_lcl, out$IDE_ucl
    )
    out$evalue_IIE_ci <- mediational_evalue_ci(
        out$IIE_RR, out$IIE_lcl, out$IIE_ucl
    )

    # ---- reliability guard: drop non-estimable / separated decompositions ---
    # A decomposition is reportable only when its POINT estimate is a finite,
    # plausible rate ratio. MASS::glm.nb can, on a sparse pseudobulk gene,
    # (a) fail to converge -- leaving a NA point estimate yet a spurious
    # bootstrap p-value that would surface as a false "significant" hit -- or
    # (b) quasi-separate, giving an absurd rate ratio (|log RR| enormous, e.g.
    # 1e51). Both are numerical artefacts, not effects. NA out the effect, its
    # CI, E-value and p-value so they can never be reported or adjusted. The
    # bound (RR in [1/100, 100]) sits far outside any credible biological rate
    # ratio and cleanly separates genuine fits from separation; adjusted
    # p-values are then computed only over the reportable decompositions.
    rr_bound <- 100
    bad_ide <- !is.finite(out$IDE_RR) | abs(log(out$IDE_RR)) > log(rr_bound)
    bad_iie <- !is.finite(out$IIE_RR) | abs(log(out$IIE_RR)) > log(rr_bound)
    out$IDE_RR[bad_ide] <- NA; out$IDE_lcl[bad_ide] <- NA; out$IDE_ucl[bad_ide] <- NA
    out$evalue_IDE[bad_ide] <- NA; out$evalue_IDE_ci[bad_ide] <- NA
    out$p_IDE[bad_ide] <- NA
    out$IIE_RR[bad_iie] <- NA; out$IIE_lcl[bad_iie] <- NA; out$IIE_ucl[bad_iie] <- NA
    out$evalue_IIE[bad_iie] <- NA; out$evalue_IIE_ci[bad_iie] <- NA
    out$p_IIE[bad_iie] <- NA
    out$OE_RR[bad_ide | bad_iie] <- NA
    out$prop_mediated[bad_ide | bad_iie] <- NA
    non_estimable_format <- paste0("run_erd_interventional('%s'): %d/%d ",
                                   "gene(s) non-estimable or separated ",
                                   "(NA'd, not reported).")
    if (any(bad_ide | bad_iie))
        message(sprintf(non_estimable_format, celltype,
                        sum(bad_ide | bad_iie), nrow(out)))

    if (adjust == "within_call") {
        out$padj_IIE <- stats::p.adjust(out$p_IIE, "BH")
        out$padj_IDE <- stats::p.adjust(out$p_IDE, "BH")
    } else {
        out$padj_IIE <- NA_real_
        out$padj_IDE <- NA_real_
    }
    attr(out, "multiplicity") <- list(
        scope = adjust,
        warning = if (adjust == "within_call") {
            "Adjusted only within this function call."
        } else {
            "No adjustment; combine the pre-specified study-wide family first."
        }
    )
    out
}
