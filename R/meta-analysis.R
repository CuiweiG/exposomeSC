# R/meta-analysis.R
# Cross-cohort meta-analysis of sc-ExWAS results

#' @include AllClasses.R
#' @importFrom stats p.adjust
NULL

#' Meta-analysis across cohorts for sc-ExWAS
#'
#' Combines independent-cohort sc-ExWAS estimates using inverse-variance
#' weighted fixed-effect synthesis or random-effects synthesis. Keys are
#' matched exactly on gene, cell type, and exposure.
#'
#' @param results_list A list of \code{DataFrame}s or data frames. Each cohort
#'   must contain unique \code{gene}-\code{celltype}-\code{exposure} keys and
#'   an effect estimate plus standard error.
#' @param cohort_names Character vector with one unique name per cohort.
#'   Defaults to \code{Cohort_1}, \code{Cohort_2}, and so forth.
#' @param method Character. \code{"fixed"} (default), or \code{"random"} for
#'   REML between-cohort variance with Hartung--Knapp inference. Random-effects
#'   synthesis requires \pkg{metafor} and at least three cohorts per key.
#' @param effect_col Optional character scalar naming the effect column. When
#'   \code{NULL}, the first available column among \code{log2FC},
#'   \code{log2FoldChange}, and \code{estimate} is used separately in each
#'   cohort.
#' @param se_col Optional character scalar naming the standard-error column.
#'   When \code{NULL}, the first available column among \code{se} and
#'   \code{lfcSE} is used separately in each cohort. edgeR QL results from
#'   \code{run_sc_exwas} are deliberately unsupported because
#'   \code{glmQLFTest} does not supply coefficient standard errors.
#'
#' @return A \code{DataFrame} with one row per exact key and columns
#'   \code{meta_effect}, \code{meta_se}, \code{meta_statistic},
#'   \code{meta_df}, \code{meta_pvalue}, \code{meta_padj}, \code{I2},
#'   \code{Q_pvalue}, \code{tau2}, \code{n_cohorts}, \code{cohorts},
#'   \code{direction}, and \code{method}. The retained \code{meta_z} column is
#'   an alias of \code{meta_statistic} for backwards compatibility; under
#'   Hartung--Knapp inference it is a t statistic rather than a z statistic.
#'
#' @details
#' Fixed-effect synthesis uses
#' \deqn{\hat{\beta} = \frac{\sum_i w_i\hat{\beta}_i}{\sum_i w_i},
#'       \quad w_i = 1/SE_i^2.}
#' Random-effects synthesis uses restricted maximum likelihood for
#' \eqn{\tau^2} and Hartung--Knapp small-sample inference. It is deliberately
#' unavailable for only two cohorts. The \eqn{I^2} and Cochran Q statistics are
#' reported as heterogeneity diagnostics; they do not establish transportability
#' or independence of the contributing cohorts.
#'
#' The function never maps gene identifiers or cell-type labels. Harmonisation
#' must occur before calling it. Within-cohort duplicate keys are an error,
#' preventing accidental double weighting after label aggregation.
#' Inverse-variance synthesis also rejects both current edgeR QL results, whose
#' standard errors are explicitly unavailable, and older results containing a
#' finite QL F-ratio proxy. Use a backend that estimates valid coefficient
#' standard errors for each independent cohort before meta-analysis.
#'
#' @references
#' Borenstein M et al. (2009). \emph{Introduction to Meta-Analysis}. Wiley.
#'
#' Viechtbauer W (2010). Conducting meta-analyses in R with the metafor
#' package. \emph{Journal of Statistical Software} 36:1--48.
#' \doi{10.18637/jss.v036.i03}
#'
#' @export
#' @examples
#' library(S4Vectors)
#' r1 <- DataFrame(
#'     gene = c("A", "B"), celltype = "Mono",
#'     log2FC = c(0.5, -0.3), se = c(0.1, 0.15),
#'     exposure = "PM2.5"
#' )
#' r2 <- DataFrame(
#'     gene = c("A", "B", "C"), celltype = "Mono",
#'     log2FC = c(0.4, -0.1, 0.2), se = c(0.12, 0.2, 0.1),
#'     exposure = "PM2.5"
#' )
#' meta <- run_meta_exwas(
#'     list(r1, r2),
#'     cohort_names = c("HELIX", "ENVIRONAGE")
#' )
#' meta
run_meta_exwas <- function(results_list,
                            cohort_names = NULL,
                            method = c("fixed", "random"),
                            effect_col = NULL,
                            se_col = NULL) {
    method <- match.arg(method)

    if (!is.list(results_list) || length(results_list) < 2L) {
        stop("results_list must be a list containing at least two cohorts.")
    }
    if (!is.null(effect_col) &&
            (!is.character(effect_col) || length(effect_col) != 1L ||
                is.na(effect_col) || !nzchar(effect_col))) {
        stop("effect_col must be NULL or one non-missing column name.")
    }
    if (!is.null(se_col) &&
            (!is.character(se_col) || length(se_col) != 1L ||
                is.na(se_col) || !nzchar(se_col))) {
        stop("se_col must be NULL or one non-missing column name.")
    }

    if (is.null(cohort_names)) {
        cohort_names <- paste0("Cohort_", seq_along(results_list))
    }
    if (!is.character(cohort_names) ||
            length(cohort_names) != length(results_list) ||
            anyNA(cohort_names) || any(!nzchar(cohort_names)) ||
            anyDuplicated(cohort_names)) {
        stop(
            "cohort_names must contain one unique, non-missing name per cohort."
        )
    }
    if (method == "random" &&
            !requireNamespace("metafor", quietly = TRUE)) {
        stop("Package 'metafor' is required for random-effects synthesis.")
    }

    resolve_column <- function(data, requested, candidates, role, cohort) {
        if (!is.null(requested)) {
            if (!requested %in% colnames(data)) {
                stop(
                    "Column '", requested, "' not found for ", role,
                    " in cohort '", cohort, "'."
                )
            }
            return(requested)
        }
        available <- candidates[candidates %in% colnames(data)]
        if (!length(available)) {
            stop(
                "No supported ", role, " column found in cohort '", cohort,
                "'. Tried: ", paste(candidates, collapse = ", "), "."
            )
        }
        available[[1]]
    }

    all_data <- vector("list", length(results_list))
    for (i in seq_along(results_list)) {
        cohort <- cohort_names[[i]]
        data <- as.data.frame(results_list[[i]])
        edgeR_ql_method <- if ("method" %in% colnames(data)) {
            method_value <- as.character(data$method)
            !is.na(method_value) & grepl(
                "^edgeR(_robust)?_QL$",
                method_value,
                ignore.case = TRUE
            )
        } else {
            rep(FALSE, nrow(data))
        }
        uncertainty_parameters <- tryCatch(
            S4Vectors::metadata(results_list[[i]])$parameters,
            error = function(error) NULL
        )
        edgeR_ql_contract <- is.list(uncertainty_parameters) &&
            identical(
                uncertainty_parameters$se_method,
                "not_available_edgeR_QL"
            )
        if (any(edgeR_ql_method) || edgeR_ql_contract) {
            stop(
                "edgeR quasi-likelihood output does not supply coefficient ",
                "standard errors and cannot enter inverse-variance synthesis ",
                "for cohort '", cohort, "'. Use independently fitted ",
                "uncertainty from a method that supplies valid coefficient ",
                "standard errors."
            )
        }
        missing <- setdiff(c("gene", "celltype"), colnames(data))
        if (length(missing)) {
            stop(
                "Missing key column(s) in cohort '", cohort, "': ",
                paste(missing, collapse = ", "), "."
            )
        }
        gene <- as.character(data$gene)
        celltype <- as.character(data$celltype)
        if (anyNA(gene) || any(!nzchar(gene)) ||
                anyNA(celltype) || any(!nzchar(celltype))) {
            stop("gene and celltype keys must be complete and non-empty.")
        }

        effect_name <- resolve_column(
            data,
            effect_col,
            c("log2FC", "log2FoldChange", "estimate"),
            "effect",
            cohort
        )
        se_name <- resolve_column(
            data,
            se_col,
            c("se", "lfcSE"),
            "standard-error",
            cohort
        )
        if (!is.numeric(data[[effect_name]]) ||
                !is.numeric(data[[se_name]])) {
            stop(
                "Effect and standard-error columns must be numeric in cohort '",
                cohort, "'."
            )
        }
        effect <- as.numeric(data[[effect_name]])
        standard_error <- as.numeric(data[[se_name]])
        if (length(standard_error) && all(is.na(standard_error))) {
            stop(
                "No coefficient standard errors are available in cohort '",
                cohort, "'; inverse-variance synthesis is undefined."
            )
        }
        if (any(!is.na(effect) & !is.finite(effect))) {
            stop("Non-finite effect estimate in cohort '", cohort, "'.")
        }
        if (any(!is.na(standard_error) &
                (!is.finite(standard_error) | standard_error <= 0))) {
            stop(
                "Standard errors must be positive and finite in cohort '",
                cohort, "'."
            )
        }

        exposure <- if ("exposure" %in% colnames(data)) {
            as.character(data$exposure)
        } else {
            rep(NA_character_, nrow(data))
        }
        exposure_key <- ifelse(
            is.na(exposure),
            "<unspecified>",
            exposure
        )
        key_data <- data.frame(
            gene = gene,
            celltype = celltype,
            exposure = exposure_key,
            stringsAsFactors = FALSE
        )
        if (anyDuplicated(key_data)) {
            duplicate <- which(
                duplicated(key_data) |
                    duplicated(key_data, fromLast = TRUE)
            )[[1]]
            stop(
                "Duplicate gene-celltype-exposure key in cohort '", cohort,
                "': ", paste(key_data[duplicate, ], collapse = " | "), "."
            )
        }

        all_data[[i]] <- data.frame(
            gene = gene,
            celltype = celltype,
            .exposure_key = exposure_key,
            exposure = exposure,
            .effect = effect,
            .se = standard_error,
            cohort = cohort,
            stringsAsFactors = FALSE
        )
    }

    combined <- do.call(rbind, all_data)
    keys <- unique(combined[
        ,
        c("gene", "celltype", ".exposure_key"),
        drop = FALSE
    ])

    results <- lapply(seq_len(nrow(keys)), function(key_index) {
        key <- keys[key_index, , drop = FALSE]
        selected <- combined$gene == key$gene &
            combined$celltype == key$celltype &
            combined$.exposure_key == key$.exposure_key
        data <- combined[selected, , drop = FALSE]
        complete <- !is.na(data$.effect) & !is.na(data$.se)
        data <- data[complete, , drop = FALSE]
        data <- data[
            order(match(data$cohort, cohort_names)),
            ,
            drop = FALSE
        ]

        if (nrow(data) < 2L) return(NULL)
        if (method == "random" && nrow(data) < 3L) {
            stop(
                "Random-effects synthesis requires at least three cohorts for ",
                key$gene, " | ", key$celltype, " | ",
                key$.exposure_key, "."
            )
        }

        effects <- data$.effect
        standard_errors <- data$.se
        n_cohorts <- length(effects)
        weights <- 1 / standard_errors^2
        fixed_effect <- sum(weights * effects) / sum(weights)
        fixed_se <- sqrt(1 / sum(weights))
        Q <- sum(weights * (effects - fixed_effect)^2)
        Q_df <- n_cohorts - 1L
        Q_pvalue <- stats::pchisq(Q, df = Q_df, lower.tail = FALSE)
        I2 <- if (Q > 0) max(0, (Q - Q_df) / Q) else 0

        if (method == "random") {
            fit <- metafor::rma.uni(
                yi = effects,
                sei = standard_errors,
                method = "REML",
                test = "knha"
            )
            meta_effect <- unname(stats::coef(fit)[[1]])
            meta_se <- unname(fit$se[[1]])
            meta_statistic <- unname(fit$zval[[1]])
            meta_df <- n_cohorts - 1L
            meta_pvalue <- unname(fit$pval[[1]])
            tau2 <- unname(fit$tau2)
            method_label <- "REML_Hartung-Knapp"
        } else {
            meta_effect <- fixed_effect
            meta_se <- fixed_se
            meta_statistic <- meta_effect / meta_se
            meta_df <- Inf
            meta_pvalue <- 2 * stats::pnorm(
                abs(meta_statistic),
                lower.tail = FALSE
            )
            tau2 <- 0
            method_label <- "inverse_variance_fixed"
        }
        meta_pvalue <- pmax(meta_pvalue, .Machine$double.xmin)
        signs <- ifelse(
            effects > 0,
            "+",
            ifelse(effects < 0, "-", "0")
        )

        data.frame(
            gene = key$gene,
            celltype = key$celltype,
            meta_effect = meta_effect,
            meta_se = meta_se,
            meta_statistic = meta_statistic,
            meta_z = meta_statistic,
            meta_df = meta_df,
            meta_pvalue = meta_pvalue,
            I2 = I2,
            Q_pvalue = Q_pvalue,
            tau2 = tau2,
            n_cohorts = n_cohorts,
            exposure = if (key$.exposure_key == "<unspecified>") {
                NA_character_
            } else {
                key$.exposure_key
            },
            cohorts = paste(data$cohort, collapse = ";"),
            direction = paste(signs, collapse = ""),
            method = method_label,
            stringsAsFactors = FALSE
        )
    })

    results <- Filter(Negate(is.null), results)
    if (!length(results)) {
        return(S4Vectors::DataFrame(
            gene = character(),
            celltype = character(),
            meta_effect = numeric(),
            meta_se = numeric(),
            meta_statistic = numeric(),
            meta_z = numeric(),
            meta_df = numeric(),
            meta_pvalue = numeric(),
            meta_padj = numeric(),
            I2 = numeric(),
            Q_pvalue = numeric(),
            tau2 = numeric(),
            n_cohorts = integer(),
            exposure = character(),
            cohorts = character(),
            direction = character(),
            method = character()
        ))
    }

    output <- do.call(rbind, results)
    output$meta_padj <- stats::p.adjust(output$meta_pvalue, method = "BH")
    output <- output[
        order(output$meta_padj, output$meta_pvalue),
        ,
        drop = FALSE
    ]
    rownames(output) <- NULL
    S4Vectors::DataFrame(output)
}
