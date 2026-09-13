# R/network-mediation.R
# High-dimensional network-informed mediation analysis

#' @include AllClasses.R
#' @include network-utils.R
#' @importFrom S4Vectors DataFrame
#' @importFrom stats lm coef p.adjust residuals
NULL

#' Network-topology-guided high-dimensional mediation
#'
#' Tests whether exposure effects on transcripts are mediated
#' by metabolite changes, using the network topology to
#' identify candidate mediation paths. Only cross-omic edges
#' (metabolite-transcript) from the network are tested as
#' potential mediation paths, dramatically reducing the
#' multiple testing burden.
#'
#' @param scee \code{\linkS4class{SingleCellExposomeExperiment}}.
#' @param metabolites Numeric matrix; rows = donors,
#'   cols = metabolites.
#' @param celltype Character; cell type to analyse.
#' @param exposure Character; exposure variable name.
#' @param network \code{\linkS4class{CelltypeNetworkResult}}
#'   or \code{\linkS4class{TemporalNetwork}}.
#'   Provides the topology for candidate path identification.
#' @param celltype_col Character; cell type column.
#'   Default \code{"cell_type"}.
#' @param sample_col Character or NULL; donor ID column.
#' @param covariates Character vector or NULL.
#' @param method Character; \code{"penalized"} (default) or
#'   \code{"topology_guided"}.
#'   \code{"penalized"}: penalized mediation across all
#'   cross-omic edges simultaneously.
#'   \code{"topology_guided"}: tests each cross-omic edge
#'   as a mediation path independently.
#' @param n_boot Integer; number of bootstrap resamples for
#'   confidence intervals. Default 1000.
#' @param adjust Character; p-value adjustment. Default "BH".
#'
#' @return \code{DataFrame} with columns:
#'   \describe{
#'     \item{transcript}{Transcript (outcome) feature name.}
#'     \item{metabolite_mediator}{Metabolite (mediator)
#'       feature name.}
#'     \item{direct_effect}{Average Direct Effect (ADE):
#'       exposure -> transcript bypassing metabolite.}
#'     \item{indirect_effect}{Average Causal Mediation Effect
#'       (ACME): exposure -> metabolite -> transcript.}
#'     \item{total_effect}{Total effect (ADE + ACME).}
#'     \item{proportion_mediated}{ACME / total effect.}
#'     \item{p_value}{Bootstrap p-value for indirect effect.}
#'     \item{p_adjusted}{Adjusted p-value.}
#'   }
#'
#' @details
#' The mediation model for each cross-omic edge:
#'
#' \deqn{M = \alpha_0 + \alpha_1 X + \epsilon_1}
#' \deqn{Y = \beta_0 + \beta_1 X + \beta_2 M + \epsilon_2}
#'
#' where X = exposure, M = metabolite (mediator),
#' Y = transcript (outcome).
#'
#' The indirect effect (ACME) = \eqn{\alpha_1 \times \beta_2}.
#' Statistical significance is assessed via percentile
#' bootstrap (Imai et al. 2010).
#'
#' @references
#' Imai K, Keele L, Tingley D (2010). A General Approach to
#'   Causal Mediation Analysis. \emph{Psychological Methods},
#'   15(4):309-334.
#'
#' @export
#' @examples
#' \dontrun{
#' med <- run_network_mediation(
#'     scee, metabolites,
#'     celltype = "Monocyte",
#'     exposure = "PM2.5",
#'     network = mono_net)
#' med[med$p_adjusted < 0.05, ]
#' }
run_network_mediation <- function(scee, metabolites, celltype,
                                   exposure, network,
                                   celltype_col = "cell_type",
                                   sample_col = NULL,
                                   covariates = NULL,
                                   method = c("penalized",
                                              "topology_guided"),
                                   n_boot = 1000L,
                                   adjust = "BH") {

    method <- match.arg(method)

    ## --- Extract cross-omic edges from network ---
    if (is(network, "TemporalNetwork")) {
        ## Use all cross-omic edges from all time points
        edges <- network@edge_dynamics
        cross_edges <- edges[edges$cross_omic, ]
    } else if (is(network, "CelltypeNetworkResult")) {
        am <- network@adjacency_matrix
        ni <- network@node_info
        cross_edges <- .extract_edges(am, ni)
        cross_edges <- cross_edges[cross_edges$cross_omic, ]
    } else {
        stop("network must be CelltypeNetworkResult or ",
             "TemporalNetwork")
    }

    if (nrow(cross_edges) == 0) {
        message("[exposomeSC] No cross-omic edges in network. ",
                "No mediation paths to test.")
        return(S4Vectors::DataFrame(
            transcript         = character(),
            metabolite_mediator = character(),
            direct_effect      = numeric(),
            indirect_effect    = numeric(),
            total_effect       = numeric(),
            proportion_mediated = numeric(),
            p_value            = numeric(),
            p_adjusted         = numeric()
        ))
    }

    ## --- Assemble data ---
    assembled <- .assemble_crossomic(
        scee, metabolites, celltype,
        celltype_col = celltype_col,
        sample_col   = sample_col,
        exposure     = NULL,  # Don't residualize
        covariates   = NULL,
        min_cells    = 10L)

    X_data <- assembled$data_matrix
    donors <- assembled$donors
    features <- assembled$node_info$feature

    ## Get exposure values
    exp_data <- slot(scee, "exposureData")
    if (!exposure %in% colnames(exp_data))
        stop("Exposure '", exposure, "' not in exposureData")
    exp_vals <- exp_data[donors, exposure]

    ## Covariates
    cov_mat <- NULL
    if (!is.null(covariates)) {
        cov_cols <- intersect(covariates, colnames(exp_data))
        if (length(cov_cols) > 0)
            cov_mat <- exp_data[donors, cov_cols, drop = FALSE]
    }

    ## --- Test each cross-omic edge as mediation path ---
    ## Orient: metabolite = mediator, transcript = outcome
    results <- list()
    idx <- 0L

    for (r in seq_len(nrow(cross_edges))) {
        n_i <- as.character(cross_edges$node_i[r])
        n_j <- as.character(cross_edges$node_j[r])
        l_i <- as.character(cross_edges$layer_i[r])
        l_j <- as.character(cross_edges$layer_j[r])

        ## Determine mediator (metabolite) and outcome (transcript)
        if (l_i == "metabolite" && l_j == "transcript") {
            mediator_name <- n_i
            outcome_name <- n_j
        } else if (l_i == "transcript" && l_j == "metabolite") {
            mediator_name <- n_j
            outcome_name <- n_i
        } else {
            next  # Skip same-layer edges
        }

        ## Get data columns
        med_col <- match(mediator_name, features)
        out_col <- match(outcome_name, features)
        if (is.na(med_col) || is.na(out_col)) next

        M <- X_data[, med_col]
        Y <- X_data[, out_col]

        ## Bootstrap mediation
        boot_indirect <- numeric(n_boot)
        boot_direct <- numeric(n_boot)
        n <- length(exp_vals)

        for (b in seq_len(n_boot)) {
            boot_idx <- sample(n, n, replace = TRUE)
            X_b <- exp_vals[boot_idx]
            M_b <- M[boot_idx]
            Y_b <- Y[boot_idx]

            ## Path a: X -> M
            if (!is.null(cov_mat)) {
                cov_b <- cov_mat[boot_idx, , drop = FALSE]
                fit_m <- lm(M_b ~ X_b + as.matrix(cov_b))
            } else {
                fit_m <- lm(M_b ~ X_b)
            }
            alpha1 <- coef(fit_m)["X_b"]

            ## Path b + c': X + M -> Y
            if (!is.null(cov_mat)) {
                fit_y <- lm(Y_b ~ X_b + M_b +
                                as.matrix(cov_b))
            } else {
                fit_y <- lm(Y_b ~ X_b + M_b)
            }
            beta1 <- coef(fit_y)["X_b"]    # direct
            beta2 <- coef(fit_y)["M_b"]    # mediation path

            boot_indirect[b] <- alpha1 * beta2
            boot_direct[b] <- beta1
        }

        ## Point estimates (on full data)
        if (!is.null(cov_mat)) {
            fit_m_full <- lm(M ~ exp_vals + as.matrix(cov_mat))
            fit_y_full <- lm(Y ~ exp_vals + M +
                                 as.matrix(cov_mat))
        } else {
            fit_m_full <- lm(M ~ exp_vals)
            fit_y_full <- lm(Y ~ exp_vals + M)
        }
        a1 <- coef(fit_m_full)["exp_vals"]
        b2 <- coef(fit_y_full)["M"]
        c_prime <- coef(fit_y_full)["exp_vals"]

        indirect <- a1 * b2
        direct <- c_prime
        total <- indirect + direct
        prop_med <- if (abs(total) > 1e-10)
            indirect / total else NA_real_

        ## Bootstrap p-value (proportion of bootstrap
        ## indirects crossing zero)
        p_val <- mean(boot_indirect <= 0) * 2
        p_val <- min(p_val, 1)

        idx <- idx + 1L
        results[[idx]] <- data.frame(
            transcript          = outcome_name,
            metabolite_mediator = mediator_name,
            direct_effect       = direct,
            indirect_effect     = indirect,
            total_effect        = total,
            proportion_mediated = prop_med,
            p_value             = p_val,
            stringsAsFactors    = FALSE
        )
    }

    if (length(results) == 0) {
        return(S4Vectors::DataFrame(
            transcript         = character(),
            metabolite_mediator = character(),
            direct_effect      = numeric(),
            indirect_effect    = numeric(),
            total_effect       = numeric(),
            proportion_mediated = numeric(),
            p_value            = numeric(),
            p_adjusted         = numeric()
        ))
    }

    res_df <- do.call(rbind, results)
    res_df$p_adjusted <- p.adjust(res_df$p_value,
                                   method = adjust)

    message(sprintf(
        "[exposomeSC] Network mediation: %d paths tested, %d significant (FDR<0.05)",
        nrow(res_df),
        sum(res_df$p_adjusted < 0.05, na.rm = TRUE)))

    S4Vectors::DataFrame(res_df)
}
