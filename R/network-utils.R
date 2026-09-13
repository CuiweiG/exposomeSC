# R/network-utils.R
# Internal helpers for cross-omic network inference

#' @include AllClasses.R
#' @include utils.R
#' @importFrom S4Vectors DataFrame
#' @importFrom SummarizedExperiment assay colData
#' @importFrom stats lm residuals var cor
NULL

# -------------------------------------------------------
# Pseudobulk + metabolite block assembly
# -------------------------------------------------------

#' Assemble cross-omic data matrix for network estimation
#'
#' Internal. Creates a donors x (transcripts + metabolites)
#' data matrix from a SCEE and metabolite matrix for a
#' specified cell type.
#'
#' @param scee SingleCellExposomeExperiment
#' @param metabolites matrix; rows = donors, cols = metabolites
#' @param celltype character; cell type to aggregate
#' @param celltype_col character; colData column for cell types
#' @param sample_col character; colData column for donor IDs
#' @param exposure character or NULL; exposure to residualize
#' @param covariates character vector or NULL
#' @param min_cells integer; minimum cells per pseudobulk
#' @param top_var_genes integer or NULL; pre-select top
#'   variable genes to reduce dimensionality
#' @param vst logical; use DESeq2 VST (TRUE) or log-CPM (FALSE)
#'
#' @return list with:
#'   - data_matrix: donors x features numeric matrix
#'   - node_info: DataFrame mapping features to omic layers
#'   - donors: character vector of donor IDs used
#'   - n_transcripts: integer
#'   - n_metabolites: integer
#'
#' @keywords internal
#' @noRd
.assemble_crossomic <- function(scee, metabolites, celltype,
                                 celltype_col = "cell_type",
                                 sample_col = NULL,
                                 exposure = NULL,
                                 covariates = NULL,
                                 min_cells = 10L,
                                 top_var_genes = NULL,
                                 vst = TRUE) {
    ## --- Resolve sample column ---
    if (is.null(sample_col)) {
        sm <- slot(scee, "sampleMap")
        if (nrow(sm) > 0) {
            sample_col <- "sample_id"
        } else {
            stop("Cannot determine sample column. ",
                 "Provide sample_col or use build_scee().")
        }
    }

    ## --- Get cell type annotations ---
    cd <- colData(scee)
    if (!celltype_col %in% colnames(cd))
        stop("'", celltype_col, "' not in colData")

    cell_types <- as.character(cd[[celltype_col]])

    ## --- Resolve sample IDs per cell ---
    sm <- slot(scee, "sampleMap")
    if (nrow(sm) > 0 && sample_col == "sample_id") {
        cell_to_sample <- setNames(
            as.character(sm$sample_id),
            as.character(sm$cell_id))
        samples <- cell_to_sample[colnames(scee)]
    } else if (sample_col %in% colnames(cd)) {
        samples <- as.character(cd[[sample_col]])
    } else {
        stop("'", sample_col, "' not in colData")
    }

    ## --- Pseudobulk aggregation ---
    counts_mat <- assay(scee, "counts")
    pb <- .pseudobulk_aggregate(
        counts_mat, samples, cell_types,
        celltype, min_cells = min_cells)

    if (is.null(pb))
        stop("No valid pseudobulk samples for cell type '",
             celltype, "' (min_cells=", min_cells, ")")

    ## --- VST or log-CPM transform ---
    tx_mat <- NULL
    if (vst && requireNamespace("DESeq2", quietly = TRUE)) {
        tx_mat <- tryCatch({
            suppressMessages({
                dds <- DESeq2::DESeqDataSetFromMatrix(
                    countData = pb$pb_mat,
                    colData = data.frame(
                        donor = pb$valid_donors),
                    design = ~ 1)
                vsd <- DESeq2::vst(dds, blind = TRUE)
                t(SummarizedExperiment::assay(vsd))
            })
        }, error = function(e) NULL)
    }
    ## log-CPM when VST is not requested, DESeq2 is absent or VST fails
    ## (for example with fewer genes than its subsampling needs)
    transform <- if (is.null(tx_mat)) "log_cpm" else "vst"
    if (is.null(tx_mat)) tx_mat <- t(.log_cpm(pb$pb_mat))
    ## tx_mat: donors x genes

    ## --- Optional: select top variable genes ---
    if (!is.null(top_var_genes) && ncol(tx_mat) > top_var_genes) {
        gene_vars <- apply(tx_mat, 2, var)
        keep <- order(gene_vars, decreasing = TRUE)[
            seq_len(top_var_genes)]
        tx_mat <- tx_mat[, keep, drop = FALSE]
    }

    ## --- Match donors between transcripts and metabolites ---
    tx_donors <- rownames(tx_mat)
    met_donors <- rownames(metabolites)
    if (is.null(met_donors))
        stop("metabolites matrix must have row names (donor IDs)")

    shared_donors <- intersect(tx_donors, met_donors)
    if (length(shared_donors) < 5L)
        stop("Only ", length(shared_donors),
             " shared donors between pseudobulk and ",
             "metabolites. Need at least 5.")

    tx_mat <- tx_mat[shared_donors, , drop = FALSE]
    met_mat <- metabolites[shared_donors, , drop = FALSE]

    ## --- Residualization (optional) ---
    if (!is.null(exposure) || !is.null(covariates)) {
        exp_data <- slot(scee, "exposureData")
        adj_vars <- c(exposure, covariates)
        adj_vars <- adj_vars[adj_vars %in% colnames(exp_data)]

        if (length(adj_vars) > 0) {
            adj_mat <- exp_data[shared_donors, adj_vars,
                                drop = FALSE]
            tx_mat <- .residualize_matrix(tx_mat, adj_mat)
            met_mat <- .residualize_matrix(met_mat, adj_mat)
        }
    }

    ## --- Assemble block matrix ---
    ## Scale columns to unit variance for balanced penalization
    tx_scaled <- scale(tx_mat)
    met_scaled <- scale(met_mat)

    ## Replace NaN from zero-variance columns
    tx_scaled[is.nan(tx_scaled)] <- 0
    met_scaled[is.nan(met_scaled)] <- 0

    data_matrix <- cbind(tx_scaled, met_scaled)

    ## --- Node info ---
    tx_names <- colnames(tx_mat)
    met_names <- colnames(met_mat)
    if (is.null(tx_names))
        tx_names <- paste0("gene_", seq_len(ncol(tx_mat)))
    if (is.null(met_names))
        met_names <- paste0("metab_", seq_len(ncol(met_mat)))

    node_info <- S4Vectors::DataFrame(
        feature    = c(tx_names, met_names),
        omic_layer = c(rep("transcript", length(tx_names)),
                       rep("metabolite", length(met_names))),
        block      = c(rep(1L, length(tx_names)),
                       rep(2L, length(met_names)))
    )

    list(
        data_matrix    = data_matrix,
        node_info      = node_info,
        donors         = shared_donors,
        transform      = transform,
        n_transcripts  = length(tx_names),
        n_metabolites  = length(met_names)
    )
}

# -------------------------------------------------------
# Residualization helper
# -------------------------------------------------------

#' Residualize a matrix against adjustment variables
#'
#' For each column in \code{Y}, fits a linear model
#' \code{Y[,j] ~ adj} and returns residuals.
#'
#' @param Y numeric matrix (n x p)
#' @param adj numeric matrix (n x q) of adjustment variables
#' @return matrix of residuals (n x p)
#'
#' @keywords internal
#' @noRd
.residualize_matrix <- function(Y, adj) {
    adj_df <- as.data.frame(adj)
    result <- vapply(seq_len(ncol(Y)), function(j) {
        fit <- lm(Y[, j] ~ ., data = adj_df)
        residuals(fit)
    }, numeric(nrow(Y)))
    dimnames(result) <- dimnames(Y)
    result
}

# -------------------------------------------------------
# Block indicator vector for coglasso
# -------------------------------------------------------

#' Create block indicator for coglasso
#'
#' @param n_transcripts integer
#' @param n_metabolites integer
#' @return integer vector: 1 for transcript, 2 for metabolite
#'
#' @keywords internal
#' @noRd
.make_block_indicator <- function(n_transcripts, n_metabolites) {
    c(rep(1L, n_transcripts), rep(2L, n_metabolites))
}

# -------------------------------------------------------
# Extract adjacency from precision matrix
# -------------------------------------------------------

#' Convert precision matrix to binary adjacency
#'
#' @param precision numeric matrix (p x p)
#' @param threshold numeric; absolute value threshold
#' @return binary matrix (p x p)
#'
#' @keywords internal
#' @noRd
.precision_to_adjacency <- function(precision,
                                     threshold = 1e-10) {
    adj <- (abs(precision) > threshold) * 1L
    diag(adj) <- 0L
    adj
}

# -------------------------------------------------------
# Count cross-omic edges
# -------------------------------------------------------

#' Count edges between different omic layers
#'
#' @param adjacency binary matrix
#' @param node_info DataFrame with omic_layer column
#' @return list with total_edges, within_transcript,
#'   within_metabolite, cross_omic counts
#'
#' @keywords internal
#' @noRd
.count_edge_types <- function(adjacency, node_info) {
    layers <- node_info$omic_layer
    p <- nrow(adjacency)
    total <- within_tx <- within_met <- cross <- 0L

    for (i in seq_len(p - 1L)) {
        for (j in (i + 1L):p) {
            if (adjacency[i, j] != 0) {
                total <- total + 1L
                if (layers[i] == layers[j]) {
                    if (layers[i] == "transcript")
                        within_tx <- within_tx + 1L
                    else
                        within_met <- within_met + 1L
                } else {
                    cross <- cross + 1L
                }
            }
        }
    }

    list(total = total,
         within_transcript  = within_tx,
         within_metabolite  = within_met,
         cross_omic         = cross)
}

# -------------------------------------------------------
# Edge list extraction
# -------------------------------------------------------

#' Extract edge list from adjacency matrix
#'
#' @param adjacency binary matrix
#' @param node_info DataFrame
#' @param precision optional precision matrix for weights
#' @param stability optional stability matrix
#' @return DataFrame with node_i, node_j, weight, stability,
#'   layer_i, layer_j, cross_omic columns
#'
#' @keywords internal
#' @noRd
.extract_edges <- function(adjacency, node_info,
                            precision = NULL,
                            stability = NULL) {
    p <- nrow(adjacency)
    edges <- list()
    idx <- 0L

    for (i in seq_len(p - 1L)) {
        for (j in (i + 1L):p) {
            if (adjacency[i, j] != 0) {
                idx <- idx + 1L
                edges[[idx]] <- list(
                    node_i     = node_info$feature[i],
                    node_j     = node_info$feature[j],
                    weight     = if (!is.null(precision))
                        precision[i, j] else NA_real_,
                    stability  = if (!is.null(stability))
                        stability[i, j] else NA_real_,
                    layer_i    = node_info$omic_layer[i],
                    layer_j    = node_info$omic_layer[j],
                    cross_omic = node_info$omic_layer[i] !=
                        node_info$omic_layer[j]
                )
            }
        }
    }

    if (length(edges) == 0) {
        return(S4Vectors::DataFrame(
            node_i     = character(),
            node_j     = character(),
            weight     = numeric(),
            stability  = numeric(),
            layer_i    = character(),
            layer_j    = character(),
            cross_omic = logical()
        ))
    }

    S4Vectors::DataFrame(do.call(rbind.data.frame,
        lapply(edges, as.data.frame,
               stringsAsFactors = FALSE)))
}

# -------------------------------------------------------
# Simple block graphical lasso (fallback when coglasso
# is not available)
# -------------------------------------------------------

#' Block graphical lasso via glasso with separate penalties
#'
#' When coglasso is not installed, estimates the precision
#' matrix using separate lambda values for within-block and
#' between-block entries. Requires the glasso package.
#'
#' @param S sample covariance matrix (p x p)
#' @param blocks integer vector of block assignments
#' @param lambda_w numeric; within-block penalty
#' @param lambda_b numeric; between-block penalty
#' @return list with wi (precision), adj (adjacency)
#'
#' @keywords internal
#' @noRd
.block_glasso <- function(S, blocks, lambda_w = 0.3,
                           lambda_b = 0.5) {
    p <- nrow(S)
    rho_mat <- matrix(lambda_b, p, p)
    for (b in unique(blocks)) {
        idx <- which(blocks == b)
        rho_mat[idx, idx] <- lambda_w
    }
    diag(rho_mat) <- 0

    if (!requireNamespace("stats", quietly = TRUE))
        stop("stats package required")

    ## Use a simple thresholded inverse
    ## (true glasso would need the glasso package)
    ## This is a conservative fallback
    S_reg <- S + diag(max(lambda_w, lambda_b), p)
    wi <- tryCatch(
        solve(S_reg),
        error = function(e) {
            diag(1 / diag(S_reg))
        }
    )

    ## Threshold off-diagonal by rho_mat
    for (i in seq_len(p)) {
        for (j in seq_len(p)) {
            if (i != j && abs(wi[i, j]) < rho_mat[i, j])
                wi[i, j] <- 0
        }
    }

    adj <- .precision_to_adjacency(wi)
    list(wi = wi, adj = adj)
}

# -------------------------------------------------------
# Base selector for stability selection: 0/1 indicator of the q
# features with the largest absolute Spearman correlation with y.
# Constant features are never selected.
# -------------------------------------------------------
#' @keywords internal
.top_marginal_features <- function(x, y, q) {
    chosen <- numeric(ncol(x))
    if (!isTRUE(stats::sd(y) > 0)) return(chosen)
    varies <- which(apply(x, 2L, function(v) isTRUE(stats::sd(v) > 0)))
    if (length(varies) == 0L) return(chosen)
    rho <- abs(stats::cor(x[, varies, drop = FALSE], y,
                          method = "spearman"))[, 1]
    ranked <- varies[order(-rho, varies)]
    chosen[ranked[seq_len(min(q, length(varies)))]] <- 1
    chosen
}
