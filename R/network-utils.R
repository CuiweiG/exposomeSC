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
#' @param x SingleCellExposomeExperiment; each cell's donor is read from
#'   its sample column
#' @param metabolites matrix; rows = donors, cols = metabolites
#' @param celltype character; cell type to aggregate
#' @param celltype_col character; colData column for cell types
#' @param exposure character or NULL; exposure to residualize
#' @param covariates character vector or NULL
#' @param min_cells integer; minimum cells per pseudobulk
#' @param top_var_genes integer or NULL; pre-select top
#'   variable genes to reduce dimensionality
#' @param vst logical; use DESeq2 VST (TRUE) or log-CPM (FALSE)
#' @param drop_constant logical; leave features that are constant across
#'   the donors out of the matrix (TRUE), or keep them as constant columns
#'   (FALSE)
#'
#' @return list with:
#'   - data_matrix: donors x features numeric matrix
#'   - node_info: DataFrame mapping features to omic layers
#'   - donors: character vector of donor IDs used
#'   - n_transcripts: integer
#'   - n_metabolites: integer
#'   - constant_features: character; features left out as constant
#'
#' @keywords internal
#' @noRd
.assemble_crossomic <- function(x, metabolites, celltype,
                                 celltype_col = "cell_type",
                                 exposure = NULL,
                                 covariates = NULL,
                                 min_cells = 10L,
                                 top_var_genes = NULL,
                                 vst = TRUE, drop_constant = TRUE) {
    ## --- Adjustment variables, checked before any computation ---
    .check_adjustment_names(x, c(exposure, covariates))

    ## --- Get cell type annotations ---
    cd <- colData(x)
    if (!celltype_col %in% colnames(cd))
        stop("'", celltype_col, "' not in colData")

    cell_types <- as.character(cd[[celltype_col]])

    ## --- Sample (donor) ID of each cell ---
    samples <- .cell_samples(x)

    ## --- Pseudobulk aggregation ---
    counts_mat <- assay(x, "counts")
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
    ## (for example with fewer genes than its subsampling needs). The switch
    ## is announced and recorded in the result's metadata.
    transform <- if (is.null(tx_mat)) "log_cpm" else "vst"
    if (is.null(tx_mat)) {
        if (vst) {
            message("[exposomeSC] DESeq2 VST ",
                    if (requireNamespace("DESeq2", quietly = TRUE)) {
                        "failed"
                    } else {
                        "is unavailable (DESeq2 is not installed)"
                    },
                    "; using log-CPM values for '", celltype, "'.")
        }
        tx_mat <- t(.log_cpm(pb$pb_mat))
    }
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
    ## A donor with a missing exposure or covariate cannot be adjusted; it is
    ## left out, as run_sc_exwas() leaves it out of its design
    adj_vars <- c(exposure, covariates)
    if (length(adj_vars)) {
        complete <- stats::complete.cases(
            exposureData(x)[shared_donors, adj_vars, drop = FALSE])
        if (!all(complete)) {
            missing_what <- if (is.null(exposure)) {
                "covariate"
            } else if (is.null(covariates)) {
                "exposure"
            } else {
                "exposure or covariate"
            }
            message("[exposomeSC] Leaving ", sum(!complete), " donor(s) ",
                    "with a missing ", missing_what, " value out of the ",
                    "network for '", celltype, "'.")
            shared_donors <- shared_donors[complete]
        }
        if (length(shared_donors) < 5L)
            stop("Only ", length(shared_donors), " of the donors shared ",
                "between pseudobulk and metabolites have values of the ",
                "adjustment variables. Need at least 5.", call. = FALSE)
    }

    tx_mat <- tx_mat[shared_donors, , drop = FALSE]
    met_mat <- metabolites[shared_donors, , drop = FALSE]
    if (ncol(met_mat) == 0L)
        stop("'metabolites' has no columns.", call. = FALSE)

    ## Unnamed columns are named before any is left out, so that the names
    ## in constant_features and in the node table agree
    colnames(tx_mat) <- .column_labels(tx_mat, "gene_")
    colnames(met_mat) <- .column_labels(met_mat, "metab_")

    ## --- Missing metabolite values ---
    ## The graphical lasso needs a finite value of every feature for every
    ## donor; the univariate screen of run_exposure_network() would take a
    ## metabolite with a missing value for a constant one and skip it, and
    ## stability selection would rank it only in the half-samples that leave
    ## out the donors whose value is missing
    met_values <- as.matrix(met_mat)
    incomplete <- colSums(is.na(met_values) | is.infinite(met_values)) > 0L
    if (any(incomplete)) {
        stop("'metabolites' has missing or infinite values for the ",
            "analysed donors in: ",
            paste(colnames(met_mat)[incomplete], collapse = ", "),
            ". Impute them or leave those metabolites out.", call. = FALSE)
    }

    ## --- Constant features ---
    ## A feature without variance across the donors (for example a gene
    ## with no counts in this cell type) has no partial correlation with
    ## anything and would stop the graphical lasso. The check comes before
    ## the adjustment below: regressing a column with the same non-zero
    ## value in every donor on the exposure or covariates leaves rounding
    ## noise, which scaling would turn into a feature of unit variance.
    ## Constant features are left out of the network and named in the
    ## result's metadata; when they are kept (for the selection stage of
    ## run_exposure_network()), they stay constant.
    tx_constant <- .constant_columns(tx_mat)
    met_constant <- .constant_columns(met_mat)
    constant_features <- character()
    if (drop_constant) {
        if (any(tx_constant) || any(met_constant)) {
            constant_features <- c(colnames(tx_mat)[tx_constant],
                colnames(met_mat)[met_constant])
            message("[exposomeSC] Leaving ", sum(tx_constant),
                    " constant transcript(s) and ", sum(met_constant),
                    " constant metabolite(s) out of the network for '",
                    celltype, "'.")
            tx_mat <- tx_mat[, !tx_constant, drop = FALSE]
            met_mat <- met_mat[, !met_constant, drop = FALSE]
            tx_constant <- logical(ncol(tx_mat))
            met_constant <- logical(ncol(met_mat))
        }
        if (ncol(tx_mat) == 0L || ncol(met_mat) == 0L) {
            stop("Every transcript or every metabolite is constant across ",
                "the donors of cell type '", celltype, "'.", call. = FALSE)
        }
    }

    ## --- Residualization (optional) ---
    if (!is.null(exposure) || !is.null(covariates)) {
        exp_data <- exposureData(x)
        adj_vars <- c(exposure, covariates)

        if (length(adj_vars) > 0) {
            adj_mat <- exp_data[shared_donors, adj_vars,
                                drop = FALSE]
            tx_mat <- .residualize_matrix(tx_mat, adj_mat)
            met_mat <- .residualize_matrix(met_mat, adj_mat)
            ## the residuals of a kept constant feature can be rounding
            ## noise; set them to 0 so that it stays constant
            tx_mat[, tx_constant] <- 0
            met_mat[, met_constant] <- 0
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
        n_metabolites  = length(met_names),
        constant_features = constant_features
    )
}

# Columns of a numeric matrix that take one value in every row; columns
# with missing values are not flagged (.assemble_crossomic() refuses
# missing metabolite values before it calls this)
# @keywords internal
# @noRd
.constant_columns <- function(m) {
    variances <- apply(m, 2L, stats::var)
    is.finite(variances) & variances == 0
}

# Column names, or a prefix and the column index where names are missing
# @keywords internal
# @noRd
.column_labels <- function(m, prefix) {
    labels <- colnames(m)
    if (is.null(labels)) {
        labels <- paste0(prefix, seq_len(ncol(m)))
    }
    labels
}

# Stop when an exposure or covariate name is not a column of exposureData(x)
# @keywords internal
# @noRd
.check_adjustment_names <- function(x, adjustment) {
    unknown <- setdiff(adjustment, colnames(exposureData(x)))
    if (length(unknown)) {
        stop("Exposure/covariate column(s) not found in exposureData: ",
            paste(unknown, collapse = ", "), call. = FALSE)
    }
    invisible(TRUE)
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
# Graphical lasso with separate within- and between-block
# penalties
# -------------------------------------------------------

#' Graphical lasso with block-specific penalties
#'
#' Fits the graphical lasso of Friedman, Hastie and Tibshirani
#' (2008) with \code{glasso::glasso()}. The penalty matrix holds
#' \code{lambda_w} for pairs of features in the same block and
#' \code{lambda_b} for pairs in different blocks; the diagonal is
#' not penalised. The covariance matrix is first scaled to a
#' correlation matrix, so the penalties do not depend on the units
#' of the features, and the returned precision matrix is on that
#' scale. It is symmetrised; an edge is present where the
#' symmetrised entry is non-zero.
#'
#' @param S sample covariance matrix (p x p)
#' @param blocks vector of block assignments, one per feature
#' @param lambda_w numeric; within-block penalty
#' @param lambda_b numeric; between-block penalty
#' @return list with wi (precision), adj (adjacency) and rho
#'   (penalty matrix)
#'
#' @keywords internal
#' @noRd
.block_glasso <- function(S, blocks, lambda_w = 0.3,
                           lambda_b = 0.5) {
    if (!requireNamespace("glasso", quietly = TRUE)) {
        stop("Package 'glasso' is required for the block graphical ",
             "lasso; install it with BiocManager::install(\"glasso\").",
             call. = FALSE)
    }
    for (penalty in list(lambda_w, lambda_b)) {
        if (!is.numeric(penalty) || length(penalty) != 1L ||
            !is.finite(penalty) || penalty < 0) {
            stop("lambda_w and lambda_b must each be a single ",
                 "non-negative number.", call. = FALSE)
        }
    }
    p <- nrow(S)
    if (length(blocks) != p) {
        stop("blocks must have one entry per feature.", call. = FALSE)
    }
    variances <- diag(S)
    if (any(!is.finite(variances)) || any(variances <= 0)) {
        stop("Every feature needs a positive, finite variance before ",
             "the graphical lasso; remove constant features first.",
             call. = FALSE)
    }

    rho_mat <- matrix(lambda_b, p, p)
    rho_mat[outer(blocks, blocks, "==")] <- lambda_w
    diag(rho_mat) <- 0

    fit <- glasso::glasso(stats::cov2cor(S), rho = rho_mat,
                          penalize.diagonal = FALSE)
    wi <- (fit$wi + t(fit$wi)) / 2

    adj <- .precision_to_adjacency(wi)
    list(wi = wi, adj = adj, rho = rho_mat)
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
