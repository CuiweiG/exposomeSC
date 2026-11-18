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
#' @param transform character; "vst" (DESeq2) or "log_cpm", the transform
#'   of the pseudobulk counts
#' @param drop_constant logical; leave features that are constant across
#'   the donors out of the matrix (TRUE), or keep them as constant columns
#'   (FALSE)
#' @param genes character or NULL; row names of x, the genes to keep after
#'   the transform, which uses all genes
#' @param keep_unadjusted logical; also return the values before the
#'   exposure and the covariates are regressed out
#'
#' @return list with:
#'   - data_matrix: donors x features numeric matrix
#'   - node_info: DataFrame mapping features to omic layers
#'   - donors: character vector of donor IDs used
#'   - n_transcripts: integer
#'   - n_metabolites: integer
#'   - constant_features: character; features left out as constant
#'   - unadjusted: the values of the features before the adjustment and the
#'     scaling, or NULL
#'
#' @keywords internal
#' @noRd
.assemble_crossomic <- function(x, metabolites, celltype,
                                 celltype_col = "cell_type",
                                 exposure = NULL,
                                 covariates = NULL,
                                 min_cells = 10L,
                                 top_var_genes = NULL,
                                 transform = "vst", drop_constant = TRUE,
                                 genes = NULL, keep_unadjusted = FALSE) {
    ## --- Adjustment variables, checked before any computation ---
    .check_adjustment_names(x, c(exposure, covariates))
    .check_network_genes(genes, rownames(x))

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

    ## --- VST or log-CPM transform; tx_mat: donors x genes ---
    tx_mat <- .network_transform(pb$pb_mat, pb$valid_donors, transform,
        celltype)

    ## --- Optional: the genes asked for, in the order given ---
    ## The transform above has used every gene, so the values of a gene do
    ## not depend on which genes are kept
    if (!is.null(genes)) {
        tx_mat <- tx_mat[, as.character(genes), drop = FALSE]
    }

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
    ## a data frame is taken as the matrix of its columns, so that every
    ## later step, the adjustment by lm.fit() included, gets a matrix
    met_mat <- as.matrix(metabolites[shared_donors, , drop = FALSE])
    if (ncol(met_mat) == 0L)
        stop("'metabolites' has no columns.", call. = FALSE)
    if (!is.numeric(met_mat))
        stop("'metabolites' must hold numbers only.", call. = FALSE)

    ## Unnamed columns are named before any is left out, so that the names
    ## in constant_features and in the node table agree
    colnames(tx_mat) <- .column_labels(tx_mat, "gene_")
    colnames(met_mat) <- .column_labels(met_mat, "metab_")

    ## Two nodes of one name could not be told apart in any table of the
    ## result: two rows of x, two columns of the metabolites, or one of each
    node_names <- c(colnames(tx_mat), colnames(met_mat))
    repeated <- unique(node_names[duplicated(node_names)])
    if (length(repeated)) {
        stop("Two nodes of a network cannot have the same name: ",
            .some_names(repeated), ". Give the rows of 'x' and the columns ",
            "of 'metabolites' distinct names, or name the genes to use with ",
            "'genes'.", call. = FALSE)
    }

    ## --- Missing metabolite values ---
    ## The graphical lasso needs a finite value of every feature for every
    ## donor; the univariate screen of run_exposure_network() would take a
    ## metabolite with a missing value for a constant one and skip it, and
    ## stability selection would rank it only in the half-samples that leave
    ## out the donors whose value is missing
    incomplete <- colSums(is.na(met_mat) | is.infinite(met_mat)) > 0L
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
            counted <- function(constant, values, what) {
                paste0(sum(constant), " constant ", what, "(s)",
                    if (any(constant)) paste0(" (",
                        .some_names(colnames(values)[constant]), ")"))
            }
            message("[exposomeSC] Leaving ",
                    counted(tx_constant, tx_mat, "transcript"), " and ",
                    counted(met_constant, met_mat, "metabolite"),
                    " out of the network for '", celltype, "'.")
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

    ## The selection stage of run_exposure_network() ranks the values as they
    ## are before the adjustment
    unadjusted <- if (keep_unadjusted) cbind(tx_mat, met_mat) else NULL

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
        constant_features = constant_features,
        unadjusted     = unadjusted
    )
}

# 'genes' of the network functions: NULL, or distinct row names of x, each
# of which names one row.
.check_network_genes <- function(genes, feature_ids) {
    if (is.null(genes)) return(invisible(NULL))
    if (!(is.character(genes) || is.factor(genes)) || !length(genes) ||
        anyNA(genes) || anyDuplicated(as.character(genes))) {
        stop("'genes' must be NULL or a vector of distinct row names of 'x'.",
            call. = FALSE)
    }
    genes <- as.character(genes)
    unknown <- setdiff(genes, feature_ids)
    if (length(unknown)) {
        stop("'genes' not found in rownames(x): ", .some_names(unknown),
            call. = FALSE)
    }
    repeated <- intersect(genes, feature_ids[duplicated(feature_ids)])
    if (length(repeated)) {
        stop("'genes' that are the name of more than one row of 'x': ",
            .some_names(repeated), call. = FALSE)
    }
    invisible(NULL)
}

# The first names of a vector and the number of the others, for a message
.some_names <- function(names, n = 5L) {
    paste0(paste(utils::head(names, n), collapse = ", "),
        if (length(names) > n) paste0(" and ", length(names) - n, " more"))
}

# The transcript values of the network functions, donors x genes, from the
# genes x donors pseudobulk counts: DESeq2's variance-stabilising
# transformation, blind to the design ("vst"), or log2 counts per million
# ("log_cpm"). The transform is the one asked for or the call stops: the
# values are not switched when DESeq2 is missing or vst() fails, which it
# does for fewer than 1000 genes with a mean normalised count above 5.
.network_transform <- function(pb_mat, donors, transform, celltype) {
    if (transform == "log_cpm") {
        return(t(.log_cpm(pb_mat)))
    }
    if (!.is_installed("DESeq2")) {
        stop("Package 'DESeq2' is required for transform = \"vst\". ",
            "Install it with BiocManager::install('DESeq2'), or use ",
            "transform = \"log_cpm\".", call. = FALSE)
    }
    tryCatch({
        .relay_vst_notes({
            dds <- DESeq2::DESeqDataSetFromMatrix(
                countData = pb_mat,
                colData = data.frame(donor = donors),
                design = ~ 1)
            vsd <- DESeq2::vst(dds, blind = TRUE)
            t(SummarizedExperiment::assay(vsd))
        }, celltype)
    }, error = function(e) {
        reason <- conditionMessage(e)
        ## DESeq2 goes on to recommend a function that is not offered here
        if (grepl("less than 'nsub' rows", reason, fixed = TRUE)) {
            reason <- paste("fewer than 1000 genes have a mean normalised",
                "count above 5")
        }
        stop("DESeq2::vst() failed for '", celltype, "': ", reason,
            ". Use transform = \"log_cpm\" for counts that vst() cannot ",
            "transform.", call. = FALSE)
    })
}

# Evaluates 'expr' without the routine messages of DESeq2, except one: when
# the parametric fit of the dispersion trend fails, DESeq2 uses a local fit
# in its place and says so in a note. The values are then those of another
# transformation, so the user is told.
.relay_vst_notes <- function(expr, celltype) {
    withCallingHandlers(expr, message = function(m) {
        if (grepl("fitType", conditionMessage(m), fixed = TRUE)) {
            message("[exposomeSC] DESeq2::vst() for '", celltype, "': the ",
                "parametric fit of the dispersion trend failed, and DESeq2 ",
                "used a local fit in its place.")
        }
        tryInvokeRestart("muffleMessage")
    })
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
#' For each column in \code{Y}, fits the linear model of the column on an
#' intercept and the columns of \code{adj} and returns the residuals. The
#' model is given to \code{lm.fit()} as a design matrix: no formula is
#' built, so the names of the adjustment variables are never evaluated and
#' can be any names.
#'
#' @param Y numeric matrix (n x p)
#' @param adj numeric matrix (n x q) of adjustment variables
#' @return matrix of residuals (n x p)
#'
#' @keywords internal
#' @noRd
.residualize_matrix <- function(Y, adj) {
    design <- cbind(1, as.matrix(adj))
    result <- stats::lm.fit(design, Y)$residuals
    ## lm.fit() returns a vector for a response of one column
    dim(result) <- dim(Y)
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
# Edges of an adjacency matrix
# -------------------------------------------------------

#' Positions of the edges of an adjacency matrix
#'
#' The edges are the non-zero entries above the diagonal. They are returned
#' by row and, within a row, by column, the order in which a loop over
#' i < j meets them; \code{which(arr.ind = TRUE)} alone gives them by column.
#'
#' @param adjacency square matrix; entries on and below the diagonal are
#'   not read
#' @return integer matrix with one row per edge and the columns i and j
#'
#' @keywords internal
#' @noRd
.upper_edges <- function(adjacency) {
    is_edge <- upper.tri(adjacency) & adjacency != 0
    if (anyNA(is_edge)) {
        stop("The adjacency matrix has a missing value above its diagonal.")
    }
    ij <- which(is_edge, arr.ind = TRUE, useNames = FALSE)
    ij[order(ij[, 1L], ij[, 2L]), , drop = FALSE]
}

# -------------------------------------------------------
# Count cross-omic edges
# -------------------------------------------------------

#' Count edges between different omic layers
#'
#' @param adjacency binary matrix
#' @param node_info DataFrame with omic_layer column
#' @return list with total, within_transcript,
#'   within_metabolite, cross_omic counts
#'
#' @keywords internal
#' @noRd
.count_edge_types <- function(adjacency, node_info) {
    layers <- node_info$omic_layer
    ij <- .upper_edges(adjacency)
    layer_i <- layers[ij[, 1L]]
    same_layer <- layer_i == layers[ij[, 2L]]
    if (is.null(layers) || anyNA(same_layer)) {
        stop("A node with an edge has no omic layer.")
    }

    list(
        total = nrow(ij),
        within_transcript = sum(same_layer & layer_i == "transcript"),
        within_metabolite = sum(same_layer & layer_i != "transcript"),
        cross_omic = sum(!same_layer)
    )
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
#'   layer_i, layer_j, cross_omic columns, one row per edge in the order
#'   of .upper_edges()
#'
#' @keywords internal
#' @noRd
.extract_edges <- function(adjacency, node_info,
                            precision = NULL,
                            stability = NULL) {
    ij <- .upper_edges(adjacency)

    if (nrow(ij) == 0) {
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

    i <- ij[, 1L]
    j <- ij[, 2L]
    unknown <- rep(NA_real_, nrow(ij))
    ## row.names = NULL: data.frame() would otherwise take the row names
    ## from names that a column of the node table may carry
    S4Vectors::DataFrame(data.frame(
        node_i     = node_info$feature[i],
        node_j     = node_info$feature[j],
        weight     = if (!is.null(precision)) precision[ij] else unknown,
        stability  = if (!is.null(stability)) stability[ij] else unknown,
        layer_i    = node_info$omic_layer[i],
        layer_j    = node_info$omic_layer[j],
        cross_omic = node_info$omic_layer[i] != node_info$omic_layer[j],
        row.names = NULL,
        stringsAsFactors = FALSE
    ))
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
# Rank correlation with the exposure, given the covariates
# -------------------------------------------------------

# Residuals of the columns of 'ranks' on the columns of 'adjustment_ranks'.
# A column that is constant, or that the adjustment reproduces, is returned
# as zeros, not as the rounding noise of its residuals.
.rank_residuals <- function(ranks, adjustment_ranks) {
    total <- colSums(sweep(ranks, 2L, colMeans(ranks))^2)
    residuals <- .residualize_matrix(ranks, adjustment_ranks)
    reproduced <- total == 0 | colSums(residuals^2) <= 1e-10 * total
    residuals[, reproduced] <- 0
    residuals
}

# Spearman correlation of each column of 'values' with 'response', and its
# two-sided p-value. With 'adjustment', a matrix of covariates, it is the
# partial rank correlation: the ranks of the covariates are regressed out of
# the ranks of the values and of the response, and the residuals are
# correlated. Regressing the covariates out of the values themselves and
# ranking the residuals would not do: values that are tied, as the zeros of
# a gene with few counts, would be ranked by the covariates. The test is the
# t test with n - 2 - k degrees of freedom for k covariates, which without
# covariates is that of cor.test(method = "spearman", exact = FALSE). A
# column or a response without variation after the adjustment gives NA.
.partial_rank_correlation <- function(values, response, adjustment = NULL) {
    k <- if (is.null(adjustment)) 0L else ncol(adjustment)
    feature_ranks <- apply(values, 2L, rank)
    response_ranks <- rank(response)
    if (k > 0L) {
        adjustment_ranks <- apply(adjustment, 2L, rank)
        ## covariates that repeat each other, or that are constant here,
        ## take no degree of freedom of their own
        k <- qr(cbind(1, adjustment_ranks))$rank - 1L
        feature_ranks <- .rank_residuals(feature_ranks, adjustment_ranks)
        response_ranks <- .rank_residuals(cbind(response_ranks),
            adjustment_ranks)[, 1L]
    }
    varying <- function(v) isTRUE(stats::sd(v) > 0)
    usable <- apply(feature_ranks, 2L, varying)
    rho <- rep(NA_real_, ncol(values))
    if (varying(response_ranks) && any(usable)) {
        rho[usable] <- stats::cor(feature_ranks[, usable, drop = FALSE],
            response_ranks)[, 1L]
    }
    rho <- pmax(pmin(rho, 1), -1)
    df <- nrow(values) - 2L - k
    statistic <- rho * sqrt(df / (1 - rho^2))
    list(rho = rho, p_value = 2 * stats::pt(-abs(statistic), df))
}

# -------------------------------------------------------
# Base selector for stability selection: 0/1 indicator of the q
# features with the largest absolute Spearman correlation with y, given the
# covariates in 'adjustment'. Constant features are never selected.
# -------------------------------------------------------
#' @keywords internal
.top_marginal_features <- function(x, y, q, adjustment = NULL) {
    chosen <- numeric(ncol(x))
    if (!isTRUE(stats::sd(y) > 0)) return(chosen)
    varies <- which(apply(x, 2L, function(v) isTRUE(stats::sd(v) > 0)))
    if (length(varies) == 0L) return(chosen)
    rho <- abs(.partial_rank_correlation(x[, varies, drop = FALSE], y,
        adjustment)$rho)
    varies <- varies[!is.na(rho)]
    rho <- rho[!is.na(rho)]
    ranked <- varies[order(-rho, varies)]
    chosen[ranked[seq_len(min(q, length(varies)))]] <- 1
    chosen
}
