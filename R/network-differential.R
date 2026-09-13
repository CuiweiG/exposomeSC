# R/network-differential.R
# Exposure-Interaction Network
# Differential precision matrix: high vs low exposure

#' @importFrom stats cor quantile median
NULL

#' Exposure-Stratified Differential Network
#'
#' Estimates separate precision matrices for donors above
#' and below the exposure median, then tests for
#' differential edges using permutation.
#'
#' @param scee A \code{SingleCellExposomeExperiment}.
#' @param exposure Character. Exposure variable.
#' @param celltype Character. Cell type to analyse.
#' @param celltype_col Character. Column with cell type labels.
#' @param sample_col Character. Column with donor IDs.
#' @param genes Character vector. Genes to include.
#' @param n_perm Integer. Permutations for edge testing.
#'   Default 100.
#' @param lambda Numeric. Graphical lasso penalty.
#'   Default 0.3.
#' @param min_cells Integer. Min cells per donor. Default 10.
#'
#' @return A list with:
#'   \describe{
#'     \item{delta}{Differential precision matrix (high - low)}
#'     \item{edges}{data.frame of differential edges with
#'       permutation p-values}
#'     \item{adj_high, adj_low}{Adjacency matrices}
#'     \item{n_high, n_low}{Number of donors per group}
#'   }
#'
#' @details
#' \strong{Mathematical framework:}
#'
#' For donors with exposure above/below median:
#' \deqn{\hat{\Omega}_{\text{high}} = \arg\min_{\Omega \succ 0}
#'   \{-\log|\Omega| + \text{tr}(S_H \Omega)
#'   + \lambda||\Omega||_1\}}
#' \deqn{\hat{\Omega}_{\text{low}} = \text{same for low group}}
#'
#' Differential edge \eqn{(i,j)}: \eqn{\Delta_{ij} =
#'   |\Omega_H|_{ij} - |\Omega_L|_{ij}}
#'
#' Significance via label permutation: shuffle exposure labels
#' across donors, re-estimate both networks, compute null
#' distribution of \eqn{\Delta_{ij}}.
#'
#' Requires the \pkg{glasso} package.
#'
#' @examples
#' donor_ids <- paste0("D", seq_len(10L))
#' exposure <- rep(c(0, 1), each = 5L)
#' cell_donor <- rep(donor_ids, each = 10L)
#' cell_ids <- paste0("cell", seq_along(cell_donor))
#' counts <- outer(
#'     seq_len(5L),
#'     seq_along(cell_ids),
#'     function(gene, cell) 3L + ((11L * gene + 7L * cell) %% 13L)
#' )
#' exposed_cells <- exposure[match(cell_donor, donor_ids)] == 1
#' counts[1L, exposed_cells] <- counts[1L, exposed_cells] + 4L
#' storage.mode(counts) <- "integer"
#' dimnames(counts) <- list(paste0("G", seq_len(5L)), cell_ids)
#' sce <- SingleCellExperiment::SingleCellExperiment(
#'     assays = list(counts = counts),
#'     colData = S4Vectors::DataFrame(
#'         cell_id = cell_ids,
#'         donor_id = cell_donor,
#'         cell_type = "Monocyte"
#'     )
#' )
#' donor_design <- matrix(exposure, ncol = 1L,
#'     dimnames = list(donor_ids, "exposure"))
#' scee <- build_scee(sce, donor_design, sample_col = "donor_id")
#' differential <- run_differential_network(
#'     scee,
#'     exposure = "exposure",
#'     celltype = "Monocyte",
#'     genes = paste0("G", seq_len(5L)),
#'     n_perm = 2L,
#'     min_cells = 5L
#' )
#' head(differential$edges)
#' @export
run_differential_network <- function(scee, exposure, celltype,
                                      celltype_col = "cell_type",
                                      sample_col = "donor_id",
                                      genes,
                                      n_perm = 100L,
                                      lambda = 0.3,
                                      min_cells = 10L) {

    stopifnot(is(scee, "SingleCellExposomeExperiment"))

    exp_data <- exposureData(scee)
    exp_vec <- exp_data[, exposure]
    donors <- rownames(exp_data)

    ## Pseudobulk
    cd <- SummarizedExperiment::colData(scee)
    counts_mat <- SummarizedExperiment::assay(scee, "counts")
    samp <- as.character(cd[[sample_col]])
    ct <- as.character(cd[[celltype_col]])

    pb <- .pseudobulk_aggregate(counts_mat, samp, ct,
        celltype, min_cells = min_cells)
    if (is.null(pb)) stop("No valid donors")

    valid <- pb$valid_donors
    pb_mat <- pb$pb_mat
    genes <- intersect(genes, rownames(pb_mat))
    if (length(genes) < 5) stop("Need >= 5 genes")

    pb_sub <- pb_mat[genes, valid, drop = FALSE]
    lcpm <- .log_cpm(pb_sub)

    ## Split by exposure median
    exp_valid <- exp_vec[valid]
    med_exp <- median(exp_valid, na.rm = TRUE)
    ## Use >= median for high, < for low (handles discrete exposures)
    high_idx <- which(exp_valid >= med_exp)
    low_idx <- which(exp_valid < med_exp)
    ## If all in one group, try strict > / <=
    if (length(high_idx) == 0 || length(low_idx) == 0) {
        high_idx <- which(exp_valid > med_exp)
        low_idx <- which(exp_valid <= med_exp)
    }
    if (length(high_idx) < 5 || length(low_idx) < 5)
        stop("Not enough donors in each group (need >= 5). ",
             "Check exposure distribution.")

    n_high <- length(high_idx)
    n_low <- length(low_idx)
    message(sprintf("Differential network: %d high, %d low donors",
        n_high, n_low))

    if (!requireNamespace("glasso", quietly = TRUE))
        stop("Package 'glasso' is required for run_differential_network(). ",
             "Install it with BiocManager::install('glasso').")

    ## Estimate precision matrices via glasso
    fit_glasso <- function(data_mat) {
        S <- cor(t(data_mat))
        fit <- glasso::glasso(S, rho = lambda)
        adj <- (abs(fit$wi) > 1e-6) * 1
        diag(adj) <- 0
        list(omega = fit$wi, adj = adj)
    }

    net_high <- fit_glasso(lcpm[, high_idx, drop = FALSE])
    net_low <- fit_glasso(lcpm[, low_idx, drop = FALSE])

    ## Differential matrix
    delta <- abs(net_high$omega) - abs(net_low$omega)
    rownames(delta) <- colnames(delta) <- genes

    ## Extract edges with |delta| > 0
    edge_idx <- which(abs(delta) > 1e-6 & upper.tri(delta),
        arr.ind = TRUE)

    if (nrow(edge_idx) == 0) {
        message("No differential edges found")
        return(list(delta = delta,
            edges = data.frame(gene1 = character(),
                gene2 = character(), delta = numeric(),
                perm_p = numeric()),
            adj_high = net_high$adj,
            adj_low = net_low$adj,
            n_high = n_high, n_low = n_low))
    }

    ## Observed delta for each edge
    obs_delta <- vapply(seq_len(nrow(edge_idx)), function(k) {
        delta[edge_idx[k, 1], edge_idx[k, 2]]
    }, numeric(1))

    ## Permutation test
    message(sprintf("Running %d permutations...", n_perm))
    null_deltas <- matrix(0, nrow = n_perm, ncol = length(obs_delta))

    for (p in seq_len(n_perm)) {
        perm_idx <- sample(ncol(lcpm))
        perm_high <- perm_idx[seq_len(n_high)]
        perm_low <- perm_idx[(n_high + 1):ncol(lcpm)]
        nh <- fit_glasso(lcpm[, perm_high, drop = FALSE])
        nl <- fit_glasso(lcpm[, perm_low, drop = FALSE])
        pd <- abs(nh$omega) - abs(nl$omega)
        null_deltas[p, ] <- vapply(seq_len(nrow(edge_idx)),
            function(k) pd[edge_idx[k, 1], edge_idx[k, 2]],
            numeric(1))
    }

    ## P-values (two-sided)
    perm_p <- vapply(seq_along(obs_delta), function(k) {
        (sum(abs(null_deltas[, k]) >= abs(obs_delta[k])) + 1) /
            (n_perm + 1)
    }, numeric(1))

    edges <- data.frame(
        gene1 = genes[edge_idx[, 1]],
        gene2 = genes[edge_idx[, 2]],
        delta = obs_delta,
        perm_p = perm_p,
        stringsAsFactors = FALSE)
    edges <- edges[order(edges$perm_p), ]
    rownames(edges) <- NULL

    list(delta = delta, edges = edges,
        adj_high = net_high$adj, adj_low = net_low$adj,
        n_high = n_high, n_low = n_low)
}
