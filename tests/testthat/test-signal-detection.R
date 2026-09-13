## Statistical validation: can the pipeline detect a KNOWN signal?
## This is the most important test for a statistical package.

library(SingleCellExperiment)
library(S4Vectors)

test_that("run_sc_exwas detects planted exposure signal", {
    set.seed(2024)
    n_donors <- 20
    n_cells_per <- 50
    n_genes <- 100
    n_total <- n_donors * n_cells_per

    donors <- paste0("D", sprintf("%02d", seq_len(n_donors)))
    donor_ids <- rep(donors, each = n_cells_per)
    cell_types <- rep(c("TypeA", "TypeB"),
                       length.out = n_total)

    ## Baseline counts
    counts <- matrix(rpois(n_genes * n_total, lambda = 10),
                      nrow = n_genes)
    rownames(counts) <- paste0("Gene", seq_len(n_genes))
    colnames(counts) <- paste0("cell_", seq_len(n_total))

    ## Exposure data
    exp_vals <- rnorm(n_donors, mean = 0, sd = 1)
    exp_mat <- matrix(exp_vals, nrow = n_donors, ncol = 1,
        dimnames = list(donors, "E1"))

    ## Plant signal: Gene1-Gene5 upregulated by E1 in TypeA
    ## ONLY. Add 3*exposure to counts for those genes.
    for (i in seq_len(n_total)) {
        d_idx <- match(donor_ids[i], donors)
        if (cell_types[i] == "TypeA") {
            boost <- pmax(0, round(3 * exp_vals[d_idx] + 3))
            counts[1:5, i] <- counts[1:5, i] + boost
        }
    }

    sce <- SingleCellExperiment(
        assays = list(counts = counts),
        colData = DataFrame(
            cell_id = colnames(counts),
            donor_id = donor_ids,
            cell_type = cell_types))
    scee <- build_scee(sce, exp_mat, sample_col = "donor_id")

    result <- run_sc_exwas(scee,
        exposure = "E1",
        celltype_col = "cell_type",
        sample_col = "donor_id",
        min_cells = 5L,
        min_donors = 5L)

    ## Check that at least some signal genes in TypeA have
    ## small p-values (top results should include Gene1-5)
    typeA <- result[result$celltype == "TypeA", ]
    typeA_sorted <- typeA[order(typeA$pvalue), ]
    top5_genes <- head(typeA_sorted$gene, 5)

    signal_genes <- paste0("Gene", 1:5)
    overlap <- length(intersect(top5_genes, signal_genes))

    ## At least 3/5 signal genes should be in top 5
    expect_true(overlap >= 3,
        info = paste("Expected >= 3 signal genes in top 5,",
                     "found", overlap, ":",
                     paste(top5_genes, collapse = ", ")))

    ## TypeB should NOT have the same signal
    typeB <- result[result$celltype == "TypeB", ]
    typeB_sig <- typeB[!is.na(typeB$pvalue) &
                       typeB$pvalue < 0.01 &
                       typeB$gene %in% signal_genes, ]
    ## Allow at most 1 false positive in TypeB
    expect_true(nrow(typeB_sig) <= 1,
        info = paste("TypeB false positives:",
                     nrow(typeB_sig)))
})
