## Create demo SCEE data for vignette and figures
## Uses realistic parameters from published sc-eQTL studies
library(SingleCellExperiment)
library(S4Vectors)

set.seed(2024)

## 20 donors, 3 cell types, 100 cells/donor, 200 genes
n_donors <- 20
n_cells_per <- 100
n_genes <- 200
donors <- paste0("D", sprintf("%02d", seq_len(n_donors)))
cell_types <- c("Monocyte", "NK", "T_CD4")

## Generate count matrix
n_total <- n_donors * n_cells_per
counts <- matrix(rpois(n_genes * n_total, lambda = 8),
                  nrow = n_genes)
rownames(counts) <- paste0("Gene", seq_len(n_genes))
colnames(counts) <- paste0("cell_", seq_len(n_total))

## Add exposure-responsive signal to first 10 genes in Monocytes
## This simulates PM2.5 upregulating inflammatory genes
donor_ids <- rep(donors, each = n_cells_per)
ct_assign <- sample(cell_types, n_total, replace = TRUE,
                     prob = c(0.3, 0.2, 0.5))

## Exposure data: 5 pollutants + 2 covariates
exp_mat <- matrix(nrow = n_donors, ncol = 7,
    dimnames = list(donors,
        c("PM2.5", "NO2", "Pb", "Cd", "BPA",
          "age", "sex")))
exp_mat[, "PM2.5"] <- rnorm(n_donors, 25, 10)
exp_mat[, "NO2"] <- rnorm(n_donors, 30, 15)
exp_mat[, "Pb"] <- rlnorm(n_donors, 1, 0.5)
exp_mat[, "Cd"] <- rlnorm(n_donors, -1, 0.8)
exp_mat[, "BPA"] <- rlnorm(n_donors, 0, 0.6)
exp_mat[, "age"] <- sample(25:75, n_donors, replace = TRUE)
exp_mat[, "sex"] <- sample(0:1, n_donors, replace = TRUE)

## Add PM2.5 signal: higher PM2.5 -> higher counts for
## Gene1-Gene10 in Monocytes only
for (i in seq_len(n_total)) {
    d <- donor_ids[i]
    pm <- exp_mat[d, "PM2.5"]
    if (ct_assign[i] == "Monocyte") {
        boost <- pmax(0, round(pm / 10))
        counts[1:10, i] <- counts[1:10, i] + boost
    }
}

sce <- SingleCellExperiment(
    assays = list(counts = counts),
    colData = DataFrame(
        cell_id = colnames(counts),
        donor_id = donor_ids,
        cell_type = ct_assign))

outdir <- "inst/extdata"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
saveRDS(sce, file.path(outdir, "demo_sce.rds"))
write.csv(exp_mat, file.path(outdir, "demo_exposures.csv"))

cat("Demo data created:\n")
cat("  SCE:", ncol(sce), "cells,", nrow(sce), "genes,",
    n_donors, "donors\n")
cat("  Exposures:", ncol(exp_mat), "variables\n")
cat("  Signal: Gene1-10 upregulated by PM2.5 in Monocytes\n")
