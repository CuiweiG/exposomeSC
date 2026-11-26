# Donor offset-pseudobulk negative-binomial sc-ExWAS

A donor-level negative-binomial alternative to
[`run_sc_exwas`](run_sc_exwas.md). Per cell type and per gene,
aggregates cells to a donor pseudobulk count and fits a
negative-binomial GLM with a library-size offset \\\log s_i =
\log(\text{donor total counts})\\. Lee and Han (2024) report that
pseudobulk models with proper offsets have the same statistical
properties as generalised linear mixed models in single-cell
case-control studies, without the pseudoreplication that inflates
cell-level tests. P-values are Wald tests with an estimated dispersion
and can be anticonservative when there are few donors. The function has
no threshold on the number of cells: every donor with at least one cell
of a cell type enters the model of that cell type, where
[`run_sc_exwas`](run_sc_exwas.md) leaves out a donor with fewer than
`min_cells` cells.

## Usage

``` r
run_sc_exwas_pb_offset(
  x,
  exposure,
  celltype_col = "cell_type",
  covariates = NULL,
  target_genes = NULL,
  BPPARAM = BiocParallel::SerialParam()
)
```

## Arguments

- x:

  A
  [`SingleCellExposomeExperiment`](SingleCellExposomeExperiment-class.md).
  Each cell's donor is read from the `colData` column named by
  [`sampleColumn`](SCEE-accessors.md)`(x)`.

- exposure:

  Character; a donor-level column of `exposureData(x)`.

- celltype_col:

  Character. Column in `colData` containing cell type annotations.
  Default `"cell_type"`.

- covariates:

  Character vector; donor covariate columns in `exposureData` (numeric),
  added as GLM main effects. The names `y`, `A` and `logoff` are used by
  the model and cannot be covariate names, and each name must be a
  syntactic R name (see
  [`make.names`](https://rdrr.io/r/base/make.names.html)) other than
  `.`, `...` and names such as `..1`, which R reads specially in a
  formula. The exposure itself cannot be a covariate.

- target_genes:

  Character; genes to test (default: all rows of `x`).

- BPPARAM:

  A `BiocParallel` param used to parallelise over genes. Defaults to
  [`BiocParallel::SerialParam()`](https://rdrr.io/pkg/BiocParallel/man/SerialParam-class.html),
  so no worker process is started unless one is asked for; pass, for
  example, `BiocParallel::SnowParam(workers = 2L)` to fit genes in
  parallel. Pass `NULL` to use a plain `lapply`.

## Value

A `DataFrame`, one row per gene x cell type: `gene`, `celltype`,
`log2FC` (exposure effect on the \\\log_2\\ rate scale), `se` (on the
same scale), `statistic` (donor-level Wald z from the NB GLM), `pvalue`,
and `padj_global` (BH across all tested gene-by-cell-type strata, as in
[`run_sc_exwas`](run_sc_exwas.md); there is no within-cell-type `padj`).
Rows ordered by `padj_global`, then `pvalue`. A gene whose GLM cannot be
fitted gets `NA` estimates; the function stops when no cell type has at
least three donors with cells.

## Note

Report empirical calibration from a donor-permutation null, with the
exposure permuted across donors, before quoting any type-I rate.

## References

Lee H & Han B (2024). Pseudobulk with proper offsets has the same
statistical properties as generalized linear mixed models in single-cell
case-control studies. *Bioinformatics* 40(8):btae498.
[doi:10.1093/bioinformatics/btae498](https://doi.org/10.1093/bioinformatics/btae498)

Squair JW et al. (2021). Confronting false discoveries in single-cell
differential expression. *Nat Commun* 12:5692.

## See also

[`run_sc_exwas_glmm`](run_sc_exwas_glmm.md) for the cell-level LMM
reference.

## Examples

``` r
donor_ids <- paste0("D", seq_len(10L))
exposure <- rep(c(0, 1), each = 5L)
cell_donor <- rep(donor_ids, each = 10L)
cell_type <- rep("TypeA", length(cell_donor))
cell_ids <- paste0("cell", seq_along(cell_donor))
pb_targets <- rbind(
    G1 = c(30, 140, 50, 180, 40, 110, 260, 70, 220, 130),
    G2 = c(200, 40, 240, 60, 300, 50, 270, 80, 330, 70),
    G3 = c(100, 120, 90, 150, 110, 130, 170, 100, 160, 140),
    G4 = c(250, 200, 300, 180, 270, 220, 320, 210, 290, 230)
)
counts <- matrix(0L, nrow = 4L, ncol = length(cell_ids))
for (donor in seq_along(donor_ids)) {
    counts[, cell_donor == donor_ids[donor]] <- pb_targets[, donor] / 10L
}
storage.mode(counts) <- "integer"
dimnames(counts) <- list(paste0("G", seq_len(4L)), cell_ids)
sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = counts),
    colData = S4Vectors::DataFrame(
        cell_id = cell_ids,
        donor_id = cell_donor,
        cell_type = cell_type
    )
)
donor_design <- matrix(exposure, ncol = 1L,
    dimnames = list(donor_ids, "exposure"))
scee <- build_scee(sce, donor_design, sample_col = "donor_id")
offset_result <- run_sc_exwas_pb_offset(
    scee,
    exposure = "exposure",
    celltype_col = "cell_type",
    target_genes = c("G1", "G2"),
    BPPARAM = BiocParallel::SerialParam()
)
head(offset_result)
#> DataFrame with 2 rows and 7 columns
#>          gene    celltype    log2FC        se statistic    pvalue padj_global
#>   <character> <character> <numeric> <numeric> <numeric> <numeric>   <numeric>
#> 1          G1       TypeA  0.468350  0.507553  0.922762  0.356131    0.426691
#> 2          G2       TypeA -0.390771  0.491618 -0.794867  0.426691    0.426691
```
