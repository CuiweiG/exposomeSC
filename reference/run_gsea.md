# Gene set enrichment on sc-ExWAS results

Runs fast gene set enrichment analysis (fgsea) on the ranked gene list
from [`run_sc_exwas`](run_sc_exwas.md), separately for each cell type
and exposure.

## Usage

``` r
run_gsea(
  exwas_result,
  gene_sets,
  rank_by = "statistic",
  min_size = 10L,
  max_size = 500L
)
```

## Arguments

- exwas_result:

  A `DataFrame` returned by [`run_sc_exwas`](run_sc_exwas.md) or
  [`run_multi_exwas`](run_multi_exwas.md).

- gene_sets:

  A named list of character vectors, one gene set each. The genes are
  named as in the `gene` column of `exwas_result`, which holds the row
  names of the object that was analysed: Ensembl identifiers when the
  rows are Ensembl identifiers, symbols when they are symbols.

- rank_by:

  Character scalar naming the numeric column used to rank genes. Default
  `"statistic"`, whose scale depends on the sc-ExWAS backend: signed
  square-root QL F for edgeR, Wald statistic for DESeq2, and moderated t
  statistic for voom-dream. `"log2FC"` ranks by effect size; any other
  numeric column of `exwas_result` may be named.

- min_size:

  Integer. Minimum gene set size. Default 10.

- max_size:

  Integer. Maximum gene set size. Default 500.

## Value

A `DataFrame` with the `fgsea` columns pathway, pval, padj (within cell
type and exposure), log2err, ES, NES, size and leadingEdge
(semicolon-separated), plus celltype, exposure and padj_global
(Benjamini-Hochberg across all cell types and exposures).

## Details

This function requires the `fgsea` package (Bioconductor). Gene sets can
come from any source, for example from MSigDB via `msigdbr::msigdbr()`,
as long as they use the identifiers of the result.

Genes are ranked by the backend-specific association statistic (default)
within each cell type, then fgsea is run per cell type. For edgeR QL
this is the signed square root of the QL F statistic and is not a Wald z
statistic. Results are combined with a final cross-cell-type FDR
adjustment.

Ranking by the signed statistic preserves direction and within-backend
evidence ordering. Statistic magnitudes are not calibrated for
comparison across different backends; use `rank_by = "log2FC"` when an
effect-size ranking is the scientific target.

`fgsea` estimates p-values by Monte Carlo sampling, so they depend on
the random number generator state. `run_gsea()` runs `fgsea` serially,
so calling [`set.seed()`](https://rdrr.io/r/base/Random.html) beforehand
makes them reproducible. When a gene appears more than once in a
stratum, its largest ranking value is kept.

## References

Korotkevich G et al. (2021). Fast gene set enrichment analysis.
*bioRxiv*. [doi:10.1101/060012](https://doi.org/10.1101/060012)

Subramanian A et al. (2005). Gene set enrichment analysis. *PNAS*
102:15545-15550.
[doi:10.1073/pnas.0506580102](https://doi.org/10.1073/pnas.0506580102)

## Examples

``` r
set.seed(1)
genes <- paste0("G", 1:200)
exwas <- S4Vectors::DataFrame(gene = genes, celltype = "Mono",
    exposure = "PM2.5",
    statistic = c(stats::rnorm(20, mean = 3), stats::rnorm(180)))
gene_sets <- list(shifted = genes[1:20], random = sample(genes, 25))
if (requireNamespace("fgsea", quietly = TRUE)) {
    run_gsea(exwas, gene_sets, min_size = 10L)
}
#> DataFrame with 2 rows and 11 columns
#>       pathway        pval        padj   log2err        ES       NES      size
#>   <character>   <numeric>   <numeric> <numeric> <numeric> <numeric> <integer>
#> 1      random 5.72581e-01 5.72581e-01  0.165657 -0.255891 -0.917942        25
#> 2     shifted 9.36476e-17 1.87295e-16  1.057464  0.965471  2.378020        20
#>              leadingEdge    celltype    exposure padj_global
#>              <character> <character> <character>   <numeric>
#> 1 G141;G168;G35;G97;G7..        Mono       PM2.5 5.72581e-01
#> 2 G4;G11;G15;G18;G19;G..        Mono       PM2.5 1.87295e-16
```
