# Run sc-ExWAS across multiple exposures

Convenience wrapper that runs [`run_sc_exwas`](run_sc_exwas.md) for each
exposure and combines results with a final global FDR correction across
all exposures and cell types.

## Usage

``` r
run_multi_exwas(
  x,
  exposures = NULL,
  celltype_col = "cell_type",
  covariates = NULL,
  min_cells = 10L,
  min_donors = 5L,
  filter_genes = TRUE,
  ...
)
```

## Arguments

- x:

  A
  [`SingleCellExposomeExperiment`](SingleCellExposomeExperiment-class.md).
  Each cell's donor is read from the `colData` column named by
  [`sampleColumn`](SCEE-accessors.md)`(x)`.

- exposures:

  Character vector. Exposure names to test. Default: all exposures in
  `exposureData` except the `covariates`.

- celltype_col:

  Character. `colData` column with cell type labels. Default
  `"cell_type"`.

- covariates:

  Character vector (optional). Columns of `exposureData` adjusted for in
  every model; they are not tested as exposures, and listing one in
  `exposures` is an error.

- min_cells:

  Integer. Default 10.

- min_donors:

  Integer. Default 5.

- filter_genes:

  Logical. Default TRUE.

- ...:

  Further arguments of `run_sc_exwas`, such as `method` or `celltypes`.
  An argument that `run_sc_exwas` does not accept is an error.

## Value

A `DataFrame` with all results combined and `padj_all` for global FDR
across all tests. The `se` and `statistic` fields retain the selected
backend's semantics: edgeR QL has numeric-`NA` coefficient standard
errors and a signed square-root QL F statistic, whereas DESeq2 and
voom-dream return their coefficient standard errors with Wald and
moderated t statistics, respectively. The columns are those of
[`run_sc_exwas`](run_sc_exwas.md) followed by `padj_all`;
`metadata(result)$per_exposure` keeps the
[`metadata()`](https://rdrr.io/pkg/S4Vectors/man/Annotated-class.html)
(parameters and diagnostics) of each exposure's `run_sc_exwas` result.

## Details

If `run_sc_exwas` stops for one of the exposures, the call stops too,
with an error that names the exposure, instead of leaving that exposure
out of `padj_all`. A cell type that `run_sc_exwas` skips for an exposure
with a warning, for example for too few donors with complete data, adds
no tests to `padj_all`; the diagnostics in
`metadata(result)$per_exposure` list the cell types fitted for each
exposure.

## Examples

``` r
library(SingleCellExperiment); library(S4Vectors)
set.seed(1)
counts <- matrix(rpois(5000, 8), nrow = 50,
    dimnames = list(paste0("G", 1:50), paste0("c", 1:100)))
sce <- SingleCellExperiment(assays = list(counts = counts),
    colData = DataFrame(cell_id = paste0("c", 1:100),
        donor_id = rep(paste0("D", 1:5), each = 20),
        cell_type = rep(c("Mono", "NK"), 50)))
exp_mat <- matrix(rnorm(15), nrow = 5,
    dimnames = list(paste0("D", 1:5), c("E1", "E2", "E3")))
scee <- build_scee(sce, exp_mat, sample_col = "donor_id")
multi <- run_multi_exwas(scee,
    exposures = c("E1", "E2"),
    celltype_col = "cell_type",
    min_cells = 5L, min_donors = 3L)
#> Running sc-ExWAS for: E1
#> Running sc-ExWAS for: E2
head(multi)
#> DataFrame with 6 rows and 18 columns
#>          gene    celltype    log2FC        se statistic      pvalue
#>   <character> <character> <numeric> <numeric> <numeric>   <numeric>
#> 1         G30          NK -0.274903        NA  -4.36521 0.000385337
#> 2         G32        Mono  0.239709        NA   4.42734 0.000431488
#> 3         G22          NK -0.335454        NA  -4.07881 0.000724244
#> 4         G42          NK  0.173579        NA   2.62443 0.017338878
#> 5          G1          NK -0.192264        NA  -2.43200 0.025856685
#> 6         G10        Mono  0.224625        NA   2.42016 0.027904356
#>   pvalue_underflow_clamped      padj padj_global    logCPM    exposure
#>                  <logical> <numeric>   <numeric> <numeric> <character>
#> 1                    FALSE 0.0181061   0.0215744   14.4285          E1
#> 2                    FALSE 0.0215744   0.0215744   14.3766          E1
#> 3                    FALSE 0.0181061   0.0241415   14.2140          E1
#> 4                    FALSE 0.2889813   0.4334720   14.4099          E1
#> 5                    FALSE 0.3232086   0.4650726   14.2384          E1
#> 6                    FALSE 0.5278257   0.4650726   14.2152          E1
#>    n_donors n_unexposed n_exposed min_cells median_cells          method
#>   <integer>   <integer> <integer> <integer>    <numeric>     <character>
#> 1         5          NA        NA        10           10 edgeR_robust_QL
#> 2         5          NA        NA        10           10 edgeR_robust_QL
#> 3         5          NA        NA        10           10 edgeR_robust_QL
#> 4         5          NA        NA        10           10 edgeR_robust_QL
#> 5         5          NA        NA        10           10 edgeR_robust_QL
#> 6         5          NA        NA        10           10 edgeR_robust_QL
#>    padj_all
#>   <numeric>
#> 1 0.0431488
#> 2 0.0431488
#> 3 0.0482830
#> 4 0.4206895
#> 5 0.4318579
#> 6 0.4318579
```
