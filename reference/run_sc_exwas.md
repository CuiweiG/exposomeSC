# Run cell-type-specific ExWAS

Performs pseudobulk differential expression analysis with a continuous
exposure variable as covariate. Aggregates counts within each
donor-celltype stratum to avoid pseudoreplication, then uses an
empirical-Bayes pseudobulk model for inference.

## Usage

``` r
run_sc_exwas(x, ...)

# S4 method for class 'SingleCellExposomeExperiment'
run_sc_exwas(
  x,
  exposure,
  celltype_col = "cell_type",
  celltypes = NULL,
  covariates = NULL,
  min_cells = 10L,
  min_donors = 5L,
  filter_genes = TRUE,
  method = c("edgeR", "DESeq2", "voom_dream"),
  min_group_donors = 5L,
  filter_min_count = 10L,
  filter_min_total_count = 15L,
  filter_method = c("fixed_support", "filterByExpr"),
  filter_min_cpm = 1,
  filter_min_donors = 5L,
  test_features = NULL,
  robust = TRUE,
  ...
)
```

## Arguments

- x:

  A
  [`SingleCellExposomeExperiment`](SingleCellExposomeExperiment-class.md).
  Each cell's donor is read from the `colData` column named by
  [`sampleColumn`](SCEE-accessors.md)`(x)`.

- ...:

  Must be empty. Any further argument is an error, so that a misspelled
  argument name cannot silently change the analysis.

- exposure:

  Character. Name of the exposure variable (column in `exposureData`).

- celltype_col:

  Character. Column in `colData` containing cell type annotations.
  Default `"cell_type"`.

- celltypes:

  Character vector (optional). Which cell types to test. Default: all.

- covariates:

  Character vector (optional). Additional covariates from
  `exposureData`. A covariate cannot be called `"exposure"`, the name of
  the model's exposure term, or be the exposure itself, and each name
  must be a syntactic R name (see
  [`make.names`](https://rdrr.io/r/base/make.names.html)) because it
  goes into the model formula; `.`, `...` and names such as `..1`, which
  R reads specially there, cannot be used.

- min_cells:

  Integer. Minimum cells per donor-celltype stratum. Default 10.

- min_donors:

  Integer. Minimum donors per cell type. Default 5, and at least 3, for
  every backend. Cell types with fewer valid donors are skipped with a
  warning.

- filter_genes:

  Logical. If `TRUE` (default), filter low-expression genes before model
  fitting. The edgeR backend defaults to exposure-label-independent
  donor-support filtering; the DESeq2 and voom-dream backends retain
  genes detected in at least 20% of donors. If `FALSE`, every backend
  leaves out, in each cell type, only the genes without any count in the
  donors analysed for it.

- method:

  Character. Backend for differential expression. `"edgeR"` (default)
  uses the edgeR v4 quasi-likelihood pipeline (see Details). `"DESeq2"`
  uses DESeq2 on aggregated counts. `"voom_dream"` uses
  variancePartition's `voomWithDreamWeights()` and `dream()` with
  precision weights from the mean-variance trend (Hoffman and Roussos
  2021). Pseudobulk data have one observation per donor and the model
  has no random effects, so `dream()` fits it with limma and the
  moderated t statistic uses the residual degrees of freedom. This
  backend does not use the dreamlet package.

- min_group_donors:

  Integer. For a binary exposure, minimum complete donors required in
  each group for the edgeR backend, at least 2. Default 5.

- filter_min_count:

  Numeric. Minimum count passed to
  [`edgeR::filterByExpr`](https://rdrr.io/pkg/edgeR/man/filterByExpr.html)
  when `filter_method = "filterByExpr"`.

- filter_min_total_count:

  Numeric. Minimum total pseudobulk count for either edgeR filtering
  method.

- filter_method:

  Character. edgeR expression-filtering rule. `"fixed_support"`
  (default) requires a fixed CPM threshold in a fixed number of
  independent donors and does not use exposure labels. `"filterByExpr"`
  retains the edgeR compatibility option.

- filter_min_cpm:

  Numeric. CPM threshold for `filter_method = "fixed_support"`. Default
  1.

- filter_min_donors:

  Integer. Minimum number of independent donors meeting
  `filter_min_cpm`. Default 5.

- test_features:

  Optional character vector of feature identifiers that defines the
  reported hypothesis family. Library-size calculation, expression
  filtering, normalisation and dispersion estimation still use all
  available features; subsetting occurs only after model fitting. Every
  requested feature must pass the pre-specified expression filter in
  every tested cell type, otherwise the call fails rather than silently
  changing the requested hypothesis family. Currently supported only by
  the edgeR backend; other backends fail explicitly rather than ignoring
  it.

- robust:

  Logical. For the edgeR backend, moderate the quasi-likelihood
  dispersions with robust empirical Bayes (`glmQLFit(robust = TRUE)`).
  Default `TRUE`.

## Value

A `DataFrame` with one row per tested gene and cell type and the same
columns, in the same order, for every backend: `gene`, `celltype`,
`log2FC`, `se`, `statistic`, `pvalue`, `pvalue_underflow_clamped`
(`TRUE` where a p-value of exactly zero was replaced by the smallest
positive double), `padj` (BH within the cell type), `padj_global` (BH
across all tested gene-by-cell-type hypotheses), `logCPM`, `exposure`,
`n_donors`, `n_unexposed` and `n_exposed` (group sizes when the exposure
takes two values, otherwise `NA`), `min_cells` and `median_cells` (cells
per donor in the cell type), and `method`. Rows are ordered by
`padj_global`, then `pvalue`.

`logCPM` is the average log2 count per million reported by the backend:
edgeR's `logCPM` for the edgeR backend, the mean of the voom log-CPM
values for voom-dream, and
[`edgeR::aveLogCPM()`](https://rdrr.io/pkg/edgeR/man/aveLogCPM.html) of
the filtered pseudobulk counts for DESeq2.

Uncertainty fields are backend-specific: edgeR QL reports `se` as
numeric `NA` and `statistic` as the directional signed square root of
the QL F statistic; DESeq2 reports its coefficient standard error and
Wald statistic; voom-dream reports the moderated standard error
(`|log2FC / t|`) and the moderated t statistic.
`metadata(result)$parameters` records the settings, the backend and its
package version, `se_method` and `statistic_type`;
`metadata(result)$diagnostics` holds, for each fitted cell type, the
donors used, the numbers of genes before and after filtering and the
design columns; for the edgeR backend also `dropped_design_columns`, the
covariate columns left out of the design as collinear (each is also
reported in a warning), and `residual_df`, the smallest adjusted
residual degrees of freedom over the fitted genes (`df.residual.adj`).

## Details

The pseudobulk approach aggregates single-cell counts within each (donor
x cell type) stratum, creating one observation per donor per cell type.
This correctly treats cells from the same donor as non-independent,
avoiding the inflated Type I error rates that plague cell-level analyses
(Squair et al. 2021, *Nat Commun*). The pseudobulk counts that enter the
models are available as `SummarizedExperiment` objects from
[`donor_pseudobulk`](donor_pseudobulk.md).

The exposure variable and specified covariates enter one donor-level
design matrix. The default edgeR backend fits the edgeR v4
quasi-likelihood model, `glmQLFit(legacy = FALSE)`, which sets the
negative-binomial dispersion internally from the most highly expressed
genes and moderates the quasi-likelihood dispersions by empirical Bayes;
the exposure coefficient is tested with a quasi-likelihood F-test (Chen
et al. 2025). No separate `estimateDisp()` step is run: from edgeR
4.10.0 the method gives the same result whether or not the data hold
such estimates, whereas earlier edgeR 4.x releases did not. The edgeR
version and `legacy = FALSE` are recorded in
`metadata(result)$parameters`. Complete-transcriptome library sizes are
retained when expression filtering or `test_features` restricts the
reported hypothesis family.

[`edgeR::glmQLFTest`](https://rdrr.io/pkg/edgeR/man/glmQLFTest.html)
does not return a coefficient standard error. Consequently, the edgeR
backend never reconstructs `se` from the ratio of the log-fold change to
the QL statistic. Its `statistic` field is a directional ranking
representation of the QL F statistic, not a Wald z or t statistic, and
must not be used to derive a coefficient confidence interval. The edgeR
result metadata records `se_method = "not_available_edgeR_QL"` and
`statistic_type = "signed_sqrt_qlf"`.

**What is delegated and what is kept.** The raw counts of each donor and
cell type are summed with
[`scrapper::aggregateAcrossCells()`](https://rdrr.io/pkg/scrapper/man/aggregateAcrossCells.html),
and the models are fitted by edgeR, DESeq2 or variancePartition and
limma. Four choices are made by exposomeSC itself, on purpose, because
they define the hypotheses of an exposome-wide study:

1.  The exposure is a donor-level variable from `exposureData(x)`,
    continuous or binary, and enters the design matrix together with the
    covariates; cells are never grouped by exposure.

2.  Gene filtering does not use the exposure. The default edgeR rule
    (`filter_method = "fixed_support"`) keeps genes with at least
    `filter_min_cpm` CPM in `filter_min_donors` donors, and the DESeq2
    and voom-dream backends keep genes detected in at least 20% of
    donors (and in at least two); which genes are tested therefore does
    not depend on the exposure values. By contrast,
    `scran::pseudoBulkDGE()`, for example, filters each label with
    [`edgeR::filterByExpr()`](https://rdrr.io/pkg/edgeR/man/filterByExpr.html)
    using the design matrix.

3.  In the edgeR backend the library sizes are the column sums of the
    complete pseudobulk transcriptome, computed before filtering, and
    are kept when filtering or `test_features` reduces the tested genes.
    The voom-dream backend computes library sizes from the retained
    genes and DESeq2 uses its own size factors.

4.  `padj_global` applies the Benjamini-Hochberg procedure once across
    all tested gene-by-cell-type hypotheses, next to the
    within-cell-type `padj`.

## References

Squair JW et al. (2021). Confronting false discoveries in single-cell
differential expression. *Nat Commun* 12:5692.
[doi:10.1038/s41467-021-25960-2](https://doi.org/10.1038/s41467-021-25960-2)

Chen Y, Chen L, Lun ATL, Baldoni PL, Smyth GK (2025). edgeR v4: powerful
differential analysis of sequencing data with expanded functionality and
improved support for small counts and larger datasets. *Nucleic Acids
Research* 53(2):gkaf018.
[doi:10.1093/nar/gkaf018](https://doi.org/10.1093/nar/gkaf018)

Hoffman GE, Roussos P (2021). Dream: powerful differential expression
analysis for repeated measures designs. *Bioinformatics* 37(2):192–201.
[doi:10.1093/bioinformatics/btaa687](https://doi.org/10.1093/bioinformatics/btaa687)

Hoffman GE, Schadt EE (2016). variancePartition: interpreting drivers of
variation in complex gene expression studies. *BMC Bioinformatics*
17:483.
[doi:10.1186/s12859-016-1323-z](https://doi.org/10.1186/s12859-016-1323-z)

## See also

[`donor_pseudobulk`](donor_pseudobulk.md) for the pseudobulk counts,
[`run_multi_exwas`](run_multi_exwas.md) for several exposures.

## Examples

``` r
library(SingleCellExperiment)
library(S4Vectors)
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
result <- run_sc_exwas(scee, exposure = "E1",
    celltype_col = "cell_type", min_cells = 5L)
head(result)
#> DataFrame with 6 rows and 17 columns
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
metadata(result)$parameters$statistic_type
#> [1] "signed_sqrt_qlf"
```
