# Cell-type-specific mixture exposure analysis

Estimates the relative contribution of multiple correlated exposures to
gene expression within a specific cell type using a simplified
quantile-scored linear model on pseudobulk data.

## Usage

``` r
run_sc_mixture(
  x,
  exposures,
  celltype,
  celltype_col = "cell_type",
  target_genes = NULL,
  covariates = NULL,
  q = 4L,
  min_cells = 10L
)
```

## Arguments

- x:

  A
  [`SingleCellExposomeExperiment`](SingleCellExposomeExperiment-class.md).
  Each cell's donor is read from the `colData` column named by
  [`sampleColumn`](SCEE-accessors.md)`(x)`.

- exposures:

  Character vector. Exposure variable names, distinct. Each name must be
  a syntactic R name (see
  [`make.names`](https://rdrr.io/r/base/make.names.html)) because it
  goes into the model formula, and cannot be `y`, the name of the
  response. Nor can it be `.`, `...` or a name such as `..1`, which R
  reads specially in a formula.

- celltype:

  Character. Which cell type to analyse.

- celltype_col:

  Character. `colData` column with cell type labels. Default
  `"cell_type"`.

- target_genes:

  Character vector (optional). Genes to test. Default: top 20 most
  variable genes.

- covariates:

  Character vector (optional). Additional covariates from
  `exposureData`; a name that is not a column of `exposureData` is an
  error, and so is an exposure listed here. Names must be syntactic and
  cannot be `y`, `.`, `...` or a name such as `..1`, as for `exposures`.

- q:

  Integer. Number of quantile bins, at least 2. Default 4. Tied exposure
  values (for example at a limit of detection) can give fewer bins.

- min_cells:

  Integer. Minimum cells per donor. Default 10.

## Value

A `DataFrame` with one row per target gene with a fitted model: `gene`
and `mixture_effect` (the sum of that gene's estimable quantile-scored
exposure coefficients; `NA` when none is estimable). `metadata(result)`
holds:

- weights:

  Named numeric vector. Relative contribution of each exposure (the
  weights that are not `NA` sum to 1): the mean over target genes of the
  absolute coefficients of the quantile-scored exposures, normalised.
  When the quantile scores of exposures are collinear, for example for
  two exposures with the same ranks,
  [`lm()`](https://rdrr.io/r/stats/lm.html) estimates a coefficient only
  for the one listed first in `exposures`; it carries their joint
  effect, and the later ones get the weight `NA`, with a warning. When
  no coefficient can be estimated, every weight is `NA`.

- mixture_coef:

  Overall mixture coefficient (sum of quantile-scored exposure
  coefficients, averaged across target genes); `NA` when no coefficient
  can be estimated.

- method:

  Character. Always `"quantile_linear"`.

- celltype:

  Character. Cell type analysed.

- n_donors:

  Integer. Donors with complete exposure and covariate data used in the
  analysis.

- n_genes:

  Integer. Genes with a fitted model.

## Details

This is a **simplified screening method**, not a formal quantile-based
g-computation (Keil et al. 2020) or weighted quantile sum regression
(Carrico et al. 2015). The approach:

1.  Donors with a missing exposure or covariate are removed, and
    exposures are scored into quantile bins (1 to `q`; fewer when
    quantiles are tied).

2.  For each target gene, a linear model is fit on pseudobulk log-CPM
    with all quantile-scored exposures and covariates as predictors.

3.  Weights are the absolute coefficients of the exposures that have
    one, averaged across target genes and normalised to sum to 1; an
    exposure without an estimable coefficient gets the weight `NA` (see
    Value). Because every exposure is on the same quantile scale, the
    coefficients are directly comparable.

**Limitations:**

- No bootstrap confidence intervals on weights.

- Does not separate positive/negative weight directions as in formal
  qgcomp.

- Log-CPM may have heteroscedasticity at low counts.

For rigorous mixture analyses suitable for publication, consider
[`qgcomp::qgcomp.noboot()`](https://rdrr.io/pkg/qgcomp/man/qgcomp.glm.noboot.html)
or `gWQS::gwqs()` applied to the pseudobulk counts returned by
[`donor_pseudobulk`](donor_pseudobulk.md), or
[`run_mixture_qgcomp`](run_mixture_qgcomp.md).

## References

Keil AP et al. (2020). A quantile-based g-computation approach to
addressing the effects of exposure mixtures. *Environ Health Perspect*
128:047004. [doi:10.1289/EHP5838](https://doi.org/10.1289/EHP5838)

Carrico C et al. (2015). Characterization of weighted quantile sum
regression for highly correlated data in a risk analysis setting. *J
Agric Biol Environ Stat* 20:100-120.

## Examples

``` r
library(SingleCellExperiment); library(S4Vectors)
set.seed(1)
counts <- matrix(rpois(5000, 8), nrow = 50,
    dimnames = list(paste0("G", 1:50), paste0("c", 1:100)))
sce <- SingleCellExperiment(assays = list(counts = counts),
    colData = DataFrame(cell_id = paste0("c", 1:100),
        donor_id = rep(paste0("D", 1:10), each = 10),
        cell_type = rep(c("Mono", "NK"), 50)))
exp_mat <- matrix(rnorm(30), nrow = 10,
    dimnames = list(paste0("D", 1:10), c("E1", "E2", "E3")))
scee <- build_scee(sce, exp_mat, sample_col = "donor_id")
mix <- run_sc_mixture(scee, exposures = c("E1", "E2", "E3"),
    celltype = "Mono", celltype_col = "cell_type",
    min_cells = 3L)
metadata(mix)$weights
#>        E1        E2        E3 
#> 0.3857350 0.2803696 0.3338954 
```
