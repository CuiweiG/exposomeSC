# Formal mixture analysis via quantile g-computation

Wraps
[`qgcomp::qgcomp.noboot()`](https://rdrr.io/pkg/qgcomp/man/qgcomp.glm.noboot.html)
or
[`qgcomp::qgcomp.boot()`](https://rdrr.io/pkg/qgcomp/man/qgcomp.glm.boot.html)
on pseudobulk data for cell-type-specific exposure mixture analysis.
Unlike [`run_sc_mixture`](run_sc_mixture.md) (simplified screening), it
reports a confidence interval and p-value for the overall mixture effect
and, without bootstrapping, the positive and negative exposure weights
estimated by qgcomp.

## Usage

``` r
run_mixture_qgcomp(
  x,
  exposures,
  celltype,
  celltype_col = "cell_type",
  target_genes = NULL,
  covariates = NULL,
  q = 4L,
  bootstrap = FALSE,
  B = 200L,
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
  goes into the model formula, and cannot be `.`, `...` or a name such
  as `..1`, which R reads specially there. The names `y`, `weights` and
  `id__` cannot be used: the response is called `y`, and qgcomp adds
  columns called `weights` and `id__` to the model data. With
  `bootstrap = FALSE` a name cannot contain `_q` either, which qgcomp
  removes from the names of the weights it returns.

- celltype:

  Character. Target cell type.

- celltype_col:

  Character. `colData` column with cell type labels. Default
  `"cell_type"`.

- target_genes:

  Character vector (optional). Genes to summarise. Default: top 20 most
  variable.

- covariates:

  Character vector (optional). Columns of `exposureData`; an unknown
  name is an error, and so is an exposure listed here. Names must be
  syntactic and cannot be `y`, `weights`, `id__`, `.`, `...` or a name
  such as `..1`, as for `exposures`.

- q:

  Integer. Number of quantile bins. Default 4.

- bootstrap:

  Logical. If `TRUE`, use `qgcomp.boot()` for bootstrap confidence
  intervals. Default `FALSE` (faster; a t interval on the residual
  degrees of freedom, see Details).

- B:

  Integer. Number of bootstrap replicates (only if `bootstrap = TRUE`).
  Default 200.

- min_cells:

  Integer. Default 10.

## Value

A one-row `DataFrame` with columns `celltype`, `mixture_effect` (the
overall mixture effect \\\psi\\), `mixture_se`, `ci_lower` and
`ci_upper` (95% confidence interval for \\\psi\\), `mixture_pvalue` (for
the interval and the p-value see Details), `n_donors`, `n_genes` and
`method` (`"qgcomp.noboot"` or `"qgcomp.boot"`). `metadata(result)`
holds:

- positive_weights:

  Named numeric. Exposures that increase the response, with relative
  weights; `NULL` when `bootstrap = TRUE`, because `qgcomp.boot()` does
  not estimate weights.

- negative_weights:

  Named numeric. Exposures that decrease the response; `NULL` when
  `bootstrap = TRUE`.

- target_genes:

  The genes whose mean log-CPM is the response.

- fit:

  The `qgcomp` model object (for plotting and further analysis).

## Details

Quantile g-computation (Keil et al. 2020) estimates the overall effect
of jointly increasing all exposures by one quantile, and without
bootstrapping decomposes it into positive and negative contributions.
The weights are interpretable only when the exposure effects are of a
similar shape (Keil et al. 2020).

The pseudobulk response is the mean log-CPM of the target genes over
donors with complete exposure and covariate data.

Without bootstrapping, `mixture_pvalue` is a two-sided t test of
\\\psi\\ on the residual degrees of freedom of the fitted model, and
`ci_lower` and `ci_upper` are the matching t interval. `qgcomp.noboot()`
itself uses the null degrees of freedom minus the number of exposures
for its p-value, which leaves the covariates out, and a normal interval;
its own values are kept in `metadata(result)$fit`. With
`bootstrap = TRUE` the p-value and the interval are those of
`qgcomp.boot()`; the call stops when `qgcomp.boot()` returns no p-value
for \\\psi\\.

**Comparison with `run_sc_mixture`:**

|                            |                |                               |
|----------------------------|----------------|-------------------------------|
| Feature                    | run_sc_mixture | run_mixture_qgcomp            |
| CI for overall effect      | No             | Yes (t interval or bootstrap) |
| Pos/neg weights            | No             | Yes (without bootstrap)       |
| P-value for overall effect | No             | Yes                           |
| Speed                      | Fast           | Moderate                      |
| Minimum n                  | 5              | 10+ recommended               |

## References

Keil AP et al. (2020). A quantile-based g-computation approach to
addressing the effects of exposure mixtures. *Environ Health Perspect*
128:047004. [doi:10.1289/EHP5838](https://doi.org/10.1289/EHP5838)

## Examples

``` r
# Requires the qgcomp package
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
if (requireNamespace("qgcomp", quietly = TRUE)) {
    mix <- run_mixture_qgcomp(scee,
        exposures = c("E1", "E2", "E3"),
        celltype = "Mono", celltype_col = "cell_type",
        min_cells = 3L)
    mix
    metadata(mix)$positive_weights
    metadata(mix)$negative_weights
}
#> E1 
#>  1 
```
