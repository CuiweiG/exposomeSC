# GAM-based dose-response analysis

Fits generalized additive models (GAMs) for cell-type- specific
dose-response curves, allowing the data to determine the shape of the
exposure-expression relationship without pre-specifying a polynomial
degree. Unlike [`run_dose_response`](run_dose_response.md), the shape is
not restricted to a polynomial.

## Usage

``` r
run_dose_response_gam(
  x,
  exposure,
  celltype,
  celltype_col = "cell_type",
  target_genes = NULL,
  k = NULL,
  min_cells = 10L,
  loocv = TRUE
)
```

## Arguments

- x:

  A
  [`SingleCellExposomeExperiment`](SingleCellExposomeExperiment-class.md).
  Each cell's donor is read from the `colData` column named by
  [`sampleColumn`](SCEE-accessors.md)`(x)`.

- exposure:

  Character. Exposure variable name.

- celltype:

  Character. Target cell type.

- celltype_col:

  Character. `colData` column with cell type labels. Default
  `"cell_type"`.

- target_genes:

  Character vector (optional). Genes to test. Default: top 50 most
  variable.

- k:

  Integer. Maximum basis dimension for the smooth term, at least 4 and
  at most the number of distinct exposure values. Default
  `min(n_donors - 1, 10)`, capped at that number. Higher k allows more
  wiggliness; the REML smoothing penalty limits overfitting. The
  exposure needs at least 4 distinct values, because the test of
  nonlinearity needs `k` of at least 4; for an exposure with 3 values
  use [`run_dose_response`](run_dose_response.md).

- min_cells:

  Integer. Min cells per donor. Default 10.

- loocv:

  Logical. If `TRUE`, compute leave-one-out cross-validated R² for the
  linear and GAM models. Default `TRUE`.

## Value

A `DataFrame` with columns: gene, edf (effective degrees of freedom for
smooth term), p_smooth (significance of smooth term vs intercept),
p_nonlinear (approximate test that a fully penalised smooth added to a
linear exposure term is zero, that is, of departure from linearity; see
Details), AIC_linear, AIC_gam, R2_linear, R2_gam, R2_loocv_linear (if
loocv=TRUE), R2_loocv_gam, deviance_explained, celltype, n_donors
(donors with at least `min_cells` cells and a non-missing exposure;
donors with a missing exposure are left out, with a message).

## Details

For each gene, two models are fit:

1.  Linear: `y ~ exposure` (OLS)

2.  GAM: `y ~ s(exposure, k = k)` (penalized thin-plate regression
    spline via [`mgcv::gam`](https://rdrr.io/pkg/mgcv/man/gam.html))

The GAM smooth is estimated using restricted maximum likelihood (REML),
which selects the smoothness penalty. The effective degrees of freedom
(edf) indicate the complexity: edf close to 1 means approximately
linear, and larger values mean more curvature.

`p_nonlinear` comes from a third model,
`y ~ exposure + s(exposure, bs = "tp", m = c(2, 0))`, whose smooth has
no unpenalised null space, so the linear trend is carried by the
parametric term and the smooth can only represent departures from
linearity. Its p-value is the approximate Wald-type test for smooth
terms of Wood (2013). Comparing the linear and GAM fits by an F-test on
their residual deviances is not valid here: the GAM's effective degrees
of freedom are estimated and often close to one.

Leave-one-out cross-validated R² (loocv) provides an honest estimate of
out-of-sample prediction accuracy, unlike apparent R² which inflates
with model complexity. This is critical for small-sample dose-response
where overfitting is a real risk.

**Why GAM over polynomial:**

- Polynomials impose a global parametric form

- GAMs adapt locally – can capture threshold effects, plateaus, U-shapes

- Penalised smoothness selection limits overfitting

- edf provides an interpretable complexity measure

## References

Wood SN (2017). Generalized Additive Models: An Introduction with R. 2nd
ed. Chapman and Hall/CRC.

Wood SN (2013). On p-values for smooth components of an extended
generalized additive model. *Biometrika* 100:221-228.
[doi:10.1093/biomet/ass048](https://doi.org/10.1093/biomet/ass048)

Vandenberg LN et al. (2012). Hormones and endocrine- disrupting
chemicals: low-dose effects and nonmonotonic dose responses. *Endocr
Rev* 33:378-455.
[doi:10.1210/er.2011-1050](https://doi.org/10.1210/er.2011-1050)

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
exp_mat <- matrix(rnorm(20), nrow = 10,
    dimnames = list(paste0("D", 1:10), c("E1", "E2")))
scee <- build_scee(sce, exp_mat, sample_col = "donor_id")
if (requireNamespace("mgcv", quietly = TRUE)) {
    dr <- run_dose_response_gam(scee, exposure = "E1",
        celltype = "Mono", celltype_col = "cell_type",
        min_cells = 3L,
        target_genes = paste0("G", 1:3))
    head(dr)
}
#> DataFrame with 3 rows and 13 columns
#>          gene       edf  p_smooth p_nonlinear AIC_linear   AIC_gam   R2_linear
#>   <character> <numeric> <numeric>   <numeric>  <numeric> <numeric>   <numeric>
#> 1          G1   1.00000 0.3073461  0.66353346   -2.90307  -2.90306 0.129359514
#> 2          G2   2.73573 0.0144659  0.00371955    6.25301  -6.40253 0.012230856
#> 3          G3   1.00001 0.9597158  0.61217941    3.01409   3.01411 0.000340272
#>       R2_gam R2_loocv_linear R2_loocv_gam deviance_explained    celltype
#>    <numeric>       <numeric>    <numeric>          <numeric> <character>
#> 1  0.0205293       -0.402119    -0.402121        0.129359722        Mono
#> 2  0.7429129       -0.688037     0.420624        0.821059614        Mono
#> 3 -0.1246175       -0.807911    -0.872429        0.000340886        Mono
#>    n_donors
#>   <integer>
#> 1        10
#> 2        10
#> 3        10
```
