# Cell-type-specific exposure dose-response analysis

Fits linear and polynomial dose-response models for each gene within a
cell type, testing whether the exposure-expression relationship is
nonlinear. Linear exposure-response models are common, but some
exposures, such as endocrine-disrupting chemicals, show non-monotonic
dose responses (Vandenberg et al. 2012).

## Usage

``` r
run_dose_response(
  x,
  exposure,
  celltype,
  celltype_col = "cell_type",
  target_genes = NULL,
  max_degree = 3L,
  min_cells = 10L
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
  variable genes.

- max_degree:

  Integer. Maximum polynomial degree to test. Default 3 (cubic). Higher
  values risk overfitting with small sample sizes.

- min_cells:

  Integer. Minimum cells per donor. Default 10.

## Value

A `DataFrame` with columns: gene, best_model ("linear", "quadratic",
"cubic"), linear_coef, AIC_linear, AIC_best, p_nonlinear (F-test of
`nonlinear_model` against the linear model), nonlinear_model, celltype,
n_donors (donors with at least `min_cells` cells and a non-missing
exposure; donors with a missing exposure are left out, with a message).

## Details

For each gene, the function fits:

1.  `y ~ exposure` (linear)

2.  `y ~ poly(exposure, 2)` (quadratic)

3.  `y ~ poly(exposure, 3)` (cubic, if max_degree \>= 3)

on pseudobulk log-CPM values. `best_model` is the model with the lowest
AIC and is descriptive.

`p_nonlinear` is an F-test of the highest-degree polynomial that could
be fitted (cubic when `max_degree = 3` and at least six donors are
available, otherwise quadratic) against the linear model. It is fixed in
advance and does not depend on the AIC selection, because testing the
model chosen by AIC would overstate the evidence for nonlinearity.

## References

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
dr <- run_dose_response(scee, exposure = "E1",
    celltype = "Mono", celltype_col = "cell_type",
    min_cells = 3L,
    max_degree = 2L)
head(dr)
#> DataFrame with 6 rows and 9 columns
#>          gene  best_model linear_coef AIC_linear  AIC_best p_nonlinear
#>   <character> <character>   <numeric>  <numeric> <numeric>   <numeric>
#> 1         G10      linear  -0.0674192   13.88972  13.88972   0.8638341
#> 2         G19      linear  -0.1881040    8.66395   8.66395   0.7741065
#> 3         G32   quadratic   0.0796200   11.05071  10.97180   0.2440552
#> 4         G13   quadratic  -0.0080711   11.55727   8.05291   0.0577554
#> 5         G49   quadratic   0.1016717    8.57050   8.37717   0.2314692
#> 6         G17      linear   0.0862627    7.92353   7.92353   0.9710102
#>   nonlinear_model    celltype  n_donors
#>       <character> <character> <integer>
#> 1       quadratic        Mono        10
#> 2       quadratic        Mono        10
#> 3       quadratic        Mono        10
#> 4       quadratic        Mono        10
#> 5       quadratic        Mono        10
#> 6       quadratic        Mono        10
```
