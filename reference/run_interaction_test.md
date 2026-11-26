# Test exposure-by-celltype interaction

Formally tests whether an exposure effect differs across cell types
using a pseudobulk interaction model. This prevents the common mistake
of concluding "the exposure affects cell type A but not B" based on
separate per- celltype p-values (Gelman & Stern 2006).

## Usage

``` r
run_interaction_test(
  x,
  exposure,
  celltype_col = "cell_type",
  target_genes = NULL,
  covariates = NULL,
  min_cells = 10L,
  min_donors = 5L
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

- celltype_col:

  Character. `colData` column with cell type labels. Default
  `"cell_type"`.

- target_genes:

  Character vector (optional). Genes to test. Default: top 50 most
  variable genes across all cell types.

- covariates:

  Character vector (optional). Columns of `exposureData` adjusted for;
  listing the exposure here is an error.

- min_cells:

  Integer. Min cells per donor per celltype. Default 10.

- min_donors:

  Integer. Min donors per celltype. Default 5.

## Value

A `DataFrame` with columns: gene, p_interaction (F-test for the
exposure-by-cell-type terms), df_interaction (numerator degrees of
freedom), padj_interaction (BH-corrected), n_celltypes, n_donors.

## Details

For each gene, pseudobulk log-CPM values of all cell types are stacked,
one row per donor and cell type, and two linear models are compared with
an F-test:


    reduced: y ~ donor + celltype + covariates:celltype
    full:    y ~ donor + celltype + covariates:celltype + exposure:celltype

Donor fixed effects absorb every donor-level quantity, including the
exposure and covariate main effects, and the shared donor component of
expression in different cell types, so the test compares exposure slopes
between cell types through within-donor contrasts rather than treating
cell types from the same donor as independent. Covariate-by-cell-type
terms allow covariate effects to differ between cell types. Donors
observed in only one cell type do not inform the interaction.

Comparing significance levels across cell types is not the same as
testing for a difference between them (Gelman & Stern 2006, Nieuwenhuis
et al. 2011); this function tests the difference directly.

## References

Gelman A, Stern H (2006). The difference between 'significant' and 'not
significant' is not itself statistically significant. *Am Stat*
60:328-331.

Nieuwenhuis S et al. (2011). Erroneous analyses of interactions in
neuroscience. *Nat Neurosci* 14:1105-1107.
[doi:10.1038/nn.2886](https://doi.org/10.1038/nn.2886)

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
ix <- run_interaction_test(scee, exposure = "E1",
    celltype_col = "cell_type", min_cells = 3L, min_donors = 3L)
head(ix)
#> DataFrame with 6 rows and 6 columns
#>          gene p_interaction df_interaction n_celltypes  n_donors
#>   <character>     <numeric>      <integer>   <integer> <integer>
#> 1         G22      0.690566              1           2        10
#> 2         G19      0.226742              1           2        10
#> 3         G30      0.361459              1           2        10
#> 4         G47      0.109298              1           2        10
#> 5         G10      0.958973              1           2        10
#> 6         G13      0.226436              1           2        10
#>   padj_interaction
#>          <numeric>
#> 1         0.917539
#> 2         0.872086
#> 3         0.890187
#> 4         0.872086
#> 5         0.958973
#> 6         0.872086
```
