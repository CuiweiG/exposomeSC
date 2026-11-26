# Cell-Level Gene-Protein Coupling Analysis

Computes within-donor Spearman partial correlations between gene
expression and protein abundance at the single-cell level, then uses
random-effects meta-regression to test whether coupling strength varies
with a donor-level exposure.

## Usage

``` r
run_cell_coupling(
  x,
  genes,
  proteins,
  exposure,
  celltype,
  celltype_col = "cell_type",
  altexp_name = "CITE",
  confounders = NULL,
  min_cells = 30L,
  min_donors = 10L,
  meta_moderators = NULL,
  var_method = c("asymptotic", "bootstrap"),
  n_boot = 200L
)
```

## Arguments

- x:

  A `SingleCellExposomeExperiment` object. Each cell's donor is read
  from the `colData` column named by
  [`sampleColumn`](SCEE-accessors.md)`(x)`. Must contain an altExp named
  `altexp_name` with cell-level protein counts.

- genes:

  Character vector of gene names to test.

- proteins:

  Character vector of protein names to test (must be present in
  `altExp`).

- exposure:

  Character; name of the exposure variable in `exposureData(x)`.

- celltype:

  Character; cell type to analyse.

- celltype_col:

  Character; column in `colData` containing cell type labels. Default
  `"cell_type"`.

- altexp_name:

  Character; name of the altExp containing protein data. Default
  `"CITE"`.

- confounders:

  Character vector of column names in `colData(x)` to partial out (e.g.,
  `c("S.Score", "G2M.Score", "percent.mt")`). Library sizes (nCount_RNA,
  nCount_protein) are always controlled regardless of this parameter.
  Each confounder is ranked within a donor, so it must be numeric,
  logical, an ordered factor or a label with at most two values; an
  unordered label with more categories, such as a cell-cycle phase, is
  an error and must be given as 0/1 indicator columns for all its
  categories but one, leaving out a category that every donor has cells
  of (for G1, S and G2M, for example S and G2M when every donor has G1
  cells). The names `gene`, `protein`, `lib_rna` and `lib_prot` are used
  for the variables of the partial correlation and cannot be confounder
  names. Nor can a confounder be a column that is not passed on as one
  column under its own name: a matrix or nested table with several
  columns or with a column name other than that of the confounder, or a
  column with the name of an argument of
  [`data.frame`](https://rdrr.io/r/base/data.frame.html) (`row.names`,
  `check.rows`, `check.names`, `fix.empty.names` or `stringsAsFactors`).
  Default `NULL` (library sizes only).

- min_cells:

  Integer; minimum cells per donor to include. Default 30.

- min_donors:

  Integer; minimum number of donors with enough cells (at least
  `min_cells`, and at least 7 plus the number of confounders, because
  the two library sizes are always adjusted for) and values of the
  exposure and of every moderator. A gene-protein pair also needs a
  partial correlation in that many donors. Default 10.

- meta_moderators:

  Character vector or one-sided formula of additional donor-level
  moderators for the meta-regression (e.g., `"age"`, `c("age", "sex")`
  or `~ age + sex`). These must be columns in `exposureData(x)` with
  syntactic R names (see
  [`make.names`](https://rdrr.io/r/base/make.names.html)) other than
  `.`, `...` and names such as `..1`, which R reads specially in a
  formula; a formula may only list such columns, without transformations
  or interactions. The exposure itself cannot be a moderator, and the
  names `yi`, `vi`, `donor_exp`, `intrcpt` and `meta_formula`, which the
  meta-regression uses for its own variables, cannot be moderator names.
  Default `NULL` (exposure only).

- var_method:

  Character. Method for estimating the sampling variance of Fisher
  z-transformed correlations. `"asymptotic"` (default) uses the formula
  \\1.06 / (n - q - 3)\\. `"bootstrap"` resamples cells within each
  donor to compute empirical variance, which is more accurate for small
  cell counts or non-normal data but slower. The bootstrap uses the
  random number generator, so set a seed for reproducible variances.

- n_boot:

  Integer. Number of bootstrap replicates per donor when
  `var_method = "bootstrap"`. Default 200.

## Value

A `DataFrame` with one row per gene-protein pair (zero rows when too few
donors qualify). A pair whose exposure does not vary among its donors,
or is collinear with the moderators, is not tested, with a warning that
counts such pairs:

- gene, protein:

  Feature names.

- n_donors:

  Number of donors in the meta-regression: donors with enough cells, a
  partial correlation and values of the exposure and of every moderator.

- beta0, se0, pval0:

  Intercept (baseline coupling).

- beta1, se1, pval1:

  Exposure coefficient: the slope, across donors, of the Fisher z of the
  coupling on the exposure.

- tau2:

  Between-donor heterogeneity variance (REML).

- I2:

  Heterogeneity percentage.

- mean_r:

  Mean partial correlation across donors.

- mean_n_cells:

  Mean cells per donor.

- n_confounders:

  Number of conditioning variables requested (library sizes plus
  confounders); variables constant within a donor are dropped for that
  donor.

## Details

For each gene-protein pair, the method proceeds in two stages:

**Stage 1 (within-donor):** For each donor \\d\\ with at least
`min_cells` cells of the specified cell type, compute the Spearman
partial correlation between gene expression and protein abundance,
controlling for library sizes and any user-specified confounders:
\$\$r_d = \rho_S(\text{gene}\_i, \text{protein}\_j \mid
\text{nCount\\RNA}, \text{nCount\\protein}, \text{confounders})\$\$ Then
apply Fisher's z-transform: \$\$z_d = \text{arctanh}(r_d), \quad
\text{Var}(z_d) \approx 1.06 / (n_d - q - 3)\$\$ where \\q\\ is the
number of conditioning variables (2 + length of `confounders`, less any
that are constant within the donor and therefore dropped). The factor
1.06 corrects for the efficiency loss of rank-based correlation (Fieller
et al. 1957).

**Stage 2 (across-donor meta-regression):** Fit a random-effects model
via REML with Knapp-Hartung adjustment: \$\$z_d = \beta_0 + \beta_1
\cdot \text{exposure}\_d + \sum_k \gamma_k \cdot
\text{moderator}\_{dk} + u_d + \varepsilon_d\$\$ where \\\varepsilon_d
\sim N(0, v_d)\\ (known sampling variance) and \\u_d \sim N(0, \tau^2)\\
(between-donor heterogeneity). The Knapp-Hartung adjustment provides
more accurate confidence intervals for small numbers of studies
(donors), using a t-distribution with \\D - p - 1\\ degrees of freedom
rather than a z-distribution.

**Confounder control:** Cell-level confounders (cell cycle,
mitochondrial content, ambient RNA scores) are partialed out in Stage 1,
ensuring that the correlation reflects genuine gene-protein coupling
rather than shared technical or biological confounds. Donor-level
moderators (age, sex) are included in Stage 2.

## References

Viechtbauer W (2010). Conducting meta-analyses in R with the metafor
package. *J Stat Softw*, 36(3), 1-48.
[doi:10.18637/jss.v036.i03](https://doi.org/10.18637/jss.v036.i03)

Fieller EC, Hartley HO, Pearson ES (1957). Tests for rank correlation
coefficients. I. *Biometrika*, 44, 470-481.

Kim S (2015). ppcor: An R Package for a Fast Calculation to Semi-partial
Correlation Coefficients. *Commun Stat Appl Methods*, 22(6), 665-674.

## Examples

``` r
set.seed(1)
donors <- sprintf("D%02d", 1:12)
donor <- rep(donors, each = 60)
gene <- matrix(stats::rpois(2 * length(donor), 10), nrow = 2,
    dimnames = list(c("Gene1", "Gene2"), paste0("c", seq_along(donor))))
protein <- matrix(stats::rpois(2 * length(donor), 5), nrow = 2,
    dimnames = list(c("Prot1", "Prot2"), colnames(gene)))
sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = gene),
    colData = S4Vectors::DataFrame(donor_id = donor, cell_type = "T"))
SingleCellExperiment::altExp(sce, "CITE") <-
    SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = protein))
exp_mat <- matrix(seq(0, 2, length.out = 12), ncol = 1,
    dimnames = list(donors, "exposure"))
scee <- build_scee(sce, exp_mat, sample_col = "donor_id")
if (requireNamespace("metafor", quietly = TRUE) &&
    requireNamespace("ppcor", quietly = TRUE)) {
    run_cell_coupling(scee, genes = "Gene1", proteins = "Prot1",
        exposure = "exposure", celltype = "T",
        min_cells = 20L, min_donors = 10L)
}
#> [exposomeSC] Cell coupling: T, 1 genes x 1 proteins = 1 pairs, 2 confounders
#> [exposomeSC] 1 pairs tested for T
#> DataFrame with 1 row and 14 columns
#>          gene     protein  n_donors      beta0       se0     pval0     beta1
#>   <character> <character> <integer>  <numeric> <numeric> <numeric> <numeric>
#> 1       Gene1       Prot1        12 -0.0991707 0.0900614  0.296631 0.0832626
#>         se1     pval1       tau2        I2     mean_r mean_n_cells
#>   <numeric> <numeric>  <numeric> <numeric>  <numeric>    <numeric>
#> 1  0.076281  0.300642 0.00823432   29.9353 -0.0150025           60
#>   n_confounders
#>       <integer>
#> 1             2
```
