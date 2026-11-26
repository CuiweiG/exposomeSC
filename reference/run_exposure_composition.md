# Test exposure associations with donor-level cell composition

Tests whether an exposure is associated with the joint cell-composition
profile while retaining donors, rather than cells or repeated biological
samples, as the independent units of inference.

## Usage

``` r
run_exposure_composition(
  x,
  exposure,
  celltype_col = "cell_type",
  replicate_col = NULL,
  covariates = NULL,
  min_cells = 10L,
  support_min_donors = 10L,
  support_min_total_cells = 100L,
  celltype_universe = NULL,
  exclude_celltypes = character(),
  zero_prior = 0.5,
  replicate_aggregation = c("equal_replicate", "pooled_counts"),
  n_permutations = 0L,
  permutation_strata = NULL,
  seed = NULL
)
```

## Arguments

- x:

  A
  [`SingleCellExposomeExperiment`](SingleCellExposomeExperiment-class.md).
  Each cell's donor is read from the `colData` column named by
  [`sampleColumn`](SCEE-accessors.md)`(x)`.

- exposure:

  Character scalar naming the exposure variable.

- celltype_col:

  Character scalar naming the cell-type column in `colData(x)`.

- replicate_col:

  Optional character scalar naming biological samples nested within
  donors. If `NULL`, each donor is treated as one sample.

- covariates:

  Character vector naming numeric design columns in `exposureData(x)`.

- min_cells:

  Minimum number of cells from the retained universe per biological
  sample.

- support_min_donors:

  Minimum number of donors in which a cell type must be observed.
  Support is evaluated without reference to exposure values.

- support_min_total_cells:

  Minimum total number of cells required for a cell type, also evaluated
  without reference to exposure values.

- celltype_universe:

  Optional pre-specified character vector of cell types eligible for
  analysis.

- exclude_celltypes:

  Character vector of annotation labels excluded before support
  filtering.

- zero_prior:

  Positive symmetric Dirichlet prior. The default 0.5 is the Jeffreys
  prior.

- replicate_aggregation:

  Either `"equal_replicate"`, the primary approach that averages
  log-ratio coordinates across samples within each donor, or
  `"pooled_counts"`, which pools cells within donor and is intended as a
  sensitivity analysis.

- n_permutations:

  Number of Freedman–Lane residual permutations for the omnibus test.
  Zero reports the asymptotic partial Pillai test only.

- permutation_strata:

  Optional exchangeability strata for permutation. Supply a named vector
  indexed by donor ID for study-stratified inference.

- seed:

  `NULL` (default) or a single integer, used only for permutation. With
  `NULL` the permutations draw from the caller's random-number stream,
  so [`set.seed()`](https://rdrr.io/r/base/Random.html) before the call
  makes them reproducible. With an integer, that seed is used inside the
  call and the caller's random-number state is restored before return.

## Value

A `DataFrame` with one HC3-adjusted cell-type-versus-rest log-ratio
contrast per supported cell type. `metadata(result)` stores the partial
Pillai omnibus test, support audit, design diagnostics, and the
donor-level model matrices used for inference. `padj` is the Holm
family-wise adjusted p-value; `padj_bh` is supplied as a secondary
false-discovery-rate summary.

## Details

Counts are formed for every biological sample by supported cell type.
Zeros are handled through posterior expected log proportions under a
symmetric Dirichlet(1/2) prior. Orthonormal isometric log-ratio
coordinates define a basis-invariant partial Pillai omnibus test.
Repeated samples are averaged within donor in log-ratio space, giving
every donor equal regression weight. Cell-type contrasts are exposure
coefficients for
`log2(cell type / geometric mean of all other supported types)` with HC3
heteroscedasticity-consistent standard errors. Nominal 95 percent
confidence intervals and Bonferroni simultaneous 95 percent family-wise
intervals are both returned and labelled explicitly.

The method estimates relative composition within the supported, retained
cell-type universe. It does not identify absolute cell abundance, and
the Dirichlet prior handles observed zeros but does not model
dissociation, capture, annotation, or tissue-sampling selection
mechanisms.

**Relation to other methods.** `speckle::propeller()` (Phipson et al.
2022) transforms the cell-type proportions of each sample (logit or
arcsine square root) and fits a limma linear model for each cell type to
compare groups of samples. `crumblr::crumblr()` (Hoffman and Roussos
2026) returns centred log-ratio transformed counts with
observation-level precision weights for weighted linear or linear mixed
models. sccomp (Mangiola et al. 2023) fits a Bayesian beta-binomial
model of composition and variability. voomCLR (Assefa et al. 2026),
distributed outside Bioconductor, applies the centred log-ratio
transform with observation-level variances from a voom-type
mean-variance trend and offers a bias correction of the fitted
coefficients. `run_exposure_composition()` differs from these in three
respects: by default (`replicate_aggregation = "equal_replicate"`)
repeated samples of a donor are averaged in log-ratio space, so that
each donor carries equal weight; the omnibus test of the whole
composition uses orthonormal isometric log-ratio coordinates and does
not depend on which basis is chosen; and the permutation test can be
restricted to exchangeability strata such as studies
(`permutation_strata`). Its per-cell-type contrasts use HC3 standard
errors, without precision weights or empirical Bayes moderation of the
variances.

## References

Aitchison J (1986). The Statistical Analysis of Compositional Data.
Chapman and Hall.

Freedman D, Lane D (1983). A nonstochastic interpretation of reported
significance levels. Journal of Business and Economic Statistics, 1,
292–298.

Phipson B, Sim CB, Porrello ER, Hewitt AW, Powell J, Oshlack A (2022).
propeller: testing for differences in cell type proportions in single
cell data. *Bioinformatics* 38(20):4720–4726.
[doi:10.1093/bioinformatics/btac582](https://doi.org/10.1093/bioinformatics/btac582)

Hoffman GE, Roussos P (2026). Fast, flexible analysis of differences in
cellular composition with crumblr. *Nat Commun* 17:9598.
[doi:10.1038/s41467-026-75681-7](https://doi.org/10.1038/s41467-026-75681-7)

Mangiola S, Roth-Schulze AJ, Trussart M, et al. (2023). sccomp: robust
differential composition and variability analysis for single-cell data.
*Proc Natl Acad Sci USA* 120(33):e2203828120.
[doi:10.1073/pnas.2203828120](https://doi.org/10.1073/pnas.2203828120)

Assefa AT, Verbist B, Van den Berge K (2026). Assessing differential
cell composition in single-cell studies using voomCLR. *Bioinformatics*
42(1):btaf637.
[doi:10.1093/bioinformatics/btaf637](https://doi.org/10.1093/bioinformatics/btaf637)

## Examples

``` r
library(SingleCellExperiment)
library(S4Vectors)
set.seed(1)
donors <- paste0("D", seq_len(12))
donor_type_counts <- t(rmultinom(
    length(donors),
    size = 30,
    prob = c(A = 0.3, B = 0.4, C = 0.3)
))
rownames(donor_type_counts) <- donors
cell_donor <- rep(donors, each = 30)
cell_type <- unlist(
    lapply(seq_along(donors), function(i) {
        rep(colnames(donor_type_counts), donor_type_counts[i, ])
    }),
    use.names = FALSE
)
counts <- matrix(rpois(20 * length(cell_donor), 5), nrow = 20,
    dimnames = list(paste0("G", seq_len(20)),
        paste0("c", seq_along(cell_donor))))
sce <- SingleCellExperiment(
    assays = list(counts = counts),
    colData = DataFrame(
        cell_id = colnames(counts),
        donor_id = cell_donor,
        cell_type = cell_type))
exp_mat <- cbind(X1 = rep(0:1, 6), age = seq(40, 62, by = 2))
rownames(exp_mat) <- donors
scee <- build_scee(sce, exp_mat, sample_col = "donor_id")
run_exposure_composition(
    scee,
    exposure = "X1",
    covariates = "age",
    min_cells = 10L,
    support_min_donors = 3L,
    support_min_total_cells = 10L
)
#> DataFrame with 3 rows and 33 columns
#>      celltype coefficient log2_ratio_change        se  statistic        df
#>   <character>   <numeric>         <numeric> <numeric>  <numeric> <integer>
#> A           A  -0.0596652        -0.0596652  0.326187 -0.1829175         9
#> B           B   0.0365780         0.0365780  0.456544  0.0801193         9
#> C           C   0.0230873         0.0230873  0.478851  0.0482139         9
#>      pvalue pvalue_underflow_clamped      padj padj_holm   padj_bh    ci_low
#>   <numeric>                <logical> <numeric> <numeric> <numeric> <numeric>
#> A  0.858917                    FALSE         1         1  0.962599 -0.797550
#> B  0.937896                    FALSE         1         1  0.962599 -0.996197
#> C  0.962599                    FALSE         1         1  0.962599 -1.060148
#>     ci_high  ci_level    ci_multiplicity simultaneous_ci_low
#>   <numeric> <numeric>        <character>           <numeric>
#> A   0.67822      0.95 nominal_unadjusted            -1.01648
#> B   1.06935      0.95 nominal_unadjusted            -1.30261
#> C   1.10632      0.95 nominal_unadjusted            -1.38154
#>   simultaneous_ci_high simultaneous_ci_method ratio_change ratio_ci_low
#>              <numeric>            <character>    <numeric>    <numeric>
#> A             0.897146 Bonferroni_95_percen..     0.959487     0.575325
#> B             1.375770 Bonferroni_95_percen..     1.025678     0.501320
#> C             1.427712 Bonferroni_95_percen..     1.016132     0.479583
#>   ratio_ci_high ratio_simultaneous_ci_low ratio_simultaneous_ci_high
#>       <numeric>                 <numeric>                  <numeric>
#> A       1.60016                  0.494322                    1.86238
#> B       2.09849                  0.405391                    2.59506
#> C       2.15296                  0.383810                    2.69020
#>      exposure  n_donors n_replicates n_donors_present zero_replicates
#>   <character> <integer>    <integer>        <numeric>       <numeric>
#> A          X1        12           12               12               0
#> B          X1        12           12               12               0
#> C          X1        12           12               12               0
#>   total_cells median_observed_proportion support_n_donors support_total_cells
#>     <numeric>                  <numeric>        <integer>           <integer>
#> A         115                   0.333333               12                 115
#> B         137                   0.366667               12                 137
#> C         108                   0.266667               12                 108
#>                   method
#>              <character>
#> A donor_logratio_HC3_w..
#> B donor_logratio_HC3_w..
#> C donor_logratio_HC3_w..
```
