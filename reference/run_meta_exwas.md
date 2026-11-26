# Meta-analysis across cohorts for sc-ExWAS

Combines independent-cohort sc-ExWAS estimates using inverse-variance
weighted fixed-effect synthesis or random-effects synthesis. Keys are
matched exactly on gene, cell type, and exposure.

## Usage

``` r
run_meta_exwas(
  results_list,
  cohort_names = NULL,
  method = c("fixed", "random"),
  effect_col = NULL,
  se_col = NULL
)
```

## Arguments

- results_list:

  A list of `DataFrame`s or data frames. Each cohort must contain unique
  `gene`-`celltype`-`exposure` keys and an effect estimate plus standard
  error.

- cohort_names:

  Character vector with one unique name per cohort. Defaults to
  `Cohort_1`, `Cohort_2`, and so forth.

- method:

  Character. `"fixed"` (default), or `"random"` for REML between-cohort
  variance with Hartung–Knapp inference. Random-effects synthesis
  requires metafor and at least three cohorts per key.

- effect_col:

  Optional character scalar naming the effect column. When `NULL`, the
  first available column among `log2FC`, `log2FoldChange`, and
  `estimate` is used separately in each cohort.

- se_col:

  Optional character scalar naming the standard-error column. When
  `NULL`, the first available column among `se` and `lfcSE` is used
  separately in each cohort. edgeR QL results from `run_sc_exwas` are
  deliberately unsupported because `glmQLFTest` does not supply
  coefficient standard errors.

## Value

A `DataFrame` with one row per exact key and columns `meta_effect`,
`meta_se`, `meta_statistic`, `meta_df`, `meta_pvalue`, `meta_padj`,
`I2`, `Q_pvalue`, `tau2`, `n_cohorts`, `cohorts`, `direction`, and
`method`. Under Hartung–Knapp inference `meta_statistic` is a t
statistic with `meta_df` degrees of freedom; under fixed-effect
synthesis it is a z statistic and `meta_df` is `Inf`.

## Details

Fixed-effect synthesis uses \$\$\hat{\beta} = \frac{\sum_i
w_i\hat{\beta}\_i}{\sum_i w_i}, \quad w_i = 1/SE_i^2.\$\$ Random-effects
synthesis uses restricted maximum likelihood for \\\tau^2\\ and
Hartung–Knapp small-sample inference. It is deliberately unavailable for
only two cohorts. The \\I^2\\ and Cochran Q statistics are reported as
heterogeneity diagnostics; they do not establish transportability or
independence of the contributing cohorts.

The function never maps gene identifiers or cell-type labels.
Harmonisation must occur before calling it. Within-cohort duplicate keys
are an error, preventing accidental double weighting after label
aggregation. Inverse-variance synthesis also rejects edgeR QL results
from [`run_sc_exwas`](run_sc_exwas.md), recognised by their `method`
label or their metadata, because the QL F test supplies no coefficient
standard errors. Use a backend that estimates valid coefficient standard
errors for each independent cohort before meta-analysis.

## References

Borenstein M et al. (2009). *Introduction to Meta-Analysis*. Wiley.

Viechtbauer W (2010). Conducting meta-analyses in R with the metafor
package. *Journal of Statistical Software* 36:1–48.
[doi:10.18637/jss.v036.i03](https://doi.org/10.18637/jss.v036.i03)

## Examples

``` r
library(S4Vectors)
r1 <- DataFrame(
    gene = c("A", "B"), celltype = "Mono",
    log2FC = c(0.5, -0.3), se = c(0.1, 0.15),
    exposure = "PM2.5"
)
r2 <- DataFrame(
    gene = c("A", "B", "C"), celltype = "Mono",
    log2FC = c(0.4, -0.1, 0.2), se = c(0.12, 0.2, 0.1),
    exposure = "PM2.5"
)
meta <- run_meta_exwas(
    list(r1, r2),
    cohort_names = c("CohortA", "CohortB")
)
meta
#> DataFrame with 2 rows and 16 columns
#>          gene    celltype meta_effect   meta_se meta_statistic   meta_df
#>   <character> <character>   <numeric> <numeric>      <numeric> <numeric>
#> 1           A        Mono    0.459016 0.0768221        5.97505       Inf
#> 2           B        Mono   -0.228000 0.1200000       -1.90000       Inf
#>   meta_pvalue        I2  Q_pvalue      tau2 n_cohorts    exposure
#>     <numeric> <numeric> <numeric> <numeric> <integer> <character>
#> 1 2.30013e-09         0  0.522053         0         2       PM2.5
#> 2 5.74331e-02         0  0.423711         0         2       PM2.5
#>           cohorts   direction                 method   meta_padj
#>       <character> <character>            <character>   <numeric>
#> 1 CohortA;CohortB          ++ inverse_variance_fixed 4.60027e-09
#> 2 CohortA;CohortB          -- inverse_variance_fixed 5.74331e-02
```
