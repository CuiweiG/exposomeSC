# Power analysis for sc-ExWAS study design

Approximates the power to detect a donor-level exposure effect on one
gene in one cell type at the pseudobulk level, for a range of donor
numbers.

## Usage

``` r
estimate_power(
  n_donors = seq(10, 100, by = 10),
  effect_size = 0.5,
  n_genes = 5000L,
  n_celltypes = 5L,
  alpha = 0.05,
  dispersion = 0.1,
  base_mean = 100,
  exposure_sd = 1
)
```

## Arguments

- n_donors:

  Integer vector. Donor sample sizes to evaluate. Default:
  `seq(10, 100, by = 10)`.

- effect_size:

  Numeric. Expected log2 fold change per unit exposure. Default 0.5.

- n_genes:

  Integer. Number of genes tested (for multiple testing correction).
  Default 5000.

- n_celltypes:

  Integer. Number of cell types tested. Default 5.

- alpha:

  Numeric. Family-wise significance level; the per-test level is
  `alpha / (n_genes * n_celltypes)`. Default 0.05.

- dispersion:

  Numeric. Negative-binomial dispersion \\\phi\\ of donor pseudobulk
  counts, whose variance is \\\mu + \phi\mu^2\\. Default 0.1.

- base_mean:

  Numeric. Mean donor pseudobulk count \\\mu\\ of the gene. Default 100;
  `Inf` gives the large-count limit.

- exposure_sd:

  Numeric. Standard deviation of the exposure variable across donors.
  Default 1 (standardised). For a binary exposure that a proportion
  \\q\\ of the donors have, the standard deviation is \\\sqrt{q(1 -
  q)}\\, 0.5 for two groups of equal size; `effect_size` is then the
  log2 fold change between the groups.

## Value

A `DataFrame` with columns: n_donors, power, effect_size, exposure_sd,
n_tests, alpha_adjusted, dispersion, base_mean.

## Details

For a negative-binomial GLM with a log link, the large-sample standard
error of the exposure coefficient on the \\\log_2\\ scale is \$\$SE =
\frac{\sqrt{1/\mu + \phi}}{\log(2)\\ s_x \sqrt{n}},\$\$ where \\s_x\\ is
`exposure_sd` and \\n\\ the number of donors. Power is that of a
two-sided Wald test at the Bonferroni per-test level. The approximation
treats the dispersion as known and uses the normal distribution, so it
is optimistic for small numbers of donors; simulate from the planned
design when an accurate figure matters. Bonferroni is conservative
relative to Benjamini-Hochberg control of the false discovery rate.

Pseudobulk inference has donors, not cells, as its independent units.
More cells per donor raise \\\mu\\ and so shrink the \\1/\mu\\ term, but
they do not increase \\n\\.

## References

Squair JW et al. (2021). Confronting false discoveries in single-cell
differential expression. *Nat Commun* 12:5692.

## Examples

``` r
pw <- estimate_power(n_donors = seq(10, 50, by = 10),
    effect_size = 0.5, n_genes = 2000, n_celltypes = 3)
pw
#> DataFrame with 5 rows and 8 columns
#>    n_donors     power effect_size exposure_sd   n_tests alpha_adjusted
#>   <numeric> <numeric>   <numeric>   <numeric> <numeric>      <numeric>
#> 1        10  0.124663         0.5           1      6000    8.33333e-06
#> 2        20  0.585803         0.5           1      6000    8.33333e-06
#> 3        30  0.897429         0.5           1      6000    8.33333e-06
#> 4        40  0.984320         0.5           1      6000    8.33333e-06
#> 5        50  0.998319         0.5           1      6000    8.33333e-06
#>   dispersion base_mean
#>    <numeric> <numeric>
#> 1        0.1       100
#> 2        0.1       100
#> 3        0.1       100
#> 4        0.1       100
#> 5        0.1       100
```
