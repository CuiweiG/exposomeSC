# Impute exposure values below the limit of detection

Environmental exposure measurements often contain values reported only
as below the limit of detection (LOD), that is, left-censored values.
This function replaces them by substitution, by regression on order
statistics or by multiple imputation from a censored log-normal model.
Missing values are kept apart from censored ones and are never imputed.

## Usage

``` r
exposure_impute_lod(
  exposure_matrix,
  lod,
  method = c("lod_sqrt2", "lod_2", "ros", "multiple"),
  censored = NULL,
  M = 5L,
  verbose = TRUE
)
```

## Arguments

- exposure_matrix:

  Numeric matrix. Rows = samples, columns = exposure variables. In the
  exposures named in `lod`, an infinite value that is not censored is an
  error.

- lod:

  Named numeric vector of positive, finite detection limits, one per
  exposure. Names must match column names of `exposure_matrix`: a name
  that is not a column is ignored with a warning, and the function stops
  when no name is. Exposures not in `lod` are left unchanged.

- method:

  Character. Imputation method:

  `"lod_sqrt2"`

  :   Substitution with LOD / sqrt(2).

  `"lod_2"`

  :   Substitution with LOD / 2.

  `"ros"`

  :   Regression on order statistics (ROS) for one detection limit per
      exposure; see Details.

  `"multiple"`

  :   Multiple imputation from a log-normal model fitted by censored
      maximum likelihood to a bootstrap sample of the data, as proposed
      by Lubin et al. (2004); returns `M` completed matrices.

- censored:

  `NULL` (default) or a logical matrix with the same dimensions as
  `exposure_matrix`, `TRUE` where a value was reported as below the LOD.
  The numeric value at a censored position is ignored and may be `NA`.
  With `NULL`, a value is censored when it is not `NA` and is at or
  below the LOD. Either way, `NA` values that are not marked as censored
  are missing and stay `NA`.

- M:

  Integer. Number of imputations (only for `method = "multiple"`).
  Default 5.

- verbose:

  Logical. Report the numbers of censored and missing values per
  exposure. Default TRUE.

## Value

For `method != "multiple"`: the numeric matrix with censored values
replaced, of the same dimensions as the input. For
`method = "multiple"`: an
[`S4Vectors::SimpleList`](https://rdrr.io/pkg/S4Vectors/man/SimpleList-class.html)
of `M` such matrices. In both cases a `data.frame` with one row per
exposure of `lod` that is a column of `exposure_matrix` (`exposure`,
`n_detected`, `n_censored`, `n_missing`, `method`) records what was
done; it is the `"lod_summary"` attribute of the matrix, or
`metadata(result)$lod_summary` for the list.

## Details

**Substitution** is simple, but it distorts variances and associations,
increasingly so as the censored fraction grows (Lubin et al. 2004).

**ROS.** With a single detection limit every censored value lies below
every detected value, so the censored values take the lowest ranks among
the non-missing values of the exposure. Blom plotting positions \\(i -
3/8) / (n + 1/4)\\ are computed for all \\n\\ non-missing values, the
logarithms of the detected values are regressed on the normal quantiles
of their plotting positions, and each censored value is replaced by the
value of the fitted line at one of the lowest plotting positions, capped
just below the LOD (Helsel 2012). ROS reproduces the distribution of the
censored values, not their individual values: the imputed values are
assigned to the censored samples in row order, so they carry no
information on which censored sample had the larger true value. Use the
ROS values for summaries of the distribution; for association analyses
prefer `method = "multiple"`. ROS is not the reverse Kaplan-Meier
estimator, which is not implemented here.

**Multiple imputation.** For each of the `M` imputations the non-missing
values of an exposure are resampled with replacement, the mean and
standard deviation of the log values are estimated by maximum likelihood
with the censored values contributing \\P(X \le \mathrm{LOD})\\, and
every censored value is drawn from the fitted log-normal distribution
truncated to \\(0, \mathrm{LOD}\]\\. The bootstrap step carries the
uncertainty of the fitted distribution into the imputations (Lubin et
al. 2004). Analyse each completed matrix and combine the estimates with
Rubin's rules (Rubin 1987). The draws use the random-number stream of
the session, so call [`set.seed()`](https://rdrr.io/r/base/Random.html)
first for reproducible imputations. If the fit fails for an exposure
even on its observed values, an error names it; impute that exposure
with another method.

ROS and multiple imputation assume a log-normal distribution and need at
least three detected values; an exposure with fewer is imputed with LOD
/ sqrt(2), with a warning, and this is recorded in `lod_summary`. An
exposure with more than 80% of its non-missing values censored gives a
warning.

Other implementations: `rexposome::ilod()` and `rexposome::imputeLOD()`
impute values below the LOD with the QRILC or MinProb methods of the
imputeLCMD package, which are not the methods offered here.

## References

Lubin JH, Colt JS, Camann D, et al. (2004). Epidemiologic evaluation of
measurement data in the presence of detection limits. *Environ Health
Perspect* 112(17):1691–1696.
[doi:10.1289/ehp.7199](https://doi.org/10.1289/ehp.7199)

Helsel DR (2012). Statistics for Censored Environmental Data Using
Minitab and R. 2nd ed. Wiley.
[doi:10.1002/9781118162729](https://doi.org/10.1002/9781118162729)

Rubin DB (1987). Multiple Imputation for Nonresponse in Surveys. Wiley.
[doi:10.1002/9780470316696](https://doi.org/10.1002/9780470316696)

## Examples

``` r
set.seed(1)
mat <- matrix(rlnorm(200, 0, 1), nrow = 40, ncol = 5,
    dimnames = list(paste0("S", 1:40),
        c("Pb", "Cd", "Hg", "As", "BPA")))
lod <- c(Pb = 0.5, Cd = 0.3, Hg = 0.2, As = 0.1, BPA = 0.4)
## Values reported as below the LOD, and one genuinely missing value
censored <- sweep(mat, 2, lod[colnames(mat)], "<")
mat[censored] <- NA
mat["S1", "Pb"] <- NA
censored["S1", "Pb"] <- FALSE
imputed <- exposure_impute_lod(mat, lod, method = "ros",
    censored = censored)
#> Pb: 6 below LOD, 1 missing, 33 detected
#> Cd: 2 below LOD, 0 missing, 38 detected
#> Hg: 0 below LOD, 0 missing, 40 detected
#> As: 0 below LOD, 0 missing, 40 detected
#> BPA: 7 below LOD, 0 missing, 33 detected
attr(imputed, "lod_summary")
#>   exposure n_detected n_censored n_missing method
#> 1       Pb         33          6         1    ros
#> 2       Cd         38          2         0    ros
#> 3       Hg         40          0         0   none
#> 4       As         40          0         0   none
#> 5      BPA         33          7         0    ros
imputations <- exposure_impute_lod(mat, lod, method = "multiple",
    censored = censored, M = 3L, verbose = FALSE)
length(imputations)
#> [1] 3
```
