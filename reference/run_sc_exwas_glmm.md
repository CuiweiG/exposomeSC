# Cell-level linear mixed model sc-ExWAS with finite-df inference

Fits, per gene, the cell-level model
`y ~ exposure * celltype + covariates + (1 | donor) + (1 | donor:celltype)`
and reports cell-type-specific exposure effects. For
`family = "gaussian"` the model is a linear mixed model fitted by
[`lmerTest::lmer`](https://rdrr.io/pkg/lmerTest/man/lmer.html) with
**REML = TRUE** (mandatory), and p-values are referred to a finite
denominator degrees of freedom – **Satterthwaite** (default) or
**Kenward-Roger** – via
[`lmerTest::contest1D`](https://rdrr.io/pkg/lmerTest/man/contest.html).
A Wald test referred to the normal distribution would be
anticonservative with few donors.

## Usage

``` r
run_sc_exwas_glmm(
  x,
  exposure,
  celltype_col = "cell_type",
  covariates = NULL,
  target_genes = NULL,
  ddf = c("Satterthwaite", "Kenward-Roger"),
  family = c("gaussian", "nbinom"),
  BPPARAM = BiocParallel::SerialParam()
)
```

## Arguments

- x:

  A
  [`SingleCellExposomeExperiment`](SingleCellExposomeExperiment-class.md).
  Each cell's donor is read from the `colData` column named by
  [`sampleColumn`](SCEE-accessors.md)`(x)`.

- exposure:

  Character; a donor-level column of `exposureData(x)`.

- celltype_col:

  Character. Column in `colData` containing cell type annotations.
  Default `"cell_type"`.

- covariates:

  Character vector; donor covariate columns in `exposureData` (numeric),
  added as fixed-effect main effects. The names `exposure`, `donor`,
  `celltype` and `y` are used by the model and cannot be covariate
  names; with `family = "nbinom"` the same holds for `logoff` and the
  `expo_<cell type>` columns. Each name must be a syntactic R name (see
  [`make.names`](https://rdrr.io/r/base/make.names.html)) other than
  `.`, `...` and names such as `..1`, which R reads specially in a
  formula, and the exposure itself cannot be a covariate.

- target_genes:

  Character; genes to test (default: all rows of `x`). A cell-level LMM
  per gene is costly – passing a focused gene set is strongly
  recommended.

- ddf:

  Denominator degrees-of-freedom method for `family = "gaussian"`:
  `"Satterthwaite"` (default) or `"Kenward-Roger"` (needs pbkrtest; use
  for \\n \lesssim 30\\ donors).

- family:

  `"gaussian"` (LMM on cell-level \\\log_2\\ CP10k, default) or
  `"nbinom"` (negative-binomial GLMM on cell counts).

- BPPARAM:

  A `BiocParallel` param used to parallelise over genes. Defaults to
  [`BiocParallel::SerialParam()`](https://rdrr.io/pkg/BiocParallel/man/SerialParam-class.html),
  so no worker process is started unless one is asked for; pass, for
  example, `BiocParallel::SnowParam(workers = 2L)` to fit genes in
  parallel.

## Value

A `DataFrame`, one row per gene x cell type: `gene`, `celltype`,
`log2FC`, `se`, `statistic` (Satterthwaite /KR \\t\\ for gaussian, LR
\\\chi^2\\ for nbinom), `df` (finite denominator df; 1 for the nbinom
LRT), `pvalue`, `padj_global` (BH across all gene-by-cell-type strata,
as in [`run_sc_exwas`](run_sc_exwas.md); there is no within-cell-type
`padj`), and `method`. Rows ordered by `padj_global`, then `pvalue`. A
gene-by-cell-type stratum that cannot be fitted is left out of the
result; with `family = "nbinom"`, a stratum whose likelihood-ratio test
cannot be computed is kept with `NA` `pvalue` and `padj_global`. Its
`statistic` is `NA` when the reduced model or
[`anova()`](https://rdrr.io/r/stats/anova.html) fails, and close to 0
when the test has no degrees of freedom, as when a covariate is a linear
function of the exposure. `padj_global` adjusts over the strata with a
p-value, and the function stops when no stratum can be fitted.

## Details

**Why REML.** Maximum likelihood biases the variance components
downward, which shrinks the residual scale and inflates the exposure
\\t\\-statistic; the Kenward-Roger adjusted covariance and the
Satterthwaite degrees-of-freedom approximation are both defined only on
the REML fit. The function therefore forces `REML = TRUE` and does not
expose a switch.

**Random effects.** The exposure varies between donors, so the
replicates for a cell-type-specific exposure effect are
donor-by-cell-type units, not cells. The model therefore includes a
donor-by-cell-type random intercept as well as a donor random intercept;
without it, variation shared by the cells of one donor and cell type is
treated as independent and the cell-type-specific tests are
anticonservative. With a single cell type the model reduces to
`y ~ exposure + covariates + (1 | donor)`. A variance component
estimated at zero (a singular fit) is expected when that source of
variation is absent and is not reported.

**Cell-type-specific effects.** Under treatment contrasts the `exposure`
coefficient is the effect at the reference cell type and each
`exposure:celltype` coefficient is a difference from it. Each cell
type's effect is therefore a linear combination of fixed effects, tested
with `contest1D` so that the degrees of freedom and standard error
account for the covariance between the main effect and the interaction.

**family = "nbinom".** Satterthwaite and Kenward-Roger are linear
-mixed-model constructions and do not extend to a generalised LMM, so
the negative-binomial path fits
[`lme4::glmer.nb`](https://rdrr.io/pkg/lme4/man/glmer.nb.html) on an
equivalent cell-type-nested reparametrisation of `exposure * celltype`
(one exposure slope column per cell type) and tests each cell type's
effect by a likelihood-ratio test (`anova` of the full versus the
column-dropped model) – again not a Wald-normal test. This path is
markedly more expensive (one `glmer.nb` refit per cell type per gene);
for count data prefer the donor-level
[`run_sc_exwas_pb_offset`](run_sc_exwas_pb_offset.md), which is
calibration- equivalent and far faster (Lee & Han 2024). The `ddf`
argument is ignored for `family = "nbinom"`.

**Missing values.** The cells of a donor whose exposure or a covariate
is missing, and cells without a cell-type label, are left out, with a
message, before either family is fitted, so that the full and reduced
models of a likelihood-ratio test use the same cells.

## Note

Calibrate against a donor-permutation null, with the exposure permuted
across donors, before quoting any type-I rate; Satterthwaite/KR control
the finite-sample error only under the working model, so the empirical
check is the headline.

## References

Kuznetsova A, Brockhoff PB & Christensen RHB (2017). lmerTest package:
tests in linear mixed effects models. *J Stat Softw* 82(13).

Kenward MG & Roger JH (1997). Small sample inference for fixed effects
from restricted maximum likelihood. *Biometrics* 53:983-997.

Halekoh U & Hojsgaard S (2014). A Kenward-Roger approximation and
parametric bootstrap methods for tests in linear mixed models – the R
package pbkrtest. *J Stat Softw* 59(9).

## See also

[`run_sc_exwas_pb_offset`](run_sc_exwas_pb_offset.md) for the
donor-level default.

## Examples

``` r
donor_ids <- paste0("D", seq_len(10L))
exposure <- rep(c(0, 1), each = 5L)
cell_donor <- rep(donor_ids, each = 12L)
cell_type <- rep(
    rep(c("TypeA", "TypeB"), each = 6L),
    times = length(donor_ids)
)
cell_ids <- paste0("cell", seq_along(cell_donor))
counts <- outer(
    seq_len(2L),
    seq_along(cell_ids),
    function(gene, cell) 2L + ((7L * gene + 5L * cell) %% 11L)
)
exposed_type_a <- exposure[match(cell_donor, donor_ids)] == 1 &
    cell_type == "TypeA"
counts[1L, exposed_type_a] <- counts[1L, exposed_type_a] + 3L
storage.mode(counts) <- "integer"
dimnames(counts) <- list(c("G1", "G2"), cell_ids)
sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = counts),
    colData = S4Vectors::DataFrame(
        cell_id = cell_ids,
        donor_id = cell_donor,
        cell_type = cell_type
    )
)
donor_design <- matrix(exposure, ncol = 1L,
    dimnames = list(donor_ids, "exposure"))
scee <- build_scee(sce, donor_design, sample_col = "donor_id")
if (requireNamespace("lme4", quietly = TRUE) &&
        requireNamespace("lmerTest", quietly = TRUE)) {
    mixed_result <- run_sc_exwas_glmm(
        scee,
        exposure = "exposure",
        celltype_col = "cell_type",
        target_genes = "G1",
        family = "gaussian",
        BPPARAM = BiocParallel::SerialParam()
    )
    mixed_result
}
#> DataFrame with 2 rows and 9 columns
#>          gene    celltype    log2FC        se statistic        df    pvalue
#>   <character> <character> <numeric> <numeric> <numeric> <numeric> <numeric>
#> 1          G1       TypeB 0.3088714  0.170915  1.807164       116 0.0733286
#> 2          G1       TypeA 0.0863299  0.170915  0.505104       116 0.6144439
#>   padj_global                 method
#>     <numeric>            <character>
#> 1    0.146657 lmerTest/Satterthwaite
#> 2    0.614444 lmerTest/Satterthwaite
```
