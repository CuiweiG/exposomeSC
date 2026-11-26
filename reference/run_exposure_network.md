# Network among exposure-associated features

Pre-selects the transcripts and metabolites associated with an exposure,
by stability selection or univariate testing, and estimates a cross-omic
network among the selected features only. Cheng et al. (2024) analysed
transcriptomic and metabolomic signatures of air pollution in two such
stages: features associated with the exposures were identified first,
and conditional multi-omic networks estimated afterwards.

## Usage

``` r
run_exposure_network(
  x,
  metabolites,
  celltype,
  exposure,
  celltype_col = "cell_type",
  covariates = NULL,
  selection_method = c("stability", "univariate"),
  fdr_threshold = 0.05,
  pi_threshold = 0.6,
  stability_q = NULL,
  min_cells = 10L,
  network_method = c("coglasso", "block_glasso"),
  lambda_w = 0.3,
  lambda_b = 0.5,
  transform = c("vst", "log_cpm"),
  genes = NULL
)
```

## Arguments

- x:

  [`SingleCellExposomeExperiment`](SingleCellExposomeExperiment-class.md).
  Each cell's donor is read from the `colData` column named by
  [`sampleColumn`](SCEE-accessors.md)`(x)`.

- metabolites:

  Numeric matrix, or a data frame of numeric columns; rows = donors,
  cols = metabolites. Row names must be donor IDs. A matrix without
  columns, or a missing or infinite value for an analysed donor, is an
  error, and so is a column with the name of another column or of a
  candidate gene: two nodes cannot have the same name. Leave the
  exposure itself out of the columns: a column that holds the values of
  the exposure, or a monotone transform of them, has a rank correlation
  of 1 or -1 with it, is always selected and becomes a node.

- celltype:

  Character; cell type to analyse.

- exposure:

  Character; exposure variable name, a column of `exposureData(x)`;
  another name is an error.

- celltype_col:

  Character; cell type column. Default `"cell_type"`.

- covariates:

  Character vector or NULL; columns of `exposureData(x)`. The selection
  is made given the covariates (see Details), and they are regressed out
  of the transcript and metabolite values before the network is
  estimated; another name is an error, and so is the exposure itself.
  Donors with a missing covariate value are left out, with a message. A
  covariate enters the selection by its ranks and the adjustment of the
  values as the number it is, so either way it is taken as ordered; give
  a covariate with more than two unordered categories as one 0/1 column
  for each category but one.

- selection_method:

  Character; `"stability"` (default) or `"univariate"`.

- fdr_threshold:

  Numeric; FDR threshold for univariate selection. Default 0.05.

- pi_threshold:

  Numeric in (0.5, 1\]; minimum proportion of subsamples in which a
  feature must be selected under `selection_method = "stability"`.
  Default 0.6.

- stability_q:

  Integer or `NULL`; number of features selected in each subsample under
  `selection_method = "stability"`. The default `NULL` uses \\\max(1,
  \lfloor\sqrt{(2\pi\_{thr} - 1)p}\rfloor)\\, the largest value for
  which the Meinshausen-Buhlmann bound on the expected number of falsely
  selected features is at most 1 (when \\(2\pi\_{thr} - 1)p \ge 1\\).
  Here \\p\\ is the number of candidate features that vary across the
  donors with a value of the exposure.

- min_cells:

  Integer; minimum cells per pseudobulk.

- network_method:

  Character; `"coglasso"` (default) or `"block_glasso"`, the estimator
  for the network on the selected features. `"coglasso"` needs at least
  one selected transcript and metabolite; when only one omic layer is
  selected the block graphical lasso is used instead, with a warning.
  Otherwise it requires the coglasso package, version 1.1.0 or later (an
  error when it is not installed), and chooses the penalties from a 10
  by 10 grid with the defaults of
  [`coglasso::xstars()`](https://drquestion.github.io/coglasso/reference/xstars.html).
  The estimator used is recorded in `metadata(result)$network_method`.

- lambda_w:

  Numeric; graphical-lasso penalty for pairs of features in the same
  omic layer, on the correlation scale, for the block graphical lasso
  (`network_method = "block_glasso"`, or the fallback when only one omic
  layer is selected). Default 0.3.

- lambda_b:

  Numeric; graphical-lasso penalty for pairs of features in different
  omic layers, on the correlation scale, for the block graphical lasso.
  Default 0.5.

- transform:

  Character; the transform of the pseudobulk counts, `"vst"` (default)
  or `"log_cpm"`, as in
  [`run_celltype_network`](run_celltype_network.md). The call stops when
  `"vst"` cannot be computed; it does not switch to the other transform.

- genes:

  Character vector or `NULL`; row names of `x`, the genes that are
  candidates for the selection, taken after the transform as in
  [`run_celltype_network`](run_celltype_network.md). Fewer candidates
  leave fewer tests to adjust for and a smaller `stability_q`. A name
  that is not a row name of `x`, or that is the name of more than one
  row, is an error. Default `NULL` (all genes).

## Value

A [`CelltypeNetworkResult`](CelltypeNetworkResult-class.md) restricted
to exposure-associated features. Its metadata record the expression
transform, the covariates, the donors of the network (`donors`) and the
number of them with a value of the exposure (`n_donors_exposure`), on
whom the selection was made, the selection settings, for stability
selection `pfer_bound`, in `constant_features` the candidates that are
constant across the donors with a value of the exposure and so could not
be selected, and in `selection` the penalties of the network: `lambda_w`
and `lambda_b` for the block graphical lasso, or the size of the grid
and the `selected_lambda_w`, `selected_lambda_b` and `selected_c` that
coglasso selected, with `stars_variability`.

## Details

Stage 1: Feature selection, among the donors with a value of the
exposure. Each transcript and metabolite is screened by its Spearman
correlation with the exposure. With `covariates` it is the partial rank
correlation: the ranks of the covariates are regressed out of the ranks
of the feature and of the exposure, and the residuals are correlated. A
feature without variation across these donors, such as a gene without a
count in the cell type, is neither tested nor ranked by either method,
and \\p\\ below counts the features that vary. `"univariate"` tests the
correlation (asymptotic *t* approximation, with one degree of freedom
fewer for each covariate that is not constant or a combination of the
others among these donors) with BH correction over these features. The
approximation is poor when both the exposure and a feature take one
value in all but a few donors, as an exposure that only a few donors
have, or that is above its limit of detection in a few, and a gene with
a count in one or two donors: the p-value of such a gene is far too
small whenever its donors are among those few, and with many such genes
among the candidates some pass the BH correction together. For an
exposure of this kind, name the candidate genes with `genes` and leave
out those with a count in few donors. `"stability"` applies stability
selection (Meinshausen and Buhlmann 2010) to a marginal screen: in each
of 100 random subsamples of half the donors with an observed exposure,
the `stability_q` features with the largest absolute correlation with
the exposure are selected, and features selected in a proportion of at
least `pi_threshold` of subsamples are retained. The screen remains
defined when features outnumber donors. Under the exchangeability and
better-than-random-guessing assumptions of Meinshausen and Buhlmann, the
expected number of falsely selected features is at most \\q^2 /
((2\pi\_{thr} - 1)p)\\, reported as `pfer_bound`. Candidates with a
count in very few donors are hardly ever among the `stability_q` highest
of a subsample, which the exchangeability assumption does not allow for;
with many of them among the candidates the bound is too small. `genes`
can leave them out. At least three features must be selected.

Transcript values are log-CPM or VST values and therefore relative to
each donor's library size. When the exposure strongly changes genes that
carry a large share of the counts, the remaining genes shift in the
opposite direction and can pass the screen; with few measured genes this
can dominate the univariate screen.

Stage 2: Network estimation on the selected features only, with
collaborative graphical lasso model selection by XStARS (`"coglasso"`)
or the block graphical lasso (`"block_glasso"`). It is one network for
all the donors analysed. The exposure is not regressed out of the
features and is not a node, so an edge can reflect the association that
two features share with the exposure, and the network does not estimate
how the dependence between features changes with the exposure.

## References

Cheng SL et al. (2024). Multiomic signatures of traffic-related air
pollution in London reveal potential short-term perturbations in gut
microbiome-related pathways. *Environ Sci Technol* 58:8771-8782.
[doi:10.1021/acs.est.3c09148](https://doi.org/10.1021/acs.est.3c09148)

Meinshausen N, Buhlmann P (2010). Stability selection. *J R Stat Soc
Series B*, 72(4):417-473.

## Examples

``` r
set.seed(3)
donors <- sprintf("D%02d", 1:20)
donor <- rep(donors, each = 40)
exposure <- stats::setNames(stats::rnorm(20), donors)
counts <- matrix(stats::rpois(20 * length(donor), 20), nrow = 20,
    dimnames = list(paste0("G", 1:20), paste0("c", seq_along(donor))))
## two genes rise and two fall with the exposure
effect <- c(0.8, 0.8, -0.8, -0.8)
for (g in 1:4)
    counts[g, ] <- stats::rpois(length(donor),
        20 * exp(effect[g] * exposure[donor]))
sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = counts),
    colData = S4Vectors::DataFrame(donor_id = donor, cell_type = "Mono"))
scee <- build_scee(sce, matrix(exposure, ncol = 1,
    dimnames = list(donors, "PM2.5")), sample_col = "donor_id")
metabolites <- matrix(stats::rnorm(20 * 3), nrow = 20,
    dimnames = list(donors, paste0("M", 1:3)))
metabolites[, "M1"] <- metabolites[, "M1"] + 1.5 * exposure
if (requireNamespace("glasso", quietly = TRUE)) {
    net <- run_exposure_network(scee, metabolites, celltype = "Mono",
        exposure = "PM2.5",
        selection_method = "univariate",
        network_method = "block_glasso",
        transform = "log_cpm")
    net
}
#> [exposomeSC] Exposure network: 5/23 features selected
#> CelltypeNetworkResult for 'Mono'
#>   Method: exposure_restricted
#>   Nodes: 5 (4 transcripts, 1 metabolites)
#>   Edges: 10 (4 cross-omic)
#>   Donors: 20
```
