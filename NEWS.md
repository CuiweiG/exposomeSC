# exposomeSC 0.99.2

## Removed from this release

The first release keeps the donor-level association workflow and the
exploratory cross-sectional networks. Functions that sit outside that core,
or whose estimators have not been validated, are removed rather than shipped
with warnings.

* `run_mediation()`, `run_causal_mediation()`, `run_network_mediation()`,
  `mediational_evalue()` and `mediational_evalue_ci()`. Their identification
  assumptions and resampling properties have not been validated, and the
  package reports associations, not causal effects.
* `compute_iers()`. Its directness component needed a decomposition table
  that no function in the package produces.
* `run_state_coupling()` and `run_spatial_exwas()`, which lie outside the
  core donor-level analysis and have not been validated.
* `run_differential_network()`, which split a continuous exposure at its
  median and ignored covariates, and `run_temporal_network()`, whose only
  tests were argument checks, with `plot_temporal_dynamics()` and the
  `TemporalNetwork` class.
* The unexported developmental function `run_erd_interventional()`.
* mediation and SpatialExperiment are no longer suggested.

## Changes that can alter results

* `run_sc_exwas(method = "edgeR")` now calls `glmQLFit(legacy = FALSE)`
  without a preceding `estimateDisp()`. Earlier edgeR 4.x releases used
  those estimates in different ways, so one call could give different
  results on different edgeR versions; from edgeR 4.10.0 they are ignored.
  edgeR 4.10.0 is now required, and the edgeR version is recorded in
  `metadata(result)$parameters`.
* `method = "block_glasso"` in `run_celltype_network()` and
  `run_exposure_network()` is now the graphical lasso (`glasso::glasso()`)
  on the correlation matrix, with separate within- and between-layer
  penalties. It used to threshold a ridge-regularised inverse, which is not
  the graphical lasso. It needs the glasso package.
* `exposure_impute_lod()` treats `NA` as missing, not as below the LOD;
  values reported as below the LOD are marked with the new `censored`
  argument (or, without it, are the values at or below the LOD). Missing
  values stay `NA`.
* `exposure_impute_lod(method = "kaplan_meier")` is now
  `method = "ros"`: it was always regression on order statistics, not the
  reverse Kaplan-Meier estimator. Where the censored values were coded as
  `NA` and nothing else was missing, passing `censored = is.na(x)` gives
  the values of the earlier version.
* `exposure_impute_lod(method = "multiple")` now follows Lubin et al.
  (2004): for each imputation it resamples the data, fits a log-normal
  distribution by censored maximum likelihood and draws the censored values
  from that distribution truncated at the LOD. It used to fit the detected
  values only, ignoring the censoring, and did not carry the uncertainty of
  the fit into the imputations. The fit starts from two points and keeps the
  better fit: both use a spread on the log scale equal to the largest of the
  spread of the detected values, half their distance to the LOD and 0.05;
  the first is centred on the detected values, the second puts the censored
  share of the distribution below the LOD. It retries with a second
  optimiser when the first does not converge and skips a bootstrap sample
  that has fewer than two distinct detected values or whose fit fails; when
  none of 100 samples can be used, it draws from the fit to the observed
  values, and it stops with a message naming the exposure only when that fit
  fails too. Imputed values are drawn on the log-probability scale, so that
  they stay inside (0, LOD]. With few detected values and many censored
  ones, they can lie orders of magnitude below the LOD: for detected values
  of 40, 50, 60 and 70 and 100 censored values below an LOD of 1, their
  median is of the order of 1e-7, where the earlier version imputed
  LOD/sqrt(2).
* `run_exposure_network(network_method = "coglasso")` stops with an
  installation hint when coglasso is not installed, instead of switching to
  the block graphical lasso.
* The `residual_df` diagnostic of the edgeR backend of `run_sc_exwas()`
  gives the adjusted residual degrees of freedom of the quasi-likelihood fit
  instead of the nominal ones (see Bug fixes).
* Networks from `run_celltype_network()`, and the features that
  `run_exposure_network()` selects with `covariates`, can change for data
  with a feature that is constant across the donors, such as a metabolite
  with one value in every donor: such a feature is now left out.
  `method = "block_glasso"` kept it as a node without edges; with `exposure`
  or `covariates`, a feature with the same non-zero value in every donor
  could instead stay in as rounding noise scaled to unit variance (see Bug
  fixes).
* `run_comparative_network(method = "fisher_z")` with more than two networks
  gives the p-value 1 instead of `NA` to an edge whose tested pairs all have
  the same partial correlation; such an edge now counts in `p_adjusted`,
  whose values for the other edges can therefore be larger (see Bug fixes).
* The permutation p-value of `run_exposure_composition()` can change when a
  stratum of `permutation_strata` holds a single donor in a model row other
  than the first (see Bug fixes).
* `run_sc_exwas_glmm(family = "nbinom")` gives each cell type its own
  exposure slope also when `make.names()` maps several cell-type labels to
  one name; such cell types shared one slope (see Bug fixes).
* The dose-response functions give other results when a donor has no dose:
  `n_donors` and the genes of the default screen change;
  `run_dose_response()` chooses the degree of the polynomial from the donors
  with a dose, so with five of them it fits a quadratic instead of a cubic
  and `nonlinear_model`, `p_nonlinear`, `best_model` and `AIC_best` change;
  and in `run_dose_response_gam()` the leave-one-out R² and, where the
  default `k` changes, the columns of the smooth change (see Bug fixes).
* `run_cell_coupling()` gives other `n_donors`, `mean_r` and `mean_n_cells`
  when a donor has no moderator value, and leaves out a pair that is then
  left with fewer than `min_donors` donors (see Bug fixes).
* The linear mixed model of `run_sc_exwas_glmm()` gives other estimates with
  one labelled cell type and unlabelled cells, and, when every cell of the
  first cell type is left out, no longer reports the next type's estimate
  under the first type's label (see Bug fixes).

## Interface

* `exposureData<-` can add, drop or reorder exposures: `exposureInfo` now
  follows the new column names, keeping the rows of retained exposures and
  adding rows of missing metadata for new ones. The two setters
  used to validate against each other, so an exposure could be replaced but
  not added or removed through the API.
* New `donor_pseudobulk()` returns the donor-by-cell-type pseudobulk counts
  as one SummarizedExperiment per cell type, with the donors' exposure data
  as column data.
* `run_sc_exwas()` returns the same columns, in the same order, from every
  backend: DESeq2 results report `logCPM` instead of `baseMean`, and
  `padj_global` now follows `padj`. DESeq2 and voom-dream results are now
  ordered by `padj_global` like edgeR results and, like them, carry
  `metadata(result)$parameters` (including the backend, its version,
  `se_method` and `statistic_type`) and `metadata(result)$diagnostics`.
  `schema_version` is now `"exposomeSC_sc_exwas_v3"` for all backends.
* `run_multi_exwas()` keeps each exposure's result metadata in
  `metadata(result)$per_exposure`.
* `run_sc_exwas_glmm()`, `run_sc_exwas_pb_offset()` and
  `run_cell_coupling()` return a DataFrame instead of a data.frame.
  `run_sc_mixture()` returns its per-gene table as a DataFrame with the
  weights and summaries in `metadata()`. `run_mixture_qgcomp()` returns a
  one-row DataFrame with the mixture effect, its standard error, confidence
  interval and p-value, with the weights and the qgcomp fit in `metadata()`.
  `exposure_impute_lod(method = "multiple")` returns a SimpleList of
  completed matrices, and every method records what it did per exposure
  (`lod_summary`).

* `run_sc_exwas(method = "dreamlet")` is now `method = "voom_dream"`, which is
  what it runs; it never used the dreamlet package. The `"pseudobulk"` alias
  for DESeq2 is removed.
* `plot_stability_surface()` is renamed `plot_edge_stability()`. The fixed
  dashed line at 0.6 is gone; pass `threshold` to draw one.
* `run_exposure_composition()` defaults to `seed = NULL` and then draws from
  the caller's random-number stream, so `set.seed()` before the call works
  as expected. An integer seed still works as before.
* `run_gsea()` ranks by the column named in `rank_by` and stops if that
  column is missing or not numeric. It no longer falls back between
  `"statistic"` and `"stat"`, or `"log2FC"` and `"log2FoldChange"`.
* `run_meta_exwas()` no longer returns `meta_z`, which duplicated
  `meta_statistic`.
* A SingleCellExposomeExperiment no longer keeps a second copy of the
  cell-to-donor mapping in a `sampleMap` slot. It stores the name of the
  donor column of `colData`, returned by the new `sampleColumn()`, and
  `cellSampleMap()` is derived from that column. `cbind()` now combines
  objects (it used to fail), renaming cells no longer makes the object
  invalid, and replacing `colData`, one of its columns or the cell names
  checks validity, so a donor without exposure data is an error at once.
  Objects saved with 0.99.1 need to be rebuilt with `build_scee()`.
* The analysis functions read each cell's donor from the stored column and
  no longer take `sample_col`. Every function that takes the container
  calls it `x` (formerly `scee` in `run_cell_coupling()`,
  `run_sc_exwas_glmm()`, `run_sc_exwas_pb_offset()`,
  `run_celltype_network()` and `run_exposure_network()`), and
  `celltype_col` defaults to `"cell_type"` throughout (formerly
  `"celltype"` in `run_cell_coupling()` and no default in `run_sc_exwas()`,
  `run_sc_exwas_glmm()` and `run_sc_exwas_pb_offset()`).
* CelltypeNetworkResult and NetworkComparison have constructors,
  `CelltypeNetworkResult()` and `NetworkComparison()`, and accessors:
  `networkPrecision()`, `networkAdjacency()`, `networkStability()`,
  `networkNodes()`, `networkCelltype()`, `networkMethod()` and `metadata()`
  for networks, and `diffEdges()`, `sharedEdges()`,
  `celltypeSpecificEdges()` and `comparisonSummary()` for comparisons.
  NetworkComparison now has a validity method, and CelltypeNetworkResult
  checks its stability matrix, node table, cell type and method. The node
  table must have one row per node, and where the matrices have row or
  column names they must equal its `feature` column, so a node table in
  another order is an error rather than a source of mislabelled edges. The
  adjacency matrix must be symmetric and hold only 0 and 1.

## Bug fixes

* `plot_edge_stability()` matched legend labels to colours by position, so a
  plot with a single edge type was labelled wrongly. Labels are now matched
  by name.
* `show()` for a SingleCellExposomeExperiment printed part of its summary
  with `message()`; all of it now goes to standard output.
* A dimension mismatch in a CelltypeNetworkResult is reported as one
  validity message instead of two fragments.
* `run_celltype_network(precomputed_network = )` named the omic layer column
  of its node table `omic`, so `show()` counted no transcripts or metabolites
  and `plot_celltype_network()` and `run_comparative_network()` failed on a
  wrapped network. The column is now `omic_layer`, with an integer `block`
  index, as for estimated networks.
* `run_sc_exwas()` refuses a covariate called `exposure`. That name is the
  model's exposure term, and a covariate of the same name silently replaced
  the exposure in the design, so the reported test was of the wrong
  variable.
* The edgeR backend of `run_sc_exwas()` warns when it drops a covariate
  column that is collinear with the other model terms from the design of a
  cell type it fits. The dropped columns were recorded in
  `metadata(result)$diagnostics` but the result otherwise looked adjusted
  for them.
* The `residual_df` diagnostic of the edgeR backend reads the adjusted
  residual degrees of freedom that `glmQLFit(legacy = FALSE)` returns
  (`df.residual.adj`). It looked for the legacy field `df.residual.zeros`
  and so reported the nominal residual degrees of freedom.
* `run_exposure_network(network_method = "coglasso")` warns when only one
  omic layer passes the selection stage and the block graphical lasso is
  used instead; the switch used to be silent.
* `plot_celltype_network(edge_alpha_by = "stability")` stops when the
  network has no stability scores, as `color_by = "stability"` already did;
  it used to fall back to the edge weights without saying so.
* `run_comparative_network()` accepts `metadata(network)$n_donors` stored
  as a double, which is how a number is normally typed; it used to fail
  with "values must be type 'integer'". A count that is not a single
  non-negative whole number is an error that names the network; `NULL` or
  `NA` mean that the count is unknown. With more than two networks, an
  edge whose tested pairs all have the same partial correlation gets the
  p-value 1, as with two networks, instead of `NA`; it now counts in
  `p_adjusted`, whose values for the other edges can therefore be larger.
* A CelltypeNetworkResult also checks that its precision matrix is
  numeric, that the adjacency matrix is numeric with a zero diagonal, that
  stability scores lie between 0 and 1, and that the node table has the
  `omic_layer` column that `show()`, the plots and the comparison read. The
  constructor and `run_celltype_network(precomputed_network = )` store a
  logical adjacency as 0/1 integers, and the wrapper zeroes the diagonal.
* `run_exposure_composition()` works when exactly two cell types pass the
  support filter. With one ILR coordinate `anova()` has no Pillai test and
  the call stopped; the exact F test of the added terms is used, and
  Pillai's trace is reported as df1 F / (df1 F + df2).
* `run_celltype_network()` leaves features that are constant across the
  donors, such as genes with no counts in the cell type, out of the network
  with a message, records them in `metadata(result)$constant_features`, and
  stops when every transcript or every metabolite is constant. A constant
  metabolite, or a constant gene when `top_var_genes` was not set, made
  coglasso stop with "'from' must be a finite number", while
  `method = "block_glasso"` kept it as a node without edges. With `exposure`
  or `covariates`, a feature with the same non-zero value in every donor,
  such as a metabolite measured at one level, was instead regressed into
  rounding noise, which scaling could turn into a feature of unit variance;
  the call then succeeded and the feature stayed in the network, where it
  could get edges. The check now runs before `exposure` and `covariates` are
  regressed out, and for the same reason `run_exposure_network()` with
  `covariates` no longer tests or selects a constant feature. A metabolite
  matrix without columns is an error in both functions; it used to fail with
  unrelated errors or, in `run_exposure_network()`, sometimes gave a network
  of transcripts only.
* `run_celltype_network()` stops when `exposure` or `covariates`, and
  `run_exposure_network()` when `covariates`, names a variable that is not
  in `exposureData(x)`; both check the names before any computation. Such
  names were dropped without a word, so a misspelt covariate left the
  network unadjusted, and `run_celltype_network()` still listed it in
  `metadata(result)$covariates`.
* `as_scee()` stops when `exposures` names an exposure that the ExposomeSet
  does not have. Unknown names used to be dropped silently, so a misspelt
  exposure was left out of the object without notice.
* `exposure_impute_lod()` stops, whatever the method, when an exposure in
  `lod` has an infinite value that is not censored. `"lod_sqrt2"` and
  `"lod_2"` returned it unchanged, as did `"ros"` and `"multiple"` for an
  exposure with fewer than three detected values or with neither `NA` nor
  censored values; otherwise `"ros"` and `"multiple"` stopped with "missing
  value where TRUE/FALSE needed". The value at a censored position is still
  ignored.
* `run_sc_exwas_glmm()` and `run_sc_exwas_pb_offset()` refuse a covariate
  whose name their model formulas use for a variable of their own
  (`exposure`, `donor`, `celltype`, `y` for the mixed model, plus `logoff`
  and the `expo_<cell type>` columns of its negative binomial family; `y`,
  `A`, `logoff` for the offset model). Nothing named the clash: `exposure`
  made the linear mixed model test the covariate instead of the exposure,
  and the negative binomial one report the covariate's slopes without a
  p-value; `donor` replaced the donor grouping factor, so donors with equal
  values were pooled; `celltype` replaced the cell-type factor, so the
  linear mixed model reported only the first cell type and the negative
  binomial one stopped with "No gene/cell-type stratum yielded an estimable
  model"; `y` (with a `model.matrix()` warning), `A` and an
  `expo_<cell type>` column were ignored, the last also leaving that cell
  type's likelihood-ratio test without a p-value; and `logoff` turned the
  library-size offset into a free coefficient. The same rule for
  `run_sc_exwas()` is now stated in its argument documentation.
* `run_celltype_network(precomputed_network = )` checks the labels and
  dimensions of the network it wraps: the adjacency must be square, the
  precision matrix must have the dimensions of the adjacency, the feature
  names must be unique, non-missing and non-empty, `block_sizes` must hold
  non-negative whole numbers named by omic layer, and `block_sizes` or
  `feature_blocks` must label every feature, with no missing or empty
  labels. A mismatch used to surface as a base R error from `rep()`,
  `dimnames<-` or `DataFrame()`, and a duplicated name or a missing block
  label was accepted silently. A 0 x 0 adjacency without feature names,
  which failed with "length of 'dimnames' [1] not equal to array extent", is
  now wrapped into an empty network.
* `plot_celltype_network(edge_alpha_by = "stability")` maps the edge
  opacity to the stability scores also when `highlight_cross_omic = FALSE`;
  without the highlighting the opacity followed the edge weights whatever
  `edge_alpha_by` said. A network in which some edges have no stability
  score is refused, instead of being drawn with those edges fully opaque.
* The source tarball keeps the vignette index `build/vignette.rds` that
  `R CMD build` writes; the allow-lists in `.Rbuildignore` removed it. A
  package installed from the tarball on a system where knitr could not be
  loaded therefore had no vignette index: `vignette(package = "exposomeSC")`
  listed nothing and `browseVignettes("exposomeSC")` found no vignettes.
* The Freedman-Lane permutation test of `run_exposure_composition()` moves
  residuals only within their `permutation_strata` stratum and leaves the
  residual of a donor that is alone in its stratum in place. For such a
  donor in model row k > 1, `sample()` read k as `1:k`, so in about
  (k - 1) / k of the permutations the donor took the residual of a donor of
  another stratum, that residual appeared twice, and the permutation
  p-value was computed partly from invalid permutations. When every
  stratum has at least two donors, the permutations, and so seeded
  results, are unchanged.
* `run_celltype_network()` leaves donors with a missing `exposure` or
  covariate value, and `run_exposure_network()` donors with a missing
  covariate value, out of the network, with a message, as `run_sc_exwas()`
  leaves them out of its design; when fewer than five donors remain, the
  error says so. One such donor used to stop both functions with the base
  R error "values must be length ..., but FUN(X[[1]]) result is length
  ...".
* Names that go into a model formula must be syntactic R names: the
  covariates of `run_sc_exwas()`, `run_sc_exwas_pb_offset()` and
  `run_sc_exwas_glmm()`, the exposures and covariates of `run_sc_mixture()`
  and `run_mixture_qgcomp()`, and the `meta_moderators` of
  `run_cell_coupling()`. A name such as `"PCB-153"` or `"age-years"` is
  refused with a message that suggests its `make.names()` form. These names
  were pasted into the formula and misread (`age-years` as `age` minus
  `years`), so the functions returned empty or all-`NA` results, mostly
  without a warning, or stopped with an error about the data or about an
  object such as `age`; a ratio such as `"MEHP/DEHP"` next to `"MEHP"` and
  `"DEHP"` silently gave a different model. A column named after an
  expression of other analysed columns, such as `"log(age)"` next to
  `"age"`, used to be computed from those columns; it is refused as well
  and must be renamed. The exposure of `run_sc_exwas()` is renamed
  internally and can still have any name.
* `run_sc_exwas_glmm(family = "nbinom")` keeps the exposure slopes of cell
  types apart whose labels `make.names()` maps to one name, such as
  `"CD16+ Mono"` and `"CD16- Mono"`. Their slope columns shared that name,
  so the model held one slope, for the last of those cell types, and
  reported it with its likelihood-ratio test for each of them: a cell type
  without an effect was given the effect of another.
* `run_cell_coupling(meta_moderators = )` accepts a one-sided formula such
  as `~ age + sex`, as documented; a formula used to fail with "unique()
  applies only to vectors". The formula may only list columns of
  `exposureData(x)`: terms such as `I(age^2)` or `age:sex` are refused
  rather than reduced to their variables.
* `run_dose_response()` and `run_dose_response_gam()` leave donors with a
  missing exposure value out, with a message, before they count donors.
  Such donors were counted in the donor checks, in the default gene screen
  and in `n_donors`, although the fits dropped them: `n_donors` was too
  large; with `max_degree = 2`, four donors of which one had no dose fitted
  a saturated quadratic (`AIC_best = -Inf`) instead of stopping; the
  leave-one-out R² of `run_dose_response_gam()` was `NA` for every gene;
  and the default `k` of its smooth counted them too, so that with two
  such donors it could exceed the donors with a dose, every smooth failed
  and the function returned an empty table without a message.
* `run_cell_coupling()` leaves donors with a missing moderator value out of
  the meta-regression, as it left out donors with a missing exposure, and
  says how many donors it leaves out. `metafor::rma()` dropped such donors
  with a warning, but `n_donors`, `mean_r`, `mean_n_cells` and the
  `min_donors` check still counted them, so a pair could be reported from
  fewer donors than `min_donors`. When fewer than `min_donors` donors with
  enough cells have values of the exposure and of every moderator, the
  function warns and returns an empty table; with missing exposures it
  used to return that table with only the message "0 pairs tested".
* `run_sc_exwas_glmm(family = "nbinom")` fits data with a single cell type.
  Its formulas began with `0 + celltype`, which `model.matrix()` refuses for
  a factor with one level, so every gene failed and the function stopped
  with "No gene/cell-type stratum yielded an estimable model". A single cell
  type now has a plain intercept.
* `run_sc_exwas_glmm()` leaves out, with a message, the cells of donors with
  a missing exposure or covariate and the cells without a cell-type label,
  before it fits either family. The fits with several cell types already
  dropped such cells, so their results are unchanged unless every cell of a
  cell type is left out. Such a type is now left out of the model: when it
  came first, the linear mixed model reported the next type's estimate under
  its label, and when only one type kept cells, both families stopped. The
  negative binomial model for a single cell type needs it: its reduced model
  has no exposure or cell-type column and would keep cells that the full
  model drops, which leaves every likelihood-ratio test without a p-value.
  The linear mixed model with one labelled cell type no longer analyses
  unlabelled cells as cells of that type.
* When DESeq2 is not installed, or its variance-stabilising transformation
  fails, the network functions say so before using log-CPM values; the
  transform used was already recorded in the result's metadata.

## Internal changes

* Pseudobulk counts are summed with `scrapper::aggregateAcrossCells()`
  instead of package code; the counts are identical. scrapper is a new
  import. (`scuttle::aggregateAcrossCells()` is deprecated in favour of the
  scrapper function.)
* The three `run_sc_exwas()` backends share their input checks and build
  their results with one internal constructor. Results are unchanged apart
  from the interface changes above.

## Documentation

* The demo data are described as what they are: independent Poisson counts
  with an added PM2.5 signal, and simulated exposures without units.
* References added or corrected: edgeR v4, the graphical lasso, the
  published coglasso paper, dream and variancePartition, and a DOI for Lee
  and Han (2024).
* `?run_sc_exwas` states which steps are delegated to other packages and
  which choices exposomeSC makes itself. `?run_exposure_composition` compares
  the method with propeller, crumblr, sccomp and voomCLR, with references.
  `?exposure_impute_lod` describes each method correctly and points to the
  LOD imputation in rexposome.
* Requires R 4.6.0. BiocCheck and statmod are no longer suggested.

# exposomeSC 0.99.1

* `run_state_coupling()` returned an `NA` p-value for its slope-heterogeneity
  test on both macOS builds of the Bioconductor staging checks, while every
  Linux and Windows build passed. The multilevel fit stops at variance
  components of zero for the example data, and whether the default optimiser
  stops there cleanly depends on the platform. The fit now tries metafor's
  optimisers in turn and keeps only a converged fit with finite coefficients.
  If the F-test's coefficient covariance is still unusable, the same question
  is put as a likelihood-ratio test between the full and the additive model,
  both fitted by maximum likelihood, and the `method` field of the
  `slope_heterogeneity` attribute names the test that was used.

# exposomeSC 0.99.0

Initial Bioconductor submission.

## Corrections made before submission

* `run_celltype_network(method = "coglasso")` selected its penalties twice.
  It built the path with `coglasso::bs()`, which already runs its own
  selection, then ran XStARS again on the result and discarded the first
  selection. It now builds the path with
  `coglasso::coglasso()` and selects once. `run_exposure_network()` had the
  same pattern and is fixed the same way.
* `subsample_ratio` never reached XStARS: it was converted into a repetition
  count, `ceiling(1 / (1 - subsample_ratio))`, so the default 0.8 meant five
  subsamples of coglasso's default size. It is now passed as
  `stars_subsample_ratio`, and the number of subsamples is the new argument
  `rep_num` (default 20, as in coglasso).
* `stability_scores` held XStARS's single variability figure for the whole
  selected model, as a one-by-one matrix, where the class documents
  edge-level selection probabilities. It now holds the edge-level selection
  frequencies at the selected penalties; the variability figure, the selected
  penalties and the subsampling settings are kept in `metadata$selection`.
  `plot_stability_surface()` had been drawing an empty plot for every
  coglasso network as a result. The documentation of `stability` now says
  what it does: it controls whether those frequencies are returned, not
  whether XStARS runs, and `block_glasso` ignores it.
* `plot_celltype_network()` drew edges that could not be seen, with the
  absolute partial precision, typically a few hundredths, used directly as
  opacity. It labelled the edge-type legend by position, so a network with
  only cross-omic edges was labelled "Within-omic"; and it described any
  network without stored edge counts as having no cross-omic edges. Opacity
  is now scaled to a visible range, labels are matched by name, and the
  subtitle counts the edges drawn.
* The accessors `sampleMap()` and `exposureNames()` are renamed
  `cellSampleMap()` and `exposureVariables()`. Both names were already
  exported as S4 generics by current Bioconductor packages,
  MultiAssayExperiment and rexposome, so whichever package was attached
  last masked the other and the masked generic failed on the other's
  objects. `as_scee()` bridges from rexposome, so its own workflow loads
  both. The slot keeps its name, and objects built earlier remain valid.
* `run_sc_exwas()` and `run_multi_exwas()` stop on arguments they do not use.
  A misspelled argument such as `covariats = "age"` previously ran an
  unadjusted analysis without any message.
* `run_sc_exwas(method = "dreamlet")` is described accurately: a
  one-row-per-donor design has no random effects, so `dream()` fits it with
  limma and the moderated t statistic uses residual degrees of freedom. Cell
  types skipped by the DESeq2 and voom-dream backends are now reported as
  warnings rather than messages.
* `estimate_power()` computed the standard error on the natural-log scale for
  a log2 effect and omitted the Poisson term, which overstated power (0.996
  where simulation gave 0.56 for 20 donors, a log2 effect of 0.5, dispersion
  0.1, a mean count of 100 and 6,000 tests). It
  now uses `sqrt(1/base_mean + dispersion) / (log(2) * exposure_sd * sqrt(n))`
  and gains `base_mean`.
* `run_interaction_test()` ignored `covariates` and treated cell types from
  the same donor as independent. It now compares exposure slopes between cell
  types with donor fixed effects and covariate-by-cell-type terms, and reports
  `df_interaction`.
* `run_dose_response()` tested the AIC-selected polynomial against the linear
  model, which overstates the evidence for nonlinearity; `p_nonlinear` is now
  the pre-specified test of the highest-degree polynomial, named in
  `nonlinear_model`.
* `run_dose_response_gam()` compared the linear and GAM fits with a deviance
  F-test that rejected a truly linear response in 68-75% of simulated data
  sets at the 5% level. `p_nonlinear` now tests a fully penalised smooth added
  to a linear term (Wood 2013), which rejected 6-7% in the same simulation.
* `run_mixture_qgcomp()` returned the intercept p-value as `mixture_pvalue`;
  it now returns the p-value of the mixture effect. The documentation no longer
  claims weights or weight intervals from `qgcomp.boot()`, which estimates
  neither.
* `run_sc_mixture()` computes weights as documented (mean absolute per-gene
  coefficients), tolerates tied exposure values, drops donors with missing
  data, and stops on unknown covariates.
* `run_mediation()` is experimental and warns once per session, consistent with
  `run_causal_mediation()`; it and the mixture functions stop on unknown
  covariates instead of dropping them.
* `run_sc_exwas_glmm()` adds a donor-by-cell-type random intercept. Without it,
  cell-type-specific exposure tests rejected up to 13% of null simulations at
  the 5% level.
* `run_state_coupling()` tests slope heterogeneity with one multilevel
  meta-regression that accounts for donors shared between bins, replacing a
  Cochran Q test that treated bins as independent, and handles tied states.
* `run_spatial_exwas()` no longer pools k-means regions across tissue sections
  with unregistered coordinates, normalises on the whole transcriptome rather
  than the reported genes, and seeds k-means.
* `as_scee()` reads the `exp` assay of an `ExposomeSet`; previously it could
  not extract exposures from one. It keeps the rows for the samples present in
  the single-cell data, so a larger exposome cohort can be used directly.
* `simulate_crossomic_network()` ran only with an explicit seed, generated
  transcripts and metabolites independently of the cross-omic edges it
  returned, and ignored `composition_confounding` and `n_cells_per_donor`. The
  simulated data now follow the returned precision matrices, and composition
  sets pseudobulk depth.
* `compute_iers()` no longer accepts the unimplemented `weights = "adaptive"`
  and ranks unrounded scores; `run_gsea()` runs `fgsea` serially so that
  `set.seed()` makes its p-values reproducible; `run_cell_coupling()` drops a
  confounder that is constant within a donor instead of adding random noise to
  it; `seurat_to_exposure()` warns when values vary within a donor; and
  `run_temporal_network()` labels edges present at the first and last time
  points but not in between as `intermittent`.
* edgeR normalisation uses `normLibSizes()`, the current name of
  `calcNormFactors()`.
* Internal helpers for a repeated-sample sensitivity analysis (sample
  pseudobulk with donor-blocked voom) and for null-calibration runs are
  removed. No exported function used them.

## Random number state

* Functions that take a `seed` argument now restore the caller's random
  number generator when they return, including its kind. Previously
  `simulate_crossomic_network()`, `run_erd_interventional()` and the
  empirical-null permutation helpers left the global seed changed, and two
  of them left the generator set to L'Ecuyer-CMRG. Results for a given
  `seed` are unchanged.
* Seeding now goes through `withr`, which is imported.

## CORRECTED: edgeR quasi-likelihood uncertainty

* `run_sc_exwas(method = "edgeR")` now reports coefficient standard errors as
  numeric `NA`, because `edgeR::glmQLFTest()` supplies a quasi-likelihood F
  test but no coefficient standard error. The earlier `abs(logFC / sqrt(F))`
  quantity was not a valid standard error and is no longer produced.
* edgeR results record `schema_version = "exposomeSC_sc_exwas_edger_v2"`,
  `se_method = "not_available_edgeR_QL"`, and
  `statistic_type = "signed_sqrt_qlf"`. The signed square-root QL statistic is
  retained for directional ranking only; it is not a Wald statistic or a
  basis for coefficient confidence intervals.
* `run_meta_exwas()` rejects edgeR QL results, including older finite SE
  proxies, because inverse-variance synthesis requires valid coefficient
  standard errors. DESeq2 Wald and voom-dream uncertainty semantics are
  unchanged.

## EXPERIMENTAL: Cross-omic network estimation

* `run_celltype_network()` estimates exploratory cross-omic precision and
  adjacency matrices per cell type. Its edges are model-dependent summaries,
  not validated conditional-independence discoveries.

* `run_exposure_network()` replaces the mislabelled `stability_lasso` option
  with `selection_method = "stability"`: stability selection over a marginal
  Spearman screen in half-samples, which stays defined when features outnumber
  donors and reports the Meinshausen-Buhlmann bound on expected false
  selections. The former option refitted ordinary least squares on all features
  and failed whenever features outnumbered 80% of donors. The univariate screen
  no longer fails on features without variation, and network metadata record
  the expression transform (VST or log-CPM).

* `run_comparative_network()` classifies edges as shared, unique or
  differential and reports Jaccard overlap, without p-values by default. Its
  `fisher_z` mode is experimental and warns once per session: networks from the
  same donors are not independent samples, and penalised partial correlations
  do not follow the distribution the test assumes. Partial correlations are
  now derived from the precision matrices before the test.

* `run_temporal_network()` records descriptive edge presence across separately
  estimated time-point networks. It does not test longitudinal rewiring.

* `run_network_mediation()` is experimental and warns once per session: its
  selection, identification, resampling and confidence-interval properties
  have not been validated. Each cross-omic edge is tested as a single-mediator
  path with a two-sided percentile bootstrap p-value.

* `run_exposure_network()` gains `network_method` to choose between the
  collaborative and the block graphical lasso for the second stage.

## NEW: Network visualisation

* `plot_celltype_network()` renders cross-omic networks with
  omic-layer colouring, stability-based edge opacity, and
  bipartite layout options.

* `plot_network_comparison()` plots how much each edge's partial
  correlation varies across the compared cell types, for the edge set
  chosen with `highlight`, which the function previously ignored.

* `plot_stability_surface()` plots the selection probability of the most
  stable edges. It is not a lambda-by-pi surface, because the result
  object stores stability at the selected penalty only.

* `plot_celltype_network(color_by = "stability")` colours nodes by the
  largest selection probability among their edges instead of by omic
  layer, and the graph theme no longer asks for a font that most
  systems lack.

* `plot_temporal_dynamics()` tracks edge presence across time
  points with optional node highlighting.

## NEW: Simulation engine

* `simulate_crossomic_network()` generates realistic multi-omics
  data with known network structure, cell-type heterogeneity,
  and composition confounding for benchmarking.

## NEW: S4 classes

* `CelltypeNetworkResult`: stores precision matrix, adjacency,
  stability scores, and node metadata per cell type.

* `NetworkComparison`: stores differential, shared, and
  cell-type-specific edges with statistical tests.

* `TemporalNetwork`: stores time-resolved networks with edge
  dynamics tracking.

## Core infrastructure

* `SingleCellExposomeExperiment` S4 class extending
  `SingleCellExperiment` with sample-level exposure data,
  exposure metadata, and cell-to-donor mapping.
* `build_scee()` constructs integrated containers.
* Accessors: `exposureData()`, `exposureInfo()`,
  `cellSampleMap()`, `exposureVariables()` with replacement methods.
* `[` subsetting preserves exposure data.

## Cell-type-specific ExWAS

* `run_sc_exwas()` provides donor-pseudobulk exposure-association workflows;
  robust edgeR quasi-likelihood is the implemented default.
* `run_sc_exwas_glmm()` provides experimental mixed-model sensitivity paths;
  calibration is design-specific and must be assessed empirically.
* `run_multi_exwas()` batches requested exposures. Fixed-family, fail-closed
  handling is required before any cross-exposure multiplicity claim.

## Exposure-Response Decomposition (ERD)

* `run_decomposed_exwas()`, `erd_evalue()` and `erd_resolution_sensitivity()`
  are removed. The coefficient contrast they were built on does not identify a
  direct, indirect or cell-abundance-mediated effect, and an expression log2
  fold change does not determine the risk ratio an E-value needs.
* `run_erd_interventional()` is no longer exported. Its outcome is the
  within-cell-type expression rate, so its indirect effect is at most a
  contextual pathway rather than an abundance contribution, and its estimator,
  resampling and identification have not been validated. A function in that
  state should not sit in the public interface. It stays in the package, with
  its tests, for development.

## Dose-response and interaction

* `run_dose_response()` polynomial dose-response with AIC.
* `run_dose_response_gam()` GAM via mgcv with LOOCV R².
* `run_interaction_test()` formal exposure × celltype
  interaction F-test (Gelman & Stern 2006).

## Mixture and composition

* `run_sc_mixture()` quantile-scored mixture screening.
* `run_mixture_qgcomp()` formal qgcomp with bootstrap CI.
* `run_exposure_composition()` analyses relative captured-cell composition;
  it does not estimate absolute tissue abundance.

## Causal inference and meta-analysis

* `run_mediation()` is excluded from confirmatory scope pending a defensible
  identification and pre-selection contract.
* `run_causal_mediation()` is experimental and warns once per session. With
  `sensitivity = TRUE` it now returns `rho_at_zero` from
  `mediation::medsens()`, which was previously dropped, and a `method` column
  records whether `mediation::mediate()` or the difference-method bootstrap
  fallback produced each row.
* `run_meta_exwas()` cross-cohort meta-analysis with
  fixed/random effects, I², and direction strings.
* `run_spatial_exwas()` region-resolved ExWAS via
  SpatialExperiment.

## Preprocessing and interoperability

* `exposure_impute_lod()` below-LOD imputation (LOD/sqrt(2),
  ROS, multiple imputation).
* `as_scee()` converts rexposome ExposomeSet to SCEE.
* `seurat_to_exposure()` extracts donor-level metadata.
* `run_gsea()` pathway enrichment via fgsea.
* `estimate_power()` power/sample size calculator.
