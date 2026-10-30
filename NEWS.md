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

## Interface

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

## Documentation

* The demo data are described as what they are: independent Poisson counts
  with an added PM2.5 signal, and simulated exposures without units.
* References added or corrected: edgeR v4, the graphical lasso, the
  published coglasso paper, dream and variancePartition, and a DOI for Lee
  and Han (2024).
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
