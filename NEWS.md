# exposomeSC 0.99.0

Initial Bioconductor submission.

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

* `plot_network_comparison()` produces side-by-side or overlay
  visualisations of cell-type-specific network differences.

* `plot_stability_surface()` displays the lambda-pi calibration
  heatmap for hyperparameter transparency.

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
  `sampleMap()`, `exposureNames()` with replacement methods.
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
