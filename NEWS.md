# exposomeSC 0.99.4

* The single-cell values of `incov_pbmc_subset` may now be used and
  redistributed under CC BY 4.0, with the permission of the authors who
  deposited E-MTAB-9357, so the whole data set is under that licence.
  `?incov_pbmc_subset` quotes their statement.

# exposomeSC 0.99.3

* Documentation only. NEWS is cut down to the changes that affect users,
  and `?incov_pbmc_subset` gives the sources and dates of the terms that
  apply to its single-cell data.

# exposomeSC 0.99.2

## Removed from this release

The first release keeps the donor-level association workflow and the
exploratory cross-sectional networks. Functions outside that core, or whose
estimators have not been validated, are removed rather than shipped with
warnings.

* `run_mediation()`, `run_causal_mediation()`, `run_network_mediation()`,
  `mediational_evalue()` and `mediational_evalue_ci()`: their identification
  assumptions and resampling properties have not been validated, and the
  package reports associations, not causal effects.
* `compute_iers()`, whose directness component needed a decomposition table
  that no function in the package produces.
* `run_state_coupling()` and `run_spatial_exwas()`, which lie outside the
  core donor-level analysis and have not been validated.
* `run_differential_network()`, which split a continuous exposure at its
  median and ignored covariates.
* `run_temporal_network()`, whose tests checked only its arguments, with
  `plot_temporal_dynamics()` and the `TemporalNetwork` class.
* `method = "fisher_z"` of `run_comparative_network()`: networks estimated
  from the same donors are not independent samples, and penalised partial
  correlations do not follow the distribution the test assumes. The
  comparison is now descriptive; the arguments `method` and `adjust`, the
  columns `p_value` and `p_adjusted` and the `method` element of
  `comparisonSummary()` are gone.
* The unexported developmental function `run_erd_interventional()`.
* The simulated files `demo_sce.rds` and `demo_exposures.csv`, with
  `?demo_data` and `create_demo_data.R`. Two real data sets replace them.

## Changes that can alter results

* `run_sc_exwas()` with the edgeR backend calls `glmQLFit(legacy = FALSE)`
  without a preceding `estimateDisp()`. edgeR 4.x releases before 4.10.0
  used those estimates in different ways, so a call could give different
  results with different edgeR versions; edgeR 4.10.0, now required,
  ignores them.
* The voom-dream backend with `filter_genes = FALSE` no longer tests genes
  without any count. Through limma's empirical Bayes prior, such genes
  changed the standard errors and p-values of the other genes in their cell
  type.
* `run_sc_exwas_glmm(family = "nbinom")` gives each cell type its own
  exposure slope also when `make.names()` maps several labels to one name;
  such cell types shared one slope. Before fitting, both families now leave
  out unlabelled cells, the cells of donors with a missing exposure or
  covariate, and any cell type left without cells. The linear mixed model
  had reported the next cell type's estimate under the label of the first
  when the first was left without cells, and with a single labelled cell
  type it had analysed the unlabelled cells as cells of that type.
* `run_multi_exwas()` with `covariates` and the default `exposures` no
  longer tests the covariates as exposures, which changes `padj_all` and,
  with the edgeR backend, the number of rows.
* An `exposure` given as a factor selects the exposure it names, not the one
  at the position of its integer code (`run_exposure_composition()` still
  refuses a factor).
* `exposure_impute_lod()` treats `NA` as missing, not as below the limit of
  detection (LOD); the new argument `censored` marks the censored values,
  by default those at or below the LOD.
* Its method `"kaplan_meier"` is now `"ros"`, the regression on order
  statistics it computed. Where the censored values were coded as `NA` and
  nothing else was missing, `censored = is.na(exposure_matrix)` gives the
  earlier values of `"lod_sqrt2"`, `"lod_2"` and `"ros"`, except with
  `"ros"` when the detected values are all equal.
* `method = "multiple"` now follows Lubin et al. (2004), a censored
  log-normal fit to a bootstrap sample, where 0.99.1 fitted the detected
  values alone. Its imputations differ from the earlier ones and, with few
  detected and many censored values, can lie orders of magnitude below the
  LOD. With fewer than three detected values it uses LOD/sqrt(2) and warns,
  where 0.99.1 drew uniformly between 0 and the LOD.
* The dose-response functions leave out donors with a missing exposure value
  before they count donors, and `run_dose_response_gam()` caps the default
  `k` at the number of distinct exposure values.
* `run_sc_mixture()` gives an exposure whose coefficient cannot be estimated
  the weight `NA`, not 0, and gives `NA` weights and effects, not equal
  weights and zero effects, when no coefficient can be estimated.
* `run_mixture_qgcomp(bootstrap = FALSE)` refers the mixture effect to a t
  distribution on the residual degrees of freedom of the fitted model: its
  intervals are wider, and with covariates its p-values are larger.
* `run_exposure_composition()` draws from the caller's random-number stream
  unless `seed` is given (the default was 1). Its permutations leave the
  residual of a donor that is alone in its `permutation_strata` stratum in
  place; such a donor often took the residual of a donor in another stratum.
* `run_cell_coupling()` leaves out donors with a missing moderator value,
  and gene-protein pairs whose exposure is constant among their donors or
  collinear with the moderators.
* The block graphical lasso of the network functions is now
  `glasso::glasso()` on the correlation matrix, with separate within- and
  between-layer penalties; it used to threshold a ridge-regularised inverse.
* The network functions stop when `transform = "vst"` cannot be computed,
  instead of using log-CPM values, and `run_exposure_network()` stops when
  coglasso is not installed, instead of using the block graphical lasso.
* `run_exposure_network()` with `covariates` selects features by partial
  rank correlation. Ranking the covariate-adjusted values had put tied
  values, such as the zeros of sparse genes, in the order of the covariates,
  so that such genes took on the correlation of the covariates with the
  exposure.
* Its stability selection, with or without covariates, counts only the
  candidates that vary among the donors with an exposure value, in the
  default `stability_q` and in the reported `pfer_bound`; constant
  candidates had made the default too large and the bound too small.
* Features that are constant across the donors are left out of a network,
  and `run_exposure_network()` with `covariates` no longer selects them. The
  block graphical lasso had kept such a feature as a node without edges and
  coglasso had failed on it; after adjustment for `exposure` or
  `covariates`, one with the same non-zero value in every donor could become
  rounding noise, scaled to unit variance, and gain edges.

## Changes to arguments and outputs

* A `SingleCellExposomeExperiment` stores the name of its donor column
  (`sampleColumn()`) instead of a `sampleMap` slot; objects saved with
  0.99.1 must be rebuilt with `build_scee()`. `cbind()` and `rbind()` have
  methods, and `exposureData<-` can add, drop or reorder exposures.
* The analysis functions take the object as `x`, no longer take
  `sample_col`, and default `celltype_col` to `"cell_type"`.
* New `donor_pseudobulk()` returns the pseudobulk counts that the models use.
* `run_sc_exwas(method = "dreamlet")` is now `"voom_dream"`, which is what it
  ran (it never used the dreamlet package), and the DESeq2 alias
  `"pseudobulk"` is gone. All backends return the same columns, sorted by
  `padj_global`; DESeq2 results report `logCPM` instead of `baseMean`.
* The `padj` column of `run_sc_exwas_glmm()` and `run_sc_exwas_pb_offset()`,
  their adjustment across all strata, is now `padj_global`, as in
  `run_sc_exwas()`.
* `run_sc_exwas_glmm()`, `run_sc_exwas_pb_offset()` and
  `run_cell_coupling()` return a DataFrame instead of a data.frame, and
  `run_sc_mixture()` and `run_mixture_qgcomp()` a DataFrame instead of a
  list, with the mixture weights in `metadata()`.
* `exposure_impute_lod(method = "multiple")` returns a SimpleList instead of
  a list; `run_meta_exwas()` no longer returns `meta_z`, an alias of
  `meta_statistic`; `estimate_power()` returns the `exposure_sd` it was
  given as a column.
* The results of `run_dose_response()` and `run_cell_coupling()` and the
  edge tables of `run_comparative_network()` have no row names.
* `run_gsea()` ranks by the `rank_by` column only.
* `plot_stability_surface()` is now `plot_edge_stability()`, which draws a
  reference line only at a `threshold` it is given, not at 0.6.
* The network functions take `transform` (`"vst"` or `"log_cpm"`) and
  `genes`, the genes to use, so that the networks of two cell types can be
  built on the same genes. `run_exposure_network()` takes `lambda_w` and
  `lambda_b` (defaults 0.3 and 0.5, the values it used before), records its
  penalties, donors and covariates, and records its method as
  `"exposure_restricted"` instead of `"exposure_driven"`. Networks wrapped
  with `precomputed_network` name their omic-layer column `omic_layer`, not
  `omic`.
* The network classes have constructors and accessors, such as
  `networkAdjacency()`, `networkNodes()` and `diffEdges()` (see
  `?networkPrecision`).
* Names that a model formula would misread (such as `"PCB-153"` or `.`),
  names that a function uses for variables of its own, an exposure that is
  also a covariate or moderator of its own analysis, and missing or infinite
  metabolite values in the network functions are now errors that name the
  problem. Such input used to give wrong or empty results or unrelated
  errors, mostly without a message.
* Other input that used to be dropped or misread without a message is now an
  error: several labels in `celltype`, whose cells were pooled; a variable
  to adjust for in the network functions, or an exposure of `as_scee()`,
  that is not in the data; and a confounder of `run_cell_coupling()` that is
  an unordered label with more than two values, which was ranked as if
  ordered.
* `run_multi_exwas()` stops when `run_sc_exwas()` fails for one of its
  exposures, instead of leaving that exposure out of the result and of
  `padj_all` with a warning.

## New data and vignettes

* `hlca_airway_subset`: raw counts of 352 airway epithelial cells from 12
  donors of the Human Lung Cell Atlas core, six heavy smokers and six never
  smokers (CC BY 4.0).
* `incov_pbmc_subset`: 1,032 blood mononuclear cells from 43 patients of the
  INCOV cohort (Su et al. 2020), with 192 antibody counts per cell, 13
  plasma compounds of exogenous origin and 182 plasma metabolites. The
  clinical and plasma tables are CC BY 4.0; for the single-cell values, from
  ArrayExpress E-MTAB-9357, the authors state no licence and the EMBL-EBI
  terms of use apply.
* `inst/scripts/` holds the scripts that build both data sets, and the HLCA
  preparation scripts, which now take their paths as arguments.
* Three vignettes, in BiocStyle and on the real data: `exposomeSC`
  (replacing `exposomeSC-quickstart` and `introduction`), `network-tutorial`
  and `permutation-and-power` (replacing `simulation-benchmark`).

## Requirements

* R (>= 4.6.0) and edgeR (>= 4.10.0). The coglasso estimator needs coglasso
  (>= 1.1.0), the block graphical lasso needs glasso, and the default
  `transform = "vst"` needs DESeq2.
* BiocGenerics and scrapper (>= 1.6.0) are new imports. mediation,
  SpatialExperiment, BiocCheck and statmod are no longer suggested.

# exposomeSC 0.99.1

* The slope-heterogeneity test of `run_state_coupling()` returned an `NA`
  p-value on the macOS builds of the Bioconductor staging checks. The fit
  now tries metafor's optimisers in turn and keeps only a converged fit; if
  the F test is still unusable, a likelihood-ratio test is used, and the
  `method` field of the `slope_heterogeneity` attribute names the test.

# exposomeSC 0.99.0

Initial Bioconductor submission.

## Container and association analysis

* `SingleCellExposomeExperiment` extends `SingleCellExperiment` with
  sample-level exposure data, exposure metadata and the cell-to-donor
  mapping; `build_scee()` builds one and `as_scee()` converts an ExposomeSet.
* The accessors are `cellSampleMap()` and `exposureVariables()`, not
  `sampleMap()` and `exposureNames()`, which MultiAssayExperiment and
  rexposome export as S4 generics; whichever package was attached last
  would have masked the other's.
* `run_sc_exwas()` tests exposures within cell types on donor pseudobulk
  counts, by default with robust edgeR quasi-likelihood models;
  `run_multi_exwas()` runs several exposures.
* `run_sc_exwas_glmm()` fits mixed models as an experimental sensitivity
  analysis; its calibration depends on the design and must be checked.
* `run_dose_response()`, `run_dose_response_gam()` and
  `run_interaction_test()` for dose-response curves and exposure-by-cell-type
  interaction; `run_sc_mixture()` and `run_mixture_qgcomp()` for mixtures.
* `run_exposure_composition()` analyses the relative composition of captured
  cells; it does not estimate absolute tissue abundance.
* `run_meta_exwas()` for meta-analysis across cohorts, `run_spatial_exwas()`
  for region-resolved ExWAS via SpatialExperiment, and the experimental
  `run_mediation()` and `run_causal_mediation()`.
* `exposure_impute_lod()` imputes values below the limit of detection,
  `seurat_to_exposure()` extracts donor metadata, `run_gsea()` tests gene
  sets with fgsea, and `estimate_power()` computes power.
* edgeR results report the coefficient standard error as `NA`:
  `glmQLFTest()` gives a quasi-likelihood F test but no standard error. The
  signed square root of F is for ranking only.
* `run_meta_exwas()` rejects edgeR results, since inverse-variance synthesis
  needs valid standard errors.
* Functions that take a `seed` restore the caller's random-number generator,
  including its kind, when they return; seeding goes through withr.

## EXPERIMENTAL: Cross-omic network estimation

* `run_celltype_network()` estimates cross-omic networks within a cell type.
  Its edges are model-dependent summaries, not validated
  conditional-independence discoveries.
* `run_exposure_network()` restricts the network to features associated with
  an exposure, by stability selection over a marginal Spearman screen or by
  univariate tests; `network_method` chooses the estimator.
* `run_comparative_network()` classifies edges as shared, unique or
  differential and reports Jaccard overlap; its `fisher_z` mode is
  experimental. `run_temporal_network()` records edge presence over time.
* `run_network_mediation()` is experimental. Plots, a simulator
  (`simulate_crossomic_network()`) and the classes `CelltypeNetworkResult`,
  `NetworkComparison` and `TemporalNetwork` complete the module.

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
