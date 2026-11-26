# Package index

## Container

The class, its construction and its accessors

- [`` `[`( ``*`<SingleCellExposomeExperiment>`*`,`*`<ANY>`*`,`*`<ANY>`*`,`*`<ANY>`*`)`](SingleCellExposomeExperiment-class.md)
  [`cbind(`*`<SingleCellExposomeExperiment>`*`)`](SingleCellExposomeExperiment-class.md)
  [`rbind(`*`<SingleCellExposomeExperiment>`*`)`](SingleCellExposomeExperiment-class.md)
  [`` `colData<-`( ``*`<SingleCellExposomeExperiment>`*`,`*`<DataFrame>`*`)`](SingleCellExposomeExperiment-class.md)
  [`` `colData<-`( ``*`<SingleCellExposomeExperiment>`*`,`*`<NULL>`*`)`](SingleCellExposomeExperiment-class.md)
  [`` `dimnames<-`( ``*`<SingleCellExposomeExperiment>`*`,`*`<list>`*`)`](SingleCellExposomeExperiment-class.md)
  [`` `dimnames<-`( ``*`<SingleCellExposomeExperiment>`*`,`*`<NULL>`*`)`](SingleCellExposomeExperiment-class.md)
  [`show(`*`<SingleCellExposomeExperiment>`*`)`](SingleCellExposomeExperiment-class.md)
  : SingleCellExposomeExperiment: Container for single-cell exposome
  data
- [`build_scee()`](build_scee.md) : Build a SingleCellExposomeExperiment
- [`as_scee()`](as_scee.md) : Convert rexposome ExposomeSet to
  exposomeSC
- [`seurat_to_exposure()`](seurat_to_exposure.md) : Convert Seurat
  metadata to exposure matrix
- [`exposureData()`](SCEE-accessors.md)
  [`exposureInfo()`](SCEE-accessors.md)
  [`cellSampleMap()`](SCEE-accessors.md)
  [`sampleColumn()`](SCEE-accessors.md)
  [`` `sampleColumn<-`() ``](SCEE-accessors.md)
  [`exposureVariables()`](SCEE-accessors.md)
  [`` `exposureData<-`() ``](SCEE-accessors.md)
  [`` `exposureInfo<-`() ``](SCEE-accessors.md) : Accessors for
  SingleCellExposomeExperiment objects
- [`donor_pseudobulk()`](donor_pseudobulk.md) : Donor pseudobulk counts
  for each cell type

## Example data

- [`hlca_airway_subset`](hlca_airway_subset.md) : Airway epithelial
  cells of smokers and never smokers from the Human Lung Cell Atlas
- [`incov_pbmc_subset`](incov_pbmc_subset.md) : Blood cells, surface
  proteins and plasma measurements of patients with COVID-19 from the
  INCOV cohort

## Association of an exposure with expression

Donor-pseudobulk models within cell types, and the models and analyses
built on them

- [`run_sc_exwas()`](run_sc_exwas.md) : Run cell-type-specific ExWAS
- [`run_multi_exwas()`](run_multi_exwas.md) : Run sc-ExWAS across
  multiple exposures
- [`run_interaction_test()`](run_interaction_test.md) : Test
  exposure-by-celltype interaction
- [`run_sc_exwas_pb_offset()`](run_sc_exwas_pb_offset.md) : Donor
  offset-pseudobulk negative-binomial sc-ExWAS
- [`run_sc_exwas_glmm()`](run_sc_exwas_glmm.md) : Cell-level linear
  mixed model sc-ExWAS with finite-df inference
- [`run_meta_exwas()`](run_meta_exwas.md) : Meta-analysis across cohorts
  for sc-ExWAS
- [`run_gsea()`](run_gsea.md) : Gene set enrichment on sc-ExWAS results

## Dose and response

- [`run_dose_response()`](run_dose_response.md) : Cell-type-specific
  exposure dose-response analysis
- [`run_dose_response_gam()`](run_dose_response_gam.md) : GAM-based
  dose-response analysis

## Mixtures of exposures

- [`run_sc_mixture()`](run_sc_mixture.md) : Cell-type-specific mixture
  exposure analysis
- [`run_mixture_qgcomp()`](run_mixture_qgcomp.md) : Formal mixture
  analysis via quantile g-computation

## Cell composition and gene-protein coupling

- [`run_exposure_composition()`](run_exposure_composition.md) : Test
  exposure associations with donor-level cell composition
- [`run_cell_coupling()`](run_cell_coupling.md) : Cell-Level
  Gene-Protein Coupling Analysis

## Exposure measurements and study planning

- [`exposure_impute_lod()`](exposure_impute_lod.md) : Impute exposure
  values below the limit of detection
- [`estimate_power()`](estimate_power.md) : Power analysis for sc-ExWAS
  study design

## Exploratory cross-omic networks

Estimation, comparison, simulation, classes, accessors and plots

- [`run_celltype_network()`](run_celltype_network.md) : Build
  cell-type-resolved cross-omic network
- [`run_exposure_network()`](run_exposure_network.md) : Network among
  exposure-associated features
- [`run_comparative_network()`](run_comparative_network.md) : Compare
  networks across cell types
- [`simulate_crossomic_network()`](simulate_crossomic_network.md) :
  Simulate cross-omic network data with known structure
- [`CelltypeNetwork-class`](CelltypeNetwork-class.md) : CelltypeNetwork:
  virtual base class for cell-type networks
- [`show(`*`<CelltypeNetworkResult>`*`)`](CelltypeNetworkResult-class.md)
  : CelltypeNetworkResult: cell-type-resolved cross-omic network
- [`show(`*`<NetworkComparison>`*`)`](NetworkComparison-class.md) :
  NetworkComparison: differential network analysis result
- [`CelltypeNetworkResult()`](CelltypeNetworkResult.md) : Create a
  CelltypeNetworkResult
- [`NetworkComparison()`](NetworkComparison.md) : Create a
  NetworkComparison
- [`networkPrecision()`](network-accessors.md)
  [`networkAdjacency()`](network-accessors.md)
  [`networkStability()`](network-accessors.md)
  [`networkNodes()`](network-accessors.md)
  [`networkCelltype()`](network-accessors.md)
  [`networkMethod()`](network-accessors.md)
  [`diffEdges()`](network-accessors.md)
  [`sharedEdges()`](network-accessors.md)
  [`celltypeSpecificEdges()`](network-accessors.md)
  [`comparisonSummary()`](network-accessors.md)
  [`metadata(`*`<CelltypeNetworkResult>`*`)`](network-accessors.md) :
  Accessors for network results
- [`plot_celltype_network()`](plot_celltype_network.md) : Plot
  cell-type-specific cross-omic network
- [`plot_network_comparison()`](plot_network_comparison.md) : Plot the
  edges of a network comparison
- [`plot_edge_stability()`](plot_edge_stability.md) : Edge stability
  plot

## Package

- [`exposomeSC-package`](exposomeSC-package.md)
  [`exposomeSC`](exposomeSC-package.md) : exposomeSC: Cell-Type-Resolved
  Exposome-Wide Association Studies for Single-Cell RNA-seq
