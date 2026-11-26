# exposomeSC: Cell-Type-Resolved Exposome-Wide Association Studies for Single-Cell RNA-seq

Tests donor-level exposures for association with gene expression within
cell types, with the donor as the unit of inference. Results are
associations; the package does not identify causal effects.

## Core analysis

- [`build_scee`](build_scee.md):

  Build a SingleCellExposomeExperiment from SCE + exposure matrix

- [`run_sc_exwas`](run_sc_exwas.md):

  Cell-type-specific ExWAS via pseudobulk (edgeR, DESeq2 or voom-dream
  backend)

- [`run_sc_exwas_pb_offset`](run_sc_exwas_pb_offset.md):

  Negative binomial model of the pseudobulk counts with a donor offset

- [`run_sc_exwas_glmm`](run_sc_exwas_glmm.md):

  Mixed model of the cells with donor random effects, linear by default

- [`run_multi_exwas`](run_multi_exwas.md):

  Batch ExWAS across multiple exposures with global FDR

- [`donor_pseudobulk`](donor_pseudobulk.md):

  Donor pseudobulk counts per cell type as SummarizedExperiment objects

## Dose-response and mixtures

- [`run_dose_response`](run_dose_response.md):

  Polynomial dose-response (AIC selection)

- [`run_dose_response_gam`](run_dose_response_gam.md):

  GAM dose-response via mgcv

- [`run_sc_mixture`](run_sc_mixture.md):

  Quantile-scored mixture screening

- [`run_mixture_qgcomp`](run_mixture_qgcomp.md):

  Quantile g-computation with qgcomp

- [`run_interaction_test`](run_interaction_test.md):

  Exposure x celltype interaction F-test

## Cell-level coupling

- [`run_cell_coupling`](run_cell_coupling.md):

  Gene-protein coupling via meta-regression with confounders

## Enrichment and meta-analysis

- [`run_gsea`](run_gsea.md):

  Pathway enrichment via fgsea

- [`run_meta_exwas`](run_meta_exwas.md):

  Cross-cohort meta-analysis (fixed/random effects)

- [`estimate_power`](estimate_power.md):

  Power/sample size calculator

## Composition

- [`run_exposure_composition`](run_exposure_composition.md):

  Exposure association with donor-level cell-type composition

## Cross-omic networks (exploratory)

- [`run_celltype_network`](run_celltype_network.md):

  Cell-type-resolved network of transcripts and donor-level molecular
  measurements

- [`run_exposure_network`](run_exposure_network.md):

  Network restricted to exposure-associated features

- [`run_comparative_network`](run_comparative_network.md):

  Descriptive comparison of networks across cell types

- [`simulate_crossomic_network`](simulate_crossomic_network.md):

  Simulated data with a known network

- [`plot_celltype_network`](plot_celltype_network.md),
  [`plot_network_comparison`](plot_network_comparison.md),
  [`plot_edge_stability`](plot_edge_stability.md):

  Plots

## Preprocessing

- [`exposure_impute_lod`](exposure_impute_lod.md):

  Below-LOD imputation

- [`as_scee`](as_scee.md):

  Convert from rexposome

- [`seurat_to_exposure`](seurat_to_exposure.md):

  Extract from Seurat

## References

Squair JW et al. (2021). Confronting false discoveries in single-cell
differential expression. *Nat Commun* 12:5692.
[doi:10.1038/s41467-021-25960-2](https://doi.org/10.1038/s41467-021-25960-2)

Aitchison J (1986). *The Statistical Analysis of Compositional Data*.
Chapman and Hall.

## See also

Useful links:

- <https://github.com/CuiweiG/exposomeSC>

- Report bugs at <https://github.com/CuiweiG/exposomeSC/issues>

## Author

**Maintainer**: Cuiwei Gao <48gaocuiwei@gmail.com>
([ORCID](https://orcid.org/0009-0002-4896-0601)) \[copyright holder\]
