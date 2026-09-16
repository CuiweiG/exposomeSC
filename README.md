<div align="center">

# exposomeSC

*Study-aware, cell-type-resolved exposure association for population
single-cell transcriptomes*

[![License: Artistic-2.0](https://img.shields.io/badge/License-Artistic--2.0-blue.svg)](https://opensource.org/licenses/Artistic-2.0)

</div>

## What it does

`exposomeSC` tests donor-level exposures against gene expression within a cell
type, for single-cell RNA-seq cohorts in which each donor contributes many
cells. The exposure is a property of the donor, so cells from one donor are
measurements of the same unit rather than independent replicates. The package
therefore aggregates raw counts by donor and cell type and fits a robust edgeR
quasi-likelihood model per cell type, with the donor as the unit of inference.

The case study the package was built around is ever- versus never-smoking in
the integrated Human Lung Cell Atlas. This README quotes no result from it:
numbers enter the documentation only from result objects written by the
analysis itself, checked against the code and input manifests that produced
them, so that no figure in the text can drift from the run behind it.

## Installation

```r
# After Bioconductor acceptance:
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("exposomeSC")

# Development version from GitHub:
BiocManager::install("CuiweiG/exposomeSC")
```

## A first analysis

The package ships a small simulated data set: 2,000 cells from 20 donors in
three cell types, with five exposures and two covariates per donor. Ten genes
were given a PM2.5-dependent increase in monocytes when the data were made.

```r
library(exposomeSC)

sce <- readRDS(system.file("extdata", "demo_sce.rds", package = "exposomeSC"))
exposures <- as.matrix(read.csv(
    system.file("extdata", "demo_exposures.csv", package = "exposomeSC"),
    row.names = 1))

scee <- build_scee(sce, exposures, sample_col = "donor_id")

fit <- run_sc_exwas(scee,
    exposure = "PM2.5",
    covariates = c("age", "sex"),
    celltype_col = "cell_type",
    sample_col = "donor_id")

head(fit[order(fit$pvalue),
         c("gene", "celltype", "log2FC", "statistic", "pvalue", "padj_global")])
```

`build_scee()` checks that the counts are raw non-negative integers and that
every donor in the cells has a row in the exposure matrix. `run_sc_exwas()`
returns one row per gene and cell type, with `padj_global` adjusted across the
whole gene-by-cell-type family. On this data set the ten calls at
`padj_global < 0.05` are the ten genes given an effect, all in monocytes, and
there are none in the other two cell types.

For the default edgeR backend, `se` is `NA` by design: the quasi-likelihood F
test supplies no coefficient standard error, and `statistic` is the signed
square root of that F statistic, a ranking quantity rather than a Wald
statistic.

## What the estimate is

The primary estimand is the adjusted association between a donor-level
exposure and relative expression *within a cell type*. The workflow:

1. requires raw, non-negative integer counts and validated donor identifiers;
2. aggregates the sparse matrix by donor and cell type;
3. retains the complete target-cell-type transcriptome library size;
4. fixes expression support without inspecting exposure labels;
5. applies TMM normalisation and robust edgeR quasi-likelihood inference;
6. checks group support, design rank and numerical estimability; and
7. adjusts across the full feature-by-cell-type hypothesis family.

Cell-type composition is a separate relative-composition outcome. It is not an
offset for the primary expression estimand, and captured cell proportions are
not absolute tissue abundance.

## Functions

| Group | Functions |
|---|---|
| Container | `build_scee()`, `as_scee()` (from a rexposome `ExposomeSet`), `seurat_to_exposure()`, `exposure_impute_lod()`; accessors `exposureData()`, `exposureInfo()`, `exposureVariables()`, `cellSampleMap()` |
| Primary analysis | `run_sc_exwas()` (edgeR, DESeq2 or voom-dream backends), `run_sc_exwas_pb_offset()`, `run_multi_exwas()`, `run_sc_exwas_glmm()`, `run_meta_exwas()` |
| Design | `estimate_power()` |
| Exposure modelling | `run_dose_response()`, `run_dose_response_gam()`, `run_interaction_test()`, `run_mixture_qgcomp()`, `run_sc_mixture()` |
| Composition | `run_exposure_composition()` |
| Downstream | `run_gsea()`, `run_cell_coupling()` |
| Cross-omic networks (exploratory) | `simulate_crossomic_network()`, `run_celltype_network()`, `run_exposure_network()`, `run_comparative_network()`, `run_differential_network()`, `run_temporal_network()` and their plots |
| Experimental | `run_mediation()`, `run_causal_mediation()`, `run_network_mediation()`, `mediational_evalue()`, `compute_iers()`, `run_state_coupling()`, `run_spatial_exwas()` |

The network functions return model-dependent summaries, not validated
conditional-independence findings. The last row is this README's own grouping:
the functions in it are the ones whose estimator, resampling or identification
properties have been examined least. Three of them, `run_mediation()`,
`run_causal_mediation()` and `run_network_mediation()`, warn once per session
to say so, as does `run_comparative_network()` in its `fisher_z` mode. The
package does not supply causal identification by itself.

## Decomposition

`run_decomposed_exwas()`, `erd_evalue()` and `erd_resolution_sensitivity()`
have been removed. Their outcome was a within-cell-type expression rate, which
cannot identify the mechanical contribution of cell-type abundance to a tissue
mixture, so the direct/indirect labels and mediational E-values they produced
did not support causal or abundance claims. `run_erd_interventional()` is kept
as an unexported developmental function while a replacement estimand is worked
out.

## HLCA data provenance

The primary study uses the integrated HLCA normal-lung collection of Sikkema
*et al.* (2023), distributed under CC BY 4.0 through CZ CELLxGENE.

| Item | Pinned source |
|---|---|
| Publication | [Sikkema *et al.*, *Nature Medicine* (2023)](https://doi.org/10.1038/s41591-023-02327-2) |
| CELLxGENE collection | [Integrated Human Lung Cell Atlas](https://cellxgene.cziscience.com/collections/6f6d381a-7701-4781-935c-db10d30de293) |
| Exact H5AD asset | [688185ad-11c2-4172-a53a-f4f1f4076860.h5ad](https://datasets.cellxgene.cziscience.com/688185ad-11c2-4172-a53a-f4f1f4076860.h5ad) |
| Count layer | `/raw/X`; `/X` is not used for count modelling |
| H5AD SHA-256 | `1cbdde1e513a31bdb6a3f70296e2d0612760d2af4032e239c94596ea87b44007` |

`data-raw/fetch_hlca.R` rejects a missing raw layer, a changed asset size, and
non-finite, negative or non-integer values. From a clone of the repository,
after the pinned H5AD has been downloaded to `data/hlca_core.h5ad`:

```r
source(file.path("data-raw", "fetch_hlca.R"))
source(file.path("data-raw", "extract_hlca_design_metadata.R"))
source(file.path("data-raw", "extract_hlca_feature_metadata.R"))
```

The scripts derive the project root from their own location and do not call
`setwd()`.

## Documentation

- [Introduction](vignettes/introduction.Rmd): the estimand and the package
  contract
- [Quick start](vignettes/exposomeSC-quickstart.Rmd): a small smoking example
  end to end
- [Network tutorial](vignettes/network-tutorial.Rmd): exploratory cross-omic
  network estimation and descriptive comparison, run on simulated data
- [Simulation protocol](vignettes/simulation-benchmark.Rmd): how calibration
  and aligned comparison are to be done, with a complete-null check and a power
  approximation run in miniature

The continuous-integration configuration lives on the `ci` branch; the default
branch holds the package source only.

## Licence and citation

Artistic-2.0. Until a release is archived with a DOI, cite the exact Git commit
used together with the primary data publications.
