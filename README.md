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

The case study the package was built around is smoking in the integrated
Human Lung Cell Atlas, and a small part of those data is shipped with it. This
README quotes no result; the tables and figures of the vignettes are computed
from the shipped data when they are built.

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

The package ships a small real data set: 352 tracheal epithelial cells from 12
donors of the Human Lung Cell Atlas, six heavy smokers and six never smokers,
with up to 10 basal, secretory and multiciliated-lineage cells per donor.
`?hlca_airway_subset` gives its source, licence and selection rule. The data
set shows how the functions are called on real counts; with 12 donors and at
most 10 cells per donor and cell type it is far too small to support
conclusions about smoking. A second real data set, `?incov_pbmc_subset`, holds
blood cells of 43 patients with their surface proteins, plasma metabolites and
plasma levels of exogenous compounds, for the functions that need those.

```r
library(exposomeSC)

scee <- readRDS(system.file("extdata", "hlca_airway_subset.rds",
    package = "exposomeSC"))
scee

fit <- run_sc_exwas(scee,
    exposure = "smoking",
    covariates = "chemistry_v3")

head(fit[order(fit$pvalue),
         c("gene", "celltype", "log2FC", "statistic", "pvalue", "padj_global")])
```

The object was made by `build_scee()`, which checks that every donor in the
cells has a row in the exposure matrix and that every row belongs to a donor
with cells; the script that made the file checks that its counts are
non-negative whole numbers. `run_sc_exwas()` returns one row per gene and cell
type, with `padj_global` adjusted across the whole gene-by-cell-type family.

For the default edgeR backend, `se` is `NA` by design: the quasi-likelihood F
test supplies no coefficient standard error, and `statistic` is the signed
square root of that F statistic, a ranking quantity rather than a Wald
statistic.

## What the estimate is

The primary estimand is the adjusted association between a donor-level
exposure and relative expression *within a cell type*. The workflow:

1. takes raw counts and the donor of each cell from a container that has
   checked the donors against the exposure table;
2. aggregates the sparse matrix by donor and cell type;
3. retains the complete target-cell-type transcriptome library size;
4. fixes expression support without inspecting exposure labels;
5. applies TMM normalisation and robust edgeR quasi-likelihood inference;
6. checks group support, design rank and numerical estimability; and
7. adjusts across the full feature-by-cell-type hypothesis family.

The functions sum the counts as they are given and do not check that they are
raw: counts that were normalised or corrected beforehand are summed like any
others, and the models are then not the ones described here.

Cell-type composition is a separate relative-composition outcome. It is not an
offset for the primary expression estimand, and captured cell proportions are
not absolute tissue abundance.

## Functions

| Group | Functions |
|---|---|
| Container | `build_scee()`, `as_scee()` (from a rexposome `ExposomeSet`), `seurat_to_exposure()`, `exposure_impute_lod()`; accessors `exposureData()`, `exposureInfo()`, `exposureVariables()`, `sampleColumn()`, `cellSampleMap()` |
| Primary analysis | `run_sc_exwas()` (edgeR, DESeq2 or voom-dream backends), `donor_pseudobulk()`, `run_sc_exwas_pb_offset()`, `run_multi_exwas()`, `run_sc_exwas_glmm()`, `run_meta_exwas()` |
| Design | `estimate_power()` |
| Exposure modelling | `run_dose_response()`, `run_dose_response_gam()`, `run_interaction_test()`, `run_mixture_qgcomp()`, `run_sc_mixture()` |
| Composition | `run_exposure_composition()` |
| Downstream | `run_gsea()`, `run_cell_coupling()` |
| Cross-omic networks (exploratory) | `simulate_crossomic_network()`, `run_celltype_network()`, `run_exposure_network()`, `run_comparative_network()` and their plots |

The network functions fit sparse Gaussian graphical models with the graphical
lasso and the collaborative graphical lasso. An edge is a summary under that
model for the donors analysed, not a tested conditional independence, and
`run_comparative_network()` describes the differences between networks
without testing them. The package estimates associations; it does not supply
causal identification.

## HLCA data provenance

The airway data come from the integrated HLCA collection of Sikkema *et al.*
(2023), distributed under CC BY 4.0 through CZ CELLxGENE.

| Item | Pinned source |
|---|---|
| Publication | [Sikkema *et al.*, *Nature Medicine* (2023)](https://doi.org/10.1038/s41591-023-02327-2) |
| CELLxGENE collection | [Integrated Human Lung Cell Atlas](https://cellxgene.cziscience.com/collections/6f6d381a-7701-4781-935c-db10d30de293) |
| Exact H5AD asset | [688185ad-11c2-4172-a53a-f4f1f4076860.h5ad](https://datasets.cellxgene.cziscience.com/688185ad-11c2-4172-a53a-f4f1f4076860.h5ad) |
| Count layer | `/raw/X`; `/X` is not used for count modelling |
| H5AD SHA-256 | `1cbdde1e513a31bdb6a3f70296e2d0612760d2af4032e239c94596ea87b44007` |

The example data set shipped with the package,
`inst/extdata/hlca_airway_subset.rds`, is taken from this file: 352 cells of
12 donors of its study `Seibold_2020`, which is
[Goldfarbmuren *et al.*, *Nature Communications* (2020)](https://doi.org/10.1038/s41467-020-16239-z).
`inst/scripts/make_hlca_airway_subset.R` checks the size and the SHA-256 of
the H5AD and makes the subset from it; cells and genes are selected and the
counts are unchanged. The subset is redistributed under the same CC BY 4.0
licence.

The scripts that prepare the full ever- versus never-smoking analysis object
are installed with the package, in the directory that
`system.file("scripts", package = "exposomeSC")` returns. `fetch_hlca.R`
rejects a missing raw layer, a changed asset size, and non-finite, negative or
non-integer values. Each script takes the paths of its inputs and of its
output as arguments; `create_hlca_input_manifest.R` takes first the directory,
here `.`, in which it writes `hlca_input_manifest.rds` and
`hlca_input_manifest.csv`. With the pinned H5AD downloaded as
`hlca_core.h5ad`:

```sh
scripts=$(Rscript -e 'cat(system.file("scripts", package = "exposomeSC"))')
Rscript "$scripts/fetch_hlca.R" hlca_core.h5ad hlca_smoking_scee.rds
Rscript "$scripts/extract_hlca_design_metadata.R" hlca_core.h5ad \
    hlca_design_metadata.rds
Rscript "$scripts/extract_hlca_feature_metadata.R" hlca_core.h5ad \
    hlca_feature_metadata.rds
Rscript "$scripts/create_hlca_input_manifest.R" . hlca_core.h5ad \
    hlca_smoking_scee.rds hlca_design_metadata.rds hlca_feature_metadata.rds
```

The scripts need rhdf5, zellkonverter, AnnotationDbi, org.Hs.eg.db and
digest, which the package itself does not depend on.

## INCOV data provenance

The second data set, `inst/extdata/incov_pbmc_subset.rds`, is taken from the
INCOV cohort of
[Su *et al.*, *Cell* (2020)](https://doi.org/10.1016/j.cell.2020.10.037):
1,032 blood mononuclear cells of 43 patients with COVID-19 at their first
blood draw, with the counts of 192 antibodies for the same cells (four of them
isotype controls), 182 plasma metabolites and the plasma levels of 13
compounds of exogenous origin.

| Item | Source |
|---|---|
| Clinical records, plasma metabolites | Table S1 of [Mendeley Data, doi:10.17632/tzydswhhb5.5](https://doi.org/10.17632/tzydswhhb5.5), version 5, CC BY 4.0 |
| Single-cell genes and surface proteins | [ArrayExpress E-MTAB-9357](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-9357), redistributed under CC BY 4.0 with the permission of the authors who deposited it |

`inst/scripts/make_incov_pbmc_subset.R` checks the size and the SHA-256 of the
88 input files, recovers the counts from the logarithms that the archive
holds, assigns monocytes and T cells by gates on surface proteins and marker
genes, and draws the cells; `?incov_pbmc_subset` describes each step. It takes
the folder with the input files and the output file as arguments and needs
data.table and R.utils.

## Documentation

- [exposomeSC](vignettes/exposomeSC.Rmd): the container, the association
  model and the further analyses, each run on one of the two real data sets
  of the package
- [Network tutorial](vignettes/network-tutorial.Rmd): exploratory cross-omic
  networks of two cell types and plasma metabolites, and their descriptive
  comparison, run on the INCOV data
- [Permutation and power](vignettes/permutation-and-power.Rmd): checks of a
  result on the airway data, by permuting the exposure among the donors and
  by a second backend, the two adjustments of the p-values, and a power
  approximation

On every push to `main`, GitHub Actions run `R CMD check` on Linux, macOS and
Windows, `BiocCheck`, the test coverage and a check of the citation, the links
and the spelling; the workflows are in `.github/workflows`.

## Licence and citation

The code and the documentation of the package are under Artistic-2.0. The two
data files in `inst/extdata` are not: they are parts of the data of others and
are redistributed under CC BY 4.0. `hlca_airway_subset.rds` keeps the licence
of the atlas. In `incov_pbmc_subset.rds` the clinical records and plasma
measurements are under CC BY 4.0 as published on Mendeley Data, and for the
single-cell counts, from ArrayExpress E-MTAB-9357, the authors who deposited
them gave this permission on 21 November 2026:

> The single-cell values from E-MTAB-9357 included in exposomeSC may be used
> and redistributed under CC BY 4.0, with citation of Su et al. (2020).

The two provenance sections above and the help pages of the data sets name the
sources and the publications to cite.

To cite the package, use `citation("exposomeSC")`: it gives the
Bioconductor DOI,
[10.18129/B9.bioc.exposomeSC](https://doi.org/10.18129/B9.bioc.exposomeSC),
and the version installed. When you use the shipped data sets, please
cite their primary publications as well.
