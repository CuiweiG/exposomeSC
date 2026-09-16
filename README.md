<div align="center">

# exposomeSC

*Study-aware, cell-type-resolved exposure association for population
single-cell transcriptomes*

[![License: Artistic-2.0](https://img.shields.io/badge/License-Artistic--2.0-blue.svg)](https://opensource.org/licenses/Artistic-2.0)

</div>

## Status

`exposomeSC` is research software under active development. The current case
study and method evaluation concern ever- versus never-smoking in the
integrated Human Lung Cell Atlas (HLCA). Primary, sensitivity, empirical-null,
composition and external-targeted-support outputs are governed by
content-addressed analysis contracts.

This README deliberately does not transcribe discovery counts, effect
estimates, calibration statistics or biological conclusions. Those quantities
will be inserted only by a reporting script that validates the source result,
code manifest, input manifest and analysis-contract digest. Earlier validation
and decomposition headlines are withdrawn and must not be cited.

## Statistical target

The exposure is measured at donor level; cells from one donor are therefore
not independent experimental replicates. The primary estimand is the adjusted
exposure association with relative expression *within a cell type*. The
canonical workflow:

1. requires raw, non-negative integer counts and validated donor identifiers;
2. aggregates the sparse matrix by donor and cell type;
3. retains the complete target-cell-type transcriptome library size;
4. fixes expression support without inspecting exposure labels;
5. applies TMM normalisation and robust edgeR quasi-likelihood inference;
6. checks group support, design rank and numerical estimability; and
7. adjusts across the analysis-contract feature-by-cell-type hypothesis family.

Cell-type composition is a separate relative-composition outcome. It is not
an offset for the primary expression estimand, and captured cell proportions
are not absolute tissue abundance.

## Minimal package interface

The following example assumes that `sce` contains raw counts and donor/cell-type
annotations, while `donor_design` is a numeric donor-by-variable matrix with
matching row names. Study indicators and spline terms must be constructed
before `build_scee()`.

```r
library(exposomeSC)

scee <- build_scee(
    sce,
    donor_design,
    sample_col = "donor_id"
)

study_terms <- grep("^study_", exposureVariables(scee), value = TRUE)
fit <- run_sc_exwas(
    scee,
    exposure = "smoking_ever",
    celltype_col = "ann_level_3",
    sample_col = "donor_id",
    covariates = c(
        "age_ns1", "age_ns2", "age_ns3", "sex_binary", study_terms
    ),
    min_cells = 20L,
    min_donors = 30L,
    min_group_donors = 10L,
    method = "edgeR",
    filter_method = "fixed_support",
    filter_min_cpm = 1,
    filter_min_donors = 10L,
    filter_min_total_count = 15L,
    robust = TRUE
)
```

This is an interface illustration. A complete analysis also constructs the
covariates, checks immutable inputs, annotates feature families, records code
identity and writes outputs atomically.

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

`data-raw/fetch_hlca.R` rejects a missing raw layer, a changed asset size,
and non-finite, negative or non-integer values.

## Building the HLCA input object

From a clone of the repository, after the pinned H5AD has been downloaded to
`data/hlca_core.h5ad`:

```r
source(file.path("data-raw", "fetch_hlca.R"))
source(file.path("data-raw", "extract_hlca_design_metadata.R"))
source(file.path("data-raw", "extract_hlca_feature_metadata.R"))
```

The scripts derive the project root from their own location; they do not
require `setwd()`.

## Decomposition and causal language

`run_decomposed_exwas()`, `erd_evalue()` and `erd_resolution_sensitivity()`
have been removed. The earlier ERD outcome was a within-cell-type expression
rate and could not identify the mechanical contribution of cell-type abundance
to a tissue mixture, so its direct/indirect labels and mediational E-values did
not support causal or abundance claims. `run_erd_interventional()` is kept as
an unexported developmental function and is not for confirmatory
interpretation.

Network comparison, network mediation, causal mediation, interventional ERD and IERS
are experimental pending estimator, resampling, selection and
identification repairs. The package does not supply causal identification by
itself.

## Reproducibility and reporting

Canonical outputs record:

- verified input byte sizes and SHA-256 digests;
- the analysis script and ordered package-code manifest;
- model parameters and design diagnostics;
- session information and deterministic random seeds;
- the analysis-contract digest; and
- output MD5 and SHA-256 digests.

Figures and tables must read frozen result objects. Numerical values must not
be copied into plotting or reporting code. Public release additionally
requires a tagged source archive, an archival DOI and a clean package check;
these release fields remain pending.

## Documentation

- [`Introduction`](vignettes/introduction.Rmd): estimand and package contract
- [`Quick start`](vignettes/exposomeSC-quickstart.Rmd): concise smoking example
- [`Network tutorial`](vignettes/network-tutorial.Rmd): exploratory cross-omic
  network estimation and descriptive comparison
- [`Simulation protocol`](vignettes/simulation-benchmark.Rmd): calibration and
  aligned-comparison requirements without copied numerical claims

The continuous-integration configuration lives on the `ci` branch; the default
branch holds the package source only.

## Licence and citation

The source is available under the Artistic-2.0 licence. A release DOI and
formal software citation are **pending the frozen public release**. Until then,
cite the exact Git commit used and the primary data publications.
