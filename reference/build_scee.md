# Build a SingleCellExposomeExperiment

Combines a `SingleCellExperiment` with sample-level exposure
measurements into a unified container for sc-ExWAS.

## Usage

``` r
build_scee(sce, exposure_matrix, sample_col, exposure_info = NULL)
```

## Arguments

- sce:

  A
  [`SingleCellExperiment`](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html).

- exposure_matrix:

  Numeric matrix. Rows = samples, columns = exposure variables, with
  unique row and column names. The row names must be exactly the sample
  IDs in `colData(sce)[[sample_col]]`: every sample of the SCE must have
  a row, and rows for samples absent from the SCE are an error. `NA` is
  allowed; `NaN` and infinite values are not.

- sample_col:

  Character. Column in `colData(sce)` containing the sample/donor ID.
  Its name is stored in the object (see
  [`sampleColumn`](SCEE-accessors.md)), and the analysis functions read
  each cell's donor from this column.

- exposure_info:

  `DataFrame` (optional). Metadata for exposures (for example family,
  unit and LOD). It must contain an `exposure` column equal to
  `colnames(exposure_matrix)`, in the same order.

## Value

A
[`SingleCellExposomeExperiment`](SingleCellExposomeExperiment-class.md).

## Examples

``` r
library(SingleCellExperiment)
library(S4Vectors)
set.seed(1)
counts <- matrix(rpois(2000, 5), nrow = 50,
    dimnames = list(paste0("G", 1:50), paste0("c", 1:40)))
sce <- SingleCellExperiment(
    assays = list(counts = counts),
    colData = DataFrame(
        cell_id = paste0("c", seq_len(40)),
        donor_id = rep(paste0("D", 1:4), each = 10),
        cell_type = rep(c("Mono", "NK"), 20)))
exp_mat <- matrix(rnorm(12), nrow = 4,
    dimnames = list(paste0("D", 1:4),
                    c("PM2.5", "Pb", "BPA")))
scee <- build_scee(sce, exp_mat, sample_col = "donor_id")
scee
#> class: SingleCellExposomeExperiment 
#> dim: 50 40 
#> metadata(0):
#> assays(1): counts
#> rownames(50): G1 G2 ... G49 G50
#> rowData names(0):
#> colnames(40): c1 c2 ... c39 c40
#> colData names(3): cell_id donor_id cell_type
#> reducedDimNames(0):
#> mainExpName: NULL
#> altExpNames(0):
#> exposureData(4 samples x 3 exposures): PM2.5, Pb, BPA
#> sampleColumn: donor_id
```
