# Convert Seurat metadata to exposure matrix

Extracts donor-level metadata from a Seurat object's `meta.data` and
formats it as an exposure matrix for [`build_scee`](build_scee.md).

## Usage

``` r
seurat_to_exposure(seurat_meta, sample_col, exposure_cols)
```

## Arguments

- seurat_meta:

  A data.frame (from `seurat_obj@meta.data`).

- sample_col:

  Character. Column with donor IDs.

- exposure_cols:

  Character vector. Columns to use as exposures.

## Value

A numeric matrix (donors x exposures). Each value is the mean over the
donor's cells; a warning names columns whose values vary within a donor,
and non-numeric columns are an error.

## Examples

``` r
meta <- data.frame(
    donor = rep(c("D1","D2","D3"), each = 10),
    PM2.5 = rep(c(12, 25, 8), each = 10),
    age = rep(c(45, 32, 67), each = 10))
exp_mat <- seurat_to_exposure(meta,
    sample_col = "donor",
    exposure_cols = c("PM2.5", "age"))
dim(exp_mat)
#> [1] 3 2
```
