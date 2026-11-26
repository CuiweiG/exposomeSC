# SingleCellExposomeExperiment: Container for single-cell exposome data

Extends
[`SingleCellExperiment`](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
with sample-level environmental exposure measurements. This enables
cell-type-specific exposome-wide association studies (sc-ExWAS).

## Usage

``` r
# S4 method for class 'SingleCellExposomeExperiment,ANY,ANY,ANY'
x[i, j, ..., drop = TRUE]

# S4 method for class 'SingleCellExposomeExperiment'
cbind(..., deparse.level = 1)

# S4 method for class 'SingleCellExposomeExperiment'
rbind(..., deparse.level = 1)

# S4 method for class 'SingleCellExposomeExperiment,DataFrame'
colData(x, ...) <- value

# S4 method for class 'SingleCellExposomeExperiment,NULL'
colData(x, ...) <- value

# S4 method for class 'SingleCellExposomeExperiment,list'
dimnames(x) <- value

# S4 method for class 'SingleCellExposomeExperiment,NULL'
dimnames(x) <- value

# S4 method for class 'SingleCellExposomeExperiment'
show(object)
```

## Arguments

- x:

  A SingleCellExposomeExperiment.

- i:

  Row (gene) subscript.

- j:

  Column (cell) subscript.

- ...:

  Additional arguments.

- drop:

  Logical. Ignored (kept for S4 compatibility).

- deparse.level:

  Ignored.

- value:

  The replacement `colData` or `dimnames`.

- object:

  A SingleCellExposomeExperiment.

## Value

The subsetting method returns a valid `SingleCellExposomeExperiment`.
When cells are subset, its `exposureData` is restricted to the donors
that keep at least one cell; feature subsetting leaves the donor-level
slots unchanged. `cbind` combines objects with the same exposure
variables, exposure metadata and sample column; donors present in more
than one object must have identical exposure values. Every argument to
`cbind` or `rbind` must be a `SingleCellExposomeExperiment`: when a
plain `SingleCellExperiment` comes first, R dispatches to the
`SingleCellExperiment` method, which returns a `SingleCellExperiment`
without the exposure data (with only a warning about ambiguous
dispatch), so convert such objects with [`build_scee`](build_scee.md)
first. `rbind` combines objects with the same cells and different
features; they must have the same sample column and `exposureInfo`, and
the same exposure values for each donor. The replacement methods return
the modified object. The `show` method is called for its display side
effect and returns `NULL` invisibly.

## Details

The donor of each cell is read from one column of `colData`, whose name
is stored in the object (see [`sampleColumn`](SCEE-accessors.md)). The
cell-to-donor mapping is not stored a second time, so renaming cells,
combining objects and editing `colData` cannot leave two copies out of
step. The object is valid when every donor in that column has exactly
one row in `exposureData` and every row of `exposureData` belongs to at
least one cell.

Replacing `colData`, or a column of it with `$<-` or `[[<-`, checks
validity again, so an edit that removes the donor column or introduces a
donor without exposure data is an error; so does replacing the cell
names with `colnames<-`. To change donor identifiers, edit
`exposureData` and `colData` together and rebuild the object with
[`build_scee`](build_scee.md).

## Slots

- `exposureData`:

  Numeric matrix. Rows = samples (donors), columns = exposure variables.

- `exposureInfo`:

  `DataFrame`. Metadata for each exposure (family, unit, LOD).

- `sampleColumn`:

  Character scalar. Name of the `colData` column holding each cell's
  sample (donor) identifier.

## Data hierarchy


    Donor (n=50-1000)
      |-- Exposures (PM2.5, Pb, BPA, ...)   <- sample level
      +-- Cells x Genes matrix
          |-- T cells (1000 cells)
          |-- NK cells (500 cells)            <- cell level
          +-- Monocytes (800 cells)

## References

van der Wijst MGP et al. (2020). The single-cell eQTLGen consortium.
*eLife* 9:e52155.
[doi:10.7554/eLife.52155](https://doi.org/10.7554/eLife.52155)

## Examples

``` r
showClass("SingleCellExposomeExperiment")
#> Class "SingleCellExposomeExperiment" [package "exposomeSC"]
#> 
#> Slots:
#>                                                                 
#> Name:                  exposureData                 exposureInfo
#> Class:                       matrix                    DataFrame
#>                                                                 
#> Name:                  sampleColumn          int_elementMetadata
#> Class:                    character                    DataFrame
#>                                                                 
#> Name:                   int_colData                 int_metadata
#> Class:                    DataFrame                         list
#>                                                                 
#> Name:                     rowRanges                      colData
#> Class: GenomicRanges_OR_GRangesList                    DataFrame
#>                                                                 
#> Name:                        assays                        NAMES
#> Class:               Assays_OR_NULL            character_OR_NULL
#>                                                                 
#> Name:               elementMetadata                     metadata
#> Class:                    DataFrame                         list
#> 
#> Extends: 
#> Class "SingleCellExperiment", directly
#> Class "RangedSummarizedExperiment", by class "SingleCellExperiment", distance 2
#> Class "SummarizedExperiment", by class "SingleCellExperiment", distance 3
#> Class "RectangularData", by class "SingleCellExperiment", distance 4
#> Class "Vector", by class "SingleCellExperiment", distance 4
#> Class "Annotated", by class "SingleCellExperiment", distance 5
#> Class "vector_OR_Vector", by class "SingleCellExperiment", distance 5
```
