# Convert rexposome ExposomeSet to exposomeSC

Bridges the bulk exposome ecosystem with single-cell analysis by
extracting the exposure matrix from an `ExposomeSet` (rexposome) and
combining it with a `SingleCellExperiment`.

## Usage

``` r
as_scee(exposome_set, sce, sample_col, exposures = NULL)
```

## Arguments

- exposome_set:

  An `ExposomeSet` from rexposome.

- sce:

  A `SingleCellExperiment` with matching donor/sample IDs.

- sample_col:

  Character. Column in `colData(sce)` with sample IDs matching
  `sampleNames(exposome_set)`.

- exposures:

  Character vector (optional). Subset of exposures to include; a name
  that the ExposomeSet does not have is an error. Default: all.

## Value

A
[`SingleCellExposomeExperiment`](SingleCellExposomeExperiment-class.md).

## Details

This function enables a common workflow: researchers have existing bulk
exposome data in rexposome format and want to add single-cell
resolution. The exposure matrix and metadata are extracted from the
`ExposomeSet` and combined with the SCE to create an SCEE.

Every sample in `colData(sce)[[sample_col]]` must be present in the
ExposomeSet; ExposomeSet samples without cells are dropped.

## References

Hernandez-Ferrer C et al. (2019). Comprehensive study of the exposome
and omic data using rexposome Bioconductor packages. *Bioinformatics*
35:5344-5345.

## Examples

``` r
if (requireNamespace("Biobase", quietly = TRUE)) {
    ## as_scee() needs only what an ExposomeSet is: a Biobase eSet whose
    ## assayData element "exp" holds exposures (rows) by samples
    ## (columns). A minimal eSet of that shape stands in here; an object
    ## from rexposome::loadExposome() is converted the same way.
    setClass("ExposomeSetExample", contains = "eSet")
    samples <- paste0("D", 1:4)
    exposures <- matrix(c(12, 25, 8, 30, 0.4, 0.9, 0.2, 1.1), nrow = 2,
        byrow = TRUE, dimnames = list(c("PM2.5", "Pb"), samples))
    es <- methods::new("ExposomeSetExample",
        assayData = Biobase::assayDataNew("environment", exp = exposures),
        phenoData = Biobase::AnnotatedDataFrame(
            data.frame(sex = c(0, 1, 0, 1), row.names = samples)),
        featureData = Biobase::AnnotatedDataFrame(
            data.frame(family = c("air", "metal"),
                       row.names = c("PM2.5", "Pb"))))
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = matrix(1L, nrow = 5, ncol = 8,
            dimnames = list(paste0("G", 1:5), paste0("c", 1:8)))),
        colData = S4Vectors::DataFrame(
            cell_id = paste0("c", 1:8),
            donor_id = rep(samples, each = 2)))
    scee <- as_scee(es, sce, sample_col = "donor_id")
    exposureData(scee)
    exposureInfo(scee)
}
#> DataFrame with 2 rows and 2 columns
#>      exposure      family
#>   <character> <character>
#> 1       PM2.5         air
#> 2          Pb       metal
```
