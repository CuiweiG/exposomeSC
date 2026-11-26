# Accessors for SingleCellExposomeExperiment objects

Get and set the donor-level exposure data of a
[`SingleCellExposomeExperiment`](SingleCellExposomeExperiment-class.md)
and the column of `colData` that links cells to donors.

## Usage

``` r
exposureData(x)

exposureInfo(x)

cellSampleMap(x)

sampleColumn(x)

sampleColumn(x) <- value

exposureVariables(x)

exposureData(x) <- value

exposureInfo(x) <- value

# S4 method for class 'SingleCellExposomeExperiment'
exposureData(x)

# S4 method for class 'SingleCellExposomeExperiment'
exposureInfo(x)

# S4 method for class 'SingleCellExposomeExperiment'
cellSampleMap(x)

# S4 method for class 'SingleCellExposomeExperiment'
sampleColumn(x)

# S4 method for class 'SingleCellExposomeExperiment'
sampleColumn(x) <- value

# S4 method for class 'SingleCellExposomeExperiment'
exposureVariables(x)

# S4 method for class 'SingleCellExposomeExperiment'
exposureData(x) <- value

# S4 method for class 'SingleCellExposomeExperiment'
exposureInfo(x) <- value
```

## Arguments

- x:

  A
  [`SingleCellExposomeExperiment`](SingleCellExposomeExperiment-class.md).

- value:

  For `exposureData<-`, a numeric matrix with one row per sample, named
  by sample ID, and one column per exposure; the rows of
  `exposureInfo(x)` follow the new columns, so an exposure can be added
  (its metadata row is filled with `NA`), dropped or reordered. For
  `exposureInfo<-`, a `DataFrame` whose `exposure` column equals
  `colnames(exposureData(x))`; for `sampleColumn<-`, the name of a
  `colData` column whose sample IDs match the rows of `exposureData(x)`.

## Value

- `exposureData()`:

  Numeric matrix, samples by exposures.

- `exposureInfo()`:

  `DataFrame` with one row per exposure.

- `exposureVariables()`:

  Character vector of exposure names.

- `sampleColumn()`:

  Character scalar: the `colData` column holding each cell's sample
  (donor) ID. The analysis functions read each cell's donor from this
  column.

- `cellSampleMap()`:

  `DataFrame` with columns `cell_id` (the cell names) and `sample_id`,
  derived from that column.

- Replacement functions:

  The modified object, after a validity check.

## Examples

``` r
library(SingleCellExperiment)
#> Loading required package: SummarizedExperiment
#> Loading required package: MatrixGenerics
#> Loading required package: matrixStats
#> 
#> Attaching package: ‘MatrixGenerics’
#> The following objects are masked from ‘package:matrixStats’:
#> 
#>     colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
#>     colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
#>     colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
#>     colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
#>     colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
#>     colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
#>     colWeightedMeans, colWeightedMedians, colWeightedSds,
#>     colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
#>     rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
#>     rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
#>     rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
#>     rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
#>     rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
#>     rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
#>     rowWeightedSds, rowWeightedVars
#> Loading required package: GenomicRanges
#> Loading required package: stats4
#> Loading required package: BiocGenerics
#> Loading required package: generics
#> 
#> Attaching package: ‘generics’
#> The following objects are masked from ‘package:base’:
#> 
#>     as.difftime, as.factor, as.ordered, intersect, is.element, setdiff,
#>     setequal, union
#> 
#> Attaching package: ‘BiocGenerics’
#> The following objects are masked from ‘package:stats’:
#> 
#>     IQR, mad, sd, var, xtabs
#> The following object is masked from ‘package:utils’:
#> 
#>     data
#> The following objects are masked from ‘package:base’:
#> 
#>     Filter, Find, Map, Position, Reduce, anyDuplicated, aperm, append,
#>     as.data.frame, basename, cbind, colnames, dirname, do.call,
#>     duplicated, eval, evalq, get, grep, grepl, is.unsorted, lapply,
#>     mapply, match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
#>     rank, rbind, rownames, sapply, saveRDS, scale, sequence, table,
#>     tapply, transform, unique, unsplit, which.max, which.min
#> Loading required package: S4Vectors
#> 
#> Attaching package: ‘S4Vectors’
#> The following object is masked from ‘package:utils’:
#> 
#>     findMatches
#> The following objects are masked from ‘package:base’:
#> 
#>     I, expand.grid, unname
#> Loading required package: IRanges
#> Loading required package: Seqinfo
#> Loading required package: Biobase
#> Welcome to Bioconductor
#> 
#>     Vignettes contain introductory material; view with
#>     'browseVignettes()'. To cite Bioconductor, see
#>     'citation("Biobase")', and for packages 'citation("pkgname")'.
#> 
#> Attaching package: ‘Biobase’
#> The following object is masked from ‘package:MatrixGenerics’:
#> 
#>     rowMedians
#> The following objects are masked from ‘package:matrixStats’:
#> 
#>     anyMissing, rowMedians
counts <- matrix(1L, nrow = 5, ncol = 12,
    dimnames = list(paste0("G", 1:5), paste0("c", 1:12)))
sce <- SingleCellExperiment(assays = list(counts = counts),
    colData = S4Vectors::DataFrame(
        cell_id = paste0("c", 1:12),
        donor_id = rep(paste0("D", 1:4), each = 3)))
exp_mat <- matrix(c(12, 25, 8, 30, 45, 32, 67, 51), nrow = 4,
    dimnames = list(paste0("D", 1:4), c("PM2.5", "age")))
scee <- build_scee(sce, exp_mat, sample_col = "donor_id")

exposureVariables(scee)
#> [1] "PM2.5" "age"  
exposureData(scee)
#>    PM2.5 age
#> D1    12  45
#> D2    25  32
#> D3     8  67
#> D4    30  51
exposureInfo(scee)
#> DataFrame with 2 rows and 3 columns
#>      exposure      family        unit
#>   <character> <character> <character>
#> 1       PM2.5          NA          NA
#> 2         age          NA          NA
sampleColumn(scee)
#> [1] "donor_id"
head(cellSampleMap(scee))
#> DataFrame with 6 rows and 2 columns
#>       cell_id   sample_id
#>   <character> <character>
#> 1          c1          D1
#> 2          c2          D1
#> 3          c3          D1
#> 4          c4          D2
#> 5          c5          D2
#> 6          c6          D2

## Replace the exposure matrix, here rescaling PM2.5
exp_new <- exposureData(scee)
exp_new[, "PM2.5"] <- exp_new[, "PM2.5"] / 10
exposureData(scee) <- exp_new
exposureData(scee)[, "PM2.5"]
#>  D1  D2  D3  D4 
#> 1.2 2.5 0.8 3.0 
```
