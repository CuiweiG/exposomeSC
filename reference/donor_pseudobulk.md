# Donor pseudobulk counts for each cell type

Sums the raw counts of each cell type within each donor and returns one
`SummarizedExperiment` per cell type, with the donors' exposure data as
column data. These are the pseudobulk counts that
[`run_sc_exwas`](run_sc_exwas.md) and the other pseudobulk functions of
the package aggregate internally, so they can be passed to other
Bioconductor tools, inspected, or plotted.

## Usage

``` r
donor_pseudobulk(
  x,
  celltype_col = "cell_type",
  celltypes = NULL,
  min_cells = 10L
)
```

## Arguments

- x:

  A
  [`SingleCellExposomeExperiment`](SingleCellExposomeExperiment-class.md)
  with a `"counts"` assay. Each cell's donor is read from the `colData`
  column named by [`sampleColumn`](SCEE-accessors.md)`(x)`.

- celltype_col:

  Character. Column in `colData` containing cell type annotations.
  Default `"cell_type"`.

- celltypes:

  Character vector (optional). Cell types to aggregate. Default: all
  non-missing labels, sorted.

- min_cells:

  Integer. A donor is included for a cell type only when it has at least
  this many cells of that type. Default 10, the default of
  [`run_sc_exwas`](run_sc_exwas.md).

## Value

A
[`SimpleList`](https://rdrr.io/pkg/S4Vectors/man/SimpleList-class.html)
named by cell type. Each element is a
[`SummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
with one column per included donor and

- an assay `"counts"` with the summed raw counts (features by donors);

- `rowData` copied from `x`;

- `colData` with the donor identifier (`donor_id`), the cell type
  (`cell_type`), the number of cells summed (`n_cells`) and one column
  per variable of `exposureData(x)`;

- `metadata` recording the cell type, `min_cells` and the `colData`
  columns used for donors and cell types.

Cell types in which no donor reaches `min_cells` are left out and listed
in `metadata(result)$skipped`.

## Details

The summation is done by
[`scrapper::aggregateAcrossCells()`](https://rdrr.io/pkg/scrapper/man/aggregateAcrossCells.html).
Donors with missing exposure values are kept here; the analysis
functions drop them, and filter genes, when they fit a model, so a
fitted model may use fewer donors and genes than the object returned
here. Model diagnostics
([`metadata()`](https://rdrr.io/pkg/S4Vectors/man/Annotated-class.html)
of the `run_sc_exwas` result) list the donors and gene counts actually
used.

## References

Lun ATL, Kancherla J (2023). Powering single-cell analyses in the
browser with WebAssembly. *Journal of Open Source Software* 8(89):5603.
[doi:10.21105/joss.05603](https://doi.org/10.21105/joss.05603)

## See also

[`run_sc_exwas`](run_sc_exwas.md)

## Examples

``` r
library(SingleCellExperiment)
library(S4Vectors)
set.seed(1)
counts <- matrix(rpois(2000, 5), nrow = 20,
    dimnames = list(paste0("G", 1:20), paste0("c", 1:100)))
sce <- SingleCellExperiment(assays = list(counts = counts),
    colData = DataFrame(donor_id = rep(paste0("D", 1:5), each = 20),
        cell_type = rep(c("Mono", "NK"), 50)))
exp_mat <- matrix(rnorm(5), ncol = 1,
    dimnames = list(paste0("D", 1:5), "E1"))
scee <- build_scee(sce, exp_mat, sample_col = "donor_id")
pb <- donor_pseudobulk(scee, min_cells = 5L)
pb
#> List of length 2
#> names(2): Mono NK
SummarizedExperiment::colData(pb[["Mono"]])
#> DataFrame with 5 rows and 4 columns
#>       donor_id   cell_type   n_cells         E1
#>    <character> <character> <integer>  <numeric>
#> D1          D1        Mono        10  1.1349651
#> D2          D2        Mono        10  1.1119318
#> D3          D3        Mono        10 -0.8707776
#> D4          D4        Mono        10  0.2107316
#> D5          D5        Mono        10  0.0693956
```
