# Create a CelltypeNetworkResult

Wraps a network estimated elsewhere, or built by hand, in a
[`CelltypeNetworkResult`](CelltypeNetworkResult-class.md) so that it can
be compared with [`run_comparative_network`](run_comparative_network.md)
and plotted. The networks returned by
[`run_celltype_network`](run_celltype_network.md) and
[`run_exposure_network`](run_exposure_network.md) are already of this
class.

## Usage

``` r
CelltypeNetworkResult(
  precision,
  adjacency,
  node_info,
  celltype,
  method = "precomputed",
  stability = NULL,
  metadata = list()
)
```

## Arguments

- precision:

  Numeric square matrix: the precision (inverse covariance) matrix, one
  row and column per feature.

- adjacency:

  Symmetric matrix of edge indicators (1 = edge, 0 = no edge) with the
  dimensions of `precision` and a zero diagonal; a logical matrix is
  stored as 0/1 integers.

- node_info:

  `DataFrame` or `data.frame` with one row per feature, in matrix order,
  a `feature` column without missing values and an `omic_layer` column,
  which the summaries, plots and comparisons read; the summaries and
  plots recognise the layers `"transcript"` and `"metabolite"`. Row and
  column names of the matrices, where present, must equal `feature`, so
  that a misordered table is an error.

- celltype:

  Character scalar: the cell type.

- method:

  Character scalar: how the network was estimated. Default
  `"precomputed"`.

- stability:

  Optional matrix of edge selection frequencies between 0 and 1 with the
  dimensions of `adjacency`; `NULL` (default) stores none.

- metadata:

  List of further information, for example `n_donors`, the number of
  donors as a whole number, which the `show` method and
  [`plot_celltype_network`](plot_celltype_network.md) print.

## Value

A valid `CelltypeNetworkResult`.

## See also

[`network-accessors`](network-accessors.md)

## Examples

``` r
features <- c("G1", "G2", "M1")
precision <- diag(3)
precision[1, 3] <- precision[3, 1] <- -0.4
adjacency <- (precision != 0) * 1L
diag(adjacency) <- 0L
dimnames(precision) <- dimnames(adjacency) <- list(features, features)
network <- CelltypeNetworkResult(
    precision = precision,
    adjacency = adjacency,
    node_info = S4Vectors::DataFrame(
        feature = features,
        omic_layer = c("transcript", "transcript", "metabolite")),
    celltype = "Monocyte",
    metadata = list(n_donors = 20L))
network
#> CelltypeNetworkResult for 'Monocyte'
#>   Method: precomputed
#>   Nodes: 3 (2 transcripts, 1 metabolites)
#>   Edges: 1 (1 cross-omic)
#>   Donors: 20
```
