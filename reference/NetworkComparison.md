# Create a NetworkComparison

Builds a [`NetworkComparison`](NetworkComparison-class.md) from edge
tables. Most users obtain one from
[`run_comparative_network`](run_comparative_network.md) instead.

## Usage

``` r
NetworkComparison(
  diff_edges = S4Vectors::DataFrame(),
  shared_edges = S4Vectors::DataFrame(),
  celltype_specific = list(),
  summary = list()
)
```

## Arguments

- diff_edges, shared_edges:

  `DataFrame`s with one row per edge and at least the columns `node_i`
  and `node_j`.
  [`plot_network_comparison()`](plot_network_comparison.md) also uses
  `category` and `max_diff`. Edges in `shared_edges` have category
  `"shared"`.

- celltype_specific:

  Named list of `DataFrame`s, one per cell type, with the edges found
  only in that cell type.

- summary:

  List with, for example, `celltypes` (the compared cell types),
  `n_edges` (named edge counts) and `jaccard` (named pairwise overlaps).

## Value

A valid `NetworkComparison`.

## See also

[`network-accessors`](network-accessors.md)

## Examples

``` r
shared <- S4Vectors::DataFrame(node_i = "G1", node_j = "M1",
    category = "shared", max_diff = 0.1)
differential <- S4Vectors::DataFrame(node_i = "G2", node_j = "M1",
    category = "differential", max_diff = 0.5)
NetworkComparison(
    diff_edges = differential,
    shared_edges = shared,
    summary = list(celltypes = c("Monocyte", "B cell")))
#> NetworkComparison
#>   Cell types compared: Monocyte, B cell
#>   Shared edges: 1
#>   Differential edges: 1
```
