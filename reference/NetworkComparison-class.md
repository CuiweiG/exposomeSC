# NetworkComparison: differential network analysis result

Stores results from
[`run_comparative_network`](run_comparative_network.md), including
differential, shared, and cell-type-specific edges. Use the accessors in
[`network-accessors`](network-accessors.md) to extract them.

## Usage

``` r
# S4 method for class 'NetworkComparison'
show(object)
```

## Arguments

- object:

  A `NetworkComparison` object.

## Value

The `show` method prints a summary and returns `NULL` invisibly.

## Slots

- `diff_edges`:

  `DataFrame`. Edges present in some but not all of the compared
  networks (including edges unique to one cell type), with the range of
  their partial correlations.

- `shared_edges`:

  `DataFrame`. Edges present in all compared cell types.

- `celltype_specific`:

  Named list of `DataFrame`s. Edges unique to each cell type.

- `summary`:

  List. Summary statistics: Jaccard similarity, edge counts, cell types
  compared.

## See also

[`network-accessors`](network-accessors.md),
[`NetworkComparison`](NetworkComparison.md)

## Examples

``` r
showClass("NetworkComparison")
#> Class "NetworkComparison" [package "exposomeSC"]
#> 
#> Slots:
#>                                                                               
#> Name:         diff_edges      shared_edges celltype_specific           summary
#> Class:         DataFrame         DataFrame              list              list
```
