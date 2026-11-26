# Plot the edges of a network comparison

Plots, for the selected edge set, how much each edge's partial
correlation varies across the compared cell types (`max_diff`), with the
per-cell-type edge counts and pairwise Jaccard overlap in the subtitle.
Up to 50 edges are shown.

## Usage

``` r
plot_network_comparison(
  comparison,
  highlight = c("differential", "shared", "all")
)
```

## Arguments

- comparison:

  [`NetworkComparison`](NetworkComparison-class.md).

- highlight:

  Character; which edges to plot: `"differential"` (default),
  `"shared"`, or `"all"`.

## Value

A `ggplot` object.

## Examples

``` r
shared_edge <- S4Vectors::DataFrame(
    node_i = "G1", node_j = "M1", category = "shared",
    max_diff = 0.1
)
differential_edge <- S4Vectors::DataFrame(
    node_i = "G2", node_j = "M1", category = "differential",
    max_diff = 0.5
)
comparison <- NetworkComparison(
    diff_edges = differential_edge,
    shared_edges = shared_edge,
    celltype_specific = list(),
    summary = list(
        celltypes = c("Monocyte", "B cell"),
        n_edges = c(Monocyte = 2L, `B cell` = 1L),
        n_shared = 1L,
        n_differential = 1L,
        jaccard = list(Monocyte_vs_B_cell = 0.5)
    )
)
if (requireNamespace("ggplot2", quietly = TRUE)) {
    comparison_plot <- plot_network_comparison(
        comparison,
        highlight = "differential"
    )
    inherits(comparison_plot, "ggplot")
}
#> [1] TRUE
```
