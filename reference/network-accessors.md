# Accessors for network results

Extract the parts of a
[`CelltypeNetworkResult`](CelltypeNetworkResult-class.md) (returned by
[`run_celltype_network`](run_celltype_network.md) and
[`run_exposure_network`](run_exposure_network.md)) and of a
[`NetworkComparison`](NetworkComparison-class.md) (returned by
[`run_comparative_network`](run_comparative_network.md)).

## Usage

``` r
networkPrecision(x)

networkAdjacency(x)

networkStability(x)

networkNodes(x)

networkCelltype(x)

networkMethod(x)

diffEdges(x)

sharedEdges(x)

celltypeSpecificEdges(x)

comparisonSummary(x)

# S4 method for class 'CelltypeNetworkResult'
networkPrecision(x)

# S4 method for class 'CelltypeNetworkResult'
networkAdjacency(x)

# S4 method for class 'CelltypeNetworkResult'
networkStability(x)

# S4 method for class 'CelltypeNetworkResult'
networkNodes(x)

# S4 method for class 'CelltypeNetworkResult'
networkCelltype(x)

# S4 method for class 'CelltypeNetworkResult'
networkMethod(x)

# S4 method for class 'CelltypeNetworkResult'
metadata(x, ...)

# S4 method for class 'NetworkComparison'
diffEdges(x)

# S4 method for class 'NetworkComparison'
sharedEdges(x)

# S4 method for class 'NetworkComparison'
celltypeSpecificEdges(x)

# S4 method for class 'NetworkComparison'
comparisonSummary(x)
```

## Arguments

- x:

  A `CelltypeNetworkResult` or, for `diffEdges()`, `sharedEdges()`,
  `celltypeSpecificEdges()` and `comparisonSummary()`, a
  `NetworkComparison`.

- ...:

  Ignored by
  [`metadata()`](https://rdrr.io/pkg/S4Vectors/man/Annotated-class.html).

## Value

- `networkPrecision()`:

  Numeric feature-by-feature precision matrix.

- `networkAdjacency()`:

  Feature-by-feature matrix of edge indicators (1 = edge, 0 = no edge).

- `networkStability()`:

  Feature-by-feature matrix of edge selection frequencies, or a 0 by 0
  matrix when none were stored.

- `networkNodes()`:

  `DataFrame` with one row per feature, in the order of the matrix rows.

- `networkCelltype()`:

  Character scalar: the cell type.

- `networkMethod()`:

  Character scalar: the estimation method.

- [`metadata()`](https://rdrr.io/pkg/S4Vectors/man/Annotated-class.html):

  List of parameters and summaries, for example `n_donors`, `transform`
  and `edge_counts`.

- `diffEdges()`:

  `DataFrame` of edges that are unique to one cell type or present in
  some but not all networks.

- `sharedEdges()`:

  `DataFrame` of edges present in every network.

- `celltypeSpecificEdges()`:

  Named list of `DataFrame`s, one per cell type with edges found only in
  that cell type.

- `comparisonSummary()`:

  List with the compared cell types, edge counts and pairwise Jaccard
  overlap.

## Examples

``` r
features <- c("G1", "G2", "M1")
make_network <- function(celltype, cross_weight) {
    precision <- diag(3)
    precision[1, 3] <- precision[3, 1] <- cross_weight
    adjacency <- (precision != 0) * 1L
    diag(adjacency) <- 0L
    dimnames(precision) <- dimnames(adjacency) <- list(features, features)
    CelltypeNetworkResult(
        precision = precision,
        adjacency = adjacency,
        node_info = S4Vectors::DataFrame(
            feature = features,
            omic_layer = c("transcript", "transcript", "metabolite")),
        celltype = celltype,
        metadata = list(n_donors = 20L))
}
mono <- make_network("Monocyte", -0.45)
networkCelltype(mono)
#> [1] "Monocyte"
networkAdjacency(mono)
#>    G1 G2 M1
#> G1  0  0  1
#> G2  0  0  0
#> M1  1  0  0
networkNodes(mono)
#> DataFrame with 3 rows and 2 columns
#>       feature  omic_layer
#>   <character> <character>
#> 1          G1  transcript
#> 2          G2  transcript
#> 3          M1  metabolite
metadata(mono)$n_donors
#> [1] 20

comparison <- run_comparative_network(mono, make_network("B cell", 0))
diffEdges(comparison)
#> DataFrame with 1 row and 7 columns
#>        node_i      node_j     layer_i     layer_j cross_omic        category
#>   <character> <character> <character> <character>  <logical>     <character>
#> 1          G1          M1  transcript  metabolite       TRUE unique_Monocyte
#>    max_diff
#>   <numeric>
#> 1      0.45
celltypeSpecificEdges(comparison)
#> $Monocyte
#> DataFrame with 1 row and 7 columns
#>        node_i      node_j     layer_i     layer_j cross_omic        category
#>   <character> <character> <character> <character>  <logical>     <character>
#> 1          G1          M1  transcript  metabolite       TRUE unique_Monocyte
#>    max_diff
#>   <numeric>
#> 1      0.45
#> 
comparisonSummary(comparison)$jaccard
#> $`Monocyte_vs_B cell`
#> [1] 0
#> 
```
