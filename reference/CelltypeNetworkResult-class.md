# CelltypeNetworkResult: cell-type-resolved cross-omic network

Stores the estimated precision matrix, adjacency matrix, stability
scores, and node metadata from
[`run_celltype_network`](run_celltype_network.md) and
[`run_exposure_network`](run_exposure_network.md). Use the accessors in
[`network-accessors`](network-accessors.md) to extract them, and
[`CelltypeNetworkResult()`](CelltypeNetworkResult.md) to wrap a network
estimated elsewhere.

## Usage

``` r
# S4 method for class 'CelltypeNetworkResult'
show(object)
```

## Arguments

- object:

  A `CelltypeNetworkResult` object.

## Value

The `show` method prints a summary and returns `NULL` invisibly.

## Slots

- `precision_matrix`:

  Numeric matrix. Estimated block precision (inverse covariance) matrix.

- `adjacency_matrix`:

  Numeric matrix. Symmetric binary edge indicators (1 = edge present, 0
  = absent) with the dimensions of `precision_matrix` and a zero
  diagonal.

- `stability_scores`:

  Numeric matrix. Edge-level selection frequencies (0-1) across the
  XStARS subsamples at the selected penalties, for coglasso networks
  built with `stability = TRUE`. Empty otherwise.

- `node_info`:

  `DataFrame`. Maps each node (row/col index) to its omic layer
  ("transcript" or "metabolite"), feature name, and optional
  module/community assignment. It has one row per node, in matrix order;
  when the matrices have row or column names, these must equal its
  `feature` column.

- `celltype`:

  Character. The cell type this network represents.

- `method`:

  Character. Estimation method: `"coglasso"` or `"block_glasso"` from
  [`run_celltype_network()`](run_celltype_network.md),
  `"exposure_restricted"` from
  [`run_exposure_network()`](run_exposure_network.md), or
  `"precomputed"` for a network supplied by the user.

- `metadata`:

  List. Additional parameters: n_donors, n_features, lambda values,
  exposure, covariates.

## See also

[`network-accessors`](network-accessors.md),
[`CelltypeNetworkResult`](CelltypeNetworkResult.md)

## Examples

``` r
showClass("CelltypeNetworkResult")
#> Class "CelltypeNetworkResult" [package "exposomeSC"]
#> 
#> Slots:
#>                                                                           
#> Name:  precision_matrix adjacency_matrix stability_scores        node_info
#> Class:           matrix           matrix           matrix        DataFrame
#>                                                          
#> Name:          celltype           method         metadata
#> Class:        character        character             list
#> 
#> Extends: "CelltypeNetwork"
```
