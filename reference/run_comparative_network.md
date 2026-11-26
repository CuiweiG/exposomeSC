# Compare networks across cell types

Classifies the union of edges from two or more cell-type-specific
networks estimated on the same features as shared, unique to one cell
type or differential, and reports pairwise Jaccard overlap. The
comparison is descriptive: no edge is tested and no p-value is reported.

## Usage

``` r
run_comparative_network(...)
```

## Arguments

- ...:

  Two or more [`CelltypeNetworkResult`](CelltypeNetworkResult-class.md)
  objects. All must share the same feature set (same genes and
  metabolites).

## Value

A [`NetworkComparison`](NetworkComparison-class.md) object. Its edge
tables contain `max_diff`, the range of the edge's partial correlation
across networks.

## Details

Partial correlations are computed from each precision matrix \\\Theta\\
as \\-\Theta\_{ij} / \sqrt{\Theta\_{ii}\Theta\_{jj}}\\.

This is the single-cell analogue of the temporal edge comparison in
Cheng et al. (ES&T 2024).

## References

Cheng SL et al. (2024). Multiomic signatures of traffic-related air
pollution in London reveal potential short-term perturbations in gut
microbiome-related pathways. *Environ Sci Technol* 58:8771-8782.
[doi:10.1021/acs.est.3c09148](https://doi.org/10.1021/acs.est.3c09148)

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
            omic_layer = c("transcript", "transcript", "metabolite")
        ),
        celltype = celltype,
        metadata = list(n_donors = 20L)
    )
}
monocyte_network <- make_network("Monocyte", -0.45)
b_cell_network <- make_network("B cell", -0.15)
comparison <- run_comparative_network(monocyte_network, b_cell_network)
comparison
#> NetworkComparison
#>   Cell types compared: Monocyte, B cell
#>   Shared edges: 1
#>   Differential edges: 0
```
