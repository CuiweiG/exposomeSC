# Simulate cross-omic network data with known structure

Generates donor-level metabolite data and cell-type-specific transcript
pseudobulk counts whose latent log-expression follows a known Gaussian
graphical model, with cell-type heterogeneity and optional
exposure-dependent composition. Intended for benchmarking edge recovery
and false-positive control.

## Usage

``` r
simulate_crossomic_network(
  n_donors = 100L,
  n_celltypes = 4L,
  n_transcripts = 30L,
  n_metabolites = 10L,
  n_cells_per_donor = 200L,
  edge_density = 0.1,
  shared_edge_frac = 0.6,
  cross_omic_frac = 0.2,
  exposure_effect = 0.5,
  composition_confounding = TRUE,
  seed = NULL
)
```

## Arguments

- n_donors:

  Integer; number of donors. Default 100.

- n_celltypes:

  Integer; number of cell types. Default 4.

- n_transcripts:

  Integer; number of transcript features. Default 30.

- n_metabolites:

  Integer; number of metabolite features. Default 10.

- n_cells_per_donor:

  Integer; cells per donor, allocated to cell types by the composition.
  Pseudobulk depth scales with the number of cells. Default 200.

- edge_density:

  Numeric; proportion of possible edges present in the base network.
  Default 0.1.

- shared_edge_frac:

  Numeric; fraction of base-network edges kept by cell types 1 and 2,
  each of which also gains the same number of new edges. Default 0.6.

- cross_omic_frac:

  Numeric; fraction of edges that are cross-omic. Default 0.2.

- exposure_effect:

  Numeric; effect of the exposure on the mean of affected features and
  on composition. Default 0.5.

- composition_confounding:

  Logical; make the proportions of cell types 1 and 2 depend on the
  exposure. Default TRUE.

- seed:

  Integer or NULL; random seed. The caller's random number generator
  state is restored on return. Default NULL.

## Value

List with:

- pseudobulk:

  List of count matrices (donors x transcripts), one per cell type.

- metabolites:

  Metabolite matrix (donors x metabolites), shared by all cell types.

- exposure:

  Named numeric vector of exposure values.

- composition:

  Donor-by-cell-type matrix of proportions.

- n_cells:

  Donor-by-cell-type matrix of simulated cell numbers.

- true_networks:

  List of precision matrices, one per cell type, of the latent
  transcript log-expression and the metabolites conditional on the
  exposure.

- true_adjacency:

  List of the corresponding adjacency matrices.

- celltype_names:

  Character vector of cell type names.

- node_info:

  DataFrame of feature metadata.

- params:

  List of simulation parameters.

## Details

Simulation design:

1.  A sparse positive-definite base precision matrix \\\Theta\\ is
    generated over transcripts and metabolites with within-omic and
    cross-omic edges.

2.  Cell types 1 and 2 keep a fraction `shared_edge_frac` of the base
    edges and gain as many new edges; cell type 3 has an independently
    generated precision matrix; cell types 4 and above have no
    transcript-transcript or cross-omic edges.

3.  Metabolites are drawn once per donor from a multivariate normal
    distribution whose precision is the metabolite block of the base
    matrix, \\P_m\\.

4.  For each cell type, latent transcript log-expression is drawn from
    its conditional distribution given the metabolites,
    \\N(-\Theta\_{tt}^{-1}\Theta\_{tm} m, \Theta\_{tt}^{-1})\\. The
    joint precision of transcripts and metabolites is then
    \\\Theta\_{tt}\\, \\\Theta\_{tm}\\ and \\P_m +
    \Theta\_{mt}\Theta\_{tt}^{-1}\Theta\_{tm}\\; it is returned in
    `true_networks`, so its metabolite block can contain edges induced
    by the cross-omic edges.

5.  The exposure shifts the mean of the first 30% of metabolites and, in
    cell types 1 and 2, of the first 20% of transcripts; it does not
    change any network.

6.  Pseudobulk counts are Poisson with mean equal to the number of cells
    times \\\exp(\text{latent} + 1)\\, which makes them overdispersed
    relative to a Poisson with a fixed mean.

7.  If `composition_confounding = TRUE`, the proportions of cell types 1
    and 2 increase with the exposure, so pseudobulk depth varies with
    the exposure.

## Examples

``` r
sim <- simulate_crossomic_network(
    n_donors = 50, n_transcripts = 20,
    n_metabolites = 5, seed = 42)
names(sim)
#>  [1] "pseudobulk"     "metabolites"    "exposure"       "composition"   
#>  [5] "n_cells"        "true_networks"  "true_adjacency" "celltype_names"
#>  [9] "node_info"      "params"        
sim$true_adjacency[[1]][1:5, 1:5]
#>        gene_1 gene_2 gene_3 gene_4 gene_5
#> gene_1      0      0      0      0      0
#> gene_2      0      0      0      0      0
#> gene_3      0      0      0      0      0
#> gene_4      0      0      0      0      0
#> gene_5      0      0      0      0      0
```
