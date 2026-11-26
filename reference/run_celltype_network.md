# Build cell-type-resolved cross-omic network

Given a SCEE with pseudobulk transcriptomics and donor-level
metabolomics, constructs a block-calibrated conditional independence
network for a specified cell type.

## Usage

``` r
run_celltype_network(
  x,
  metabolites = NULL,
  celltype,
  celltype_col = "cell_type",
  exposure = NULL,
  covariates = NULL,
  method = c("coglasso", "block_glasso"),
  stability = TRUE,
  nlambda_w = 15L,
  nlambda_b = 15L,
  subsample_ratio = 0.8,
  rep_num = 20L,
  min_cells = 10L,
  top_var_genes = NULL,
  lambda_w = 0.3,
  lambda_b = 0.5,
  precomputed_network = NULL,
  transform = c("vst", "log_cpm"),
  genes = NULL
)
```

## Arguments

- x:

  [`SingleCellExposomeExperiment`](SingleCellExposomeExperiment-class.md).
  Each cell's donor is read from the `colData` column named by
  [`sampleColumn`](SCEE-accessors.md)`(x)`.

- metabolites:

  Numeric matrix; rows = donors, cols = metabolites. Row names must be
  donor IDs. A matrix without columns, or a missing or infinite value
  for an analysed donor, is an error, and so is a column with the name
  of another column or of a gene that becomes a node: two nodes cannot
  have the same name. The matrix is needed to estimate a network; the
  default, `NULL`, is for a call with `precomputed_network`, which does
  not use it.

- celltype:

  Character; which cell type to build the network for.

- celltype_col:

  Character; column in `colData` identifying cell types. Default
  `"cell_type"`.

- exposure:

  Character or NULL; exposure variable to condition on (regressed out),
  a column of `exposureData(x)`; another name is an error. Default NULL.

- covariates:

  Character vector or NULL; adjustment covariates, columns of
  `exposureData(x)`; another name is an error. Donors with a missing
  value of `exposure` or of a covariate are left out, with a message.
  Default NULL.

- method:

  Character; `"coglasso"` (default) or `"block_glasso"`. The coglasso
  method uses collaborative graphical lasso with block-calibrated
  penalties (Albanese et al. 2026) and needs the coglasso package,
  version 1.1.0 or later.

- stability:

  Logical; for `method = "coglasso"`, whether to return the edge-level
  selection frequencies from XStARS in `stability_scores`. The coglasso
  penalties are chosen by XStARS either way. `block_glasso` uses the
  fixed penalties `lambda_w` and `lambda_b`, performs no stability
  selection and ignores this argument. Default TRUE.

- nlambda_w:

  Integer; within-omic lambda grid size. Default 15.

- nlambda_b:

  Integer; between-omic lambda grid size. Default 15.

- subsample_ratio:

  Numeric strictly between 0 and 1; proportion of donors drawn in each
  XStARS subsample, passed to
  [`coglasso::xstars()`](https://drquestion.github.io/coglasso/reference/xstars.html)
  as `stars_subsample_ratio`. Default 0.8.

- rep_num:

  Integer, at least 2; number of XStARS subsamples. Default 20, as in
  [`coglasso::xstars()`](https://drquestion.github.io/coglasso/reference/xstars.html).

- min_cells:

  Integer; minimum cells per pseudobulk sample. Default 10.

- top_var_genes:

  Integer or NULL; restrict to top N most variable genes for
  dimensionality reduction. Default NULL (all genes). Recommended:
  50-200.

- lambda_w:

  Numeric; graphical-lasso penalty for pairs of features in the same
  omic layer, on the correlation scale, for the block_glasso method.
  Default 0.3.

- lambda_b:

  Numeric; graphical-lasso penalty for pairs of features in different
  omic layers, on the correlation scale, for the block_glasso method.
  Default 0.5.

- precomputed_network:

  Optional named list holding a pre-estimated network to wrap into a
  CelltypeNetworkResult without re-estimation, with elements `adjacency`
  (symmetric binary matrix, required), `precision` (numeric matrix,
  optional), `stability` (numeric matrix, optional), `feature_names`
  (character vector of unique, non-missing and non-empty names,
  optional; taken from the column names of `adjacency` when absent), and
  either `block_sizes` (a named vector giving the size of each omic
  block) or `feature_blocks` (a character vector labelling each
  feature), both optional. The block labels become the `omic_layer`
  column of the node table, so use `"transcript"` and `"metabolite"` for
  the summaries and plots to recognise the layers. A logical adjacency
  is stored as 0/1 and its diagonal is set to zero. When provided,
  nothing is estimated, and of the other arguments only `celltype` is
  used. `x` must still be a `SingleCellExposomeExperiment`.

- transform:

  Character; the transform of the pseudobulk counts. `"vst"` (default)
  is the variance-stabilising transformation of
  [`DESeq2::vst()`](https://rdrr.io/pkg/DESeq2/man/vst.html), blind to
  the design, and needs DESeq2 and at least 1000 genes with a mean
  normalised count above 5; `"log_cpm"` is the base-2 logarithm of the
  counts per million plus one. The call stops when `"vst"` cannot be
  computed; it does not switch to the other transform. One choice stays
  with DESeq2: when the parametric fit of the dispersion trend fails,
  `vst()` uses a local fit, and a message says so.

- genes:

  Character vector or `NULL`; row names of `x`, the genes that become
  the transcript nodes of the network, in the order given. The transform
  is computed from all genes of `x` before these are taken, so the
  values of a gene do not depend on which genes are kept, as they would
  after subsetting `x`. When `top_var_genes` is smaller than their
  number, the most variable of these genes are kept, in the order of
  their variances. Networks of several cell types built with the same
  `genes` have the same transcripts, which
  [`run_comparative_network`](run_comparative_network.md) needs, unless
  one of the genes is constant in a cell type and is left out there (see
  Details). A name that is not a row name of `x`, or that is the name of
  more than one row, is an error. Default `NULL` (all genes).

## Value

A [`CelltypeNetworkResult`](CelltypeNetworkResult-class.md) object.

## Details

The workflow proceeds as follows:

1.  Pseudobulk aggregation for the specified cell type, reusing existing
    exposomeSC infrastructure.

2.  Transformation of the pseudobulk counts, chosen by `transform`:
    DESeq2's variance-stabilising transformation or log-CPM. The
    transform is recorded in `metadata(result)$transform`.

3.  Constant features: a transcript or metabolite that takes one value
    in every donor (for example a gene with no counts in the cell type)
    is left out, with a message that names up to five transcripts and
    five metabolites; all are recorded in
    `metadata(result)$constant_features`; the call stops when every
    transcript or every metabolite is constant. Networks of different
    cell types can therefore differ in their features, which
    [`run_comparative_network`](run_comparative_network.md) does not
    accept. For a comparison, name the genes with `genes` and leave a
    gene that is constant in one of the cell types out of `genes` for
    all of them. `top_var_genes` on its own selects different genes in
    different cell types.

4.  Adjustment: if `exposure` or `covariates` are specified, regresses
    them out from both transcriptomic and metabolomic matrices. The
    resulting network captures the exposure-adjusted conditional
    dependency structure.

5.  Column-binding of the transcript values and the metabolite values
    (the residuals of both, where something was regressed out), each
    column scaled to unit variance.

6.  Block-calibrated graphical lasso estimation. With `coglasso`:
    collaborative graphical lasso with separate within/between-omic
    penalties and XStARS stability selection. With `block_glasso`: the
    graphical lasso
    ([`glasso::glasso()`](https://rdrr.io/pkg/glasso/man/glasso.html))
    on the correlation matrix, with the fixed penalties `lambda_w`
    within and `lambda_b` between omic layers, the diagonal unpenalised
    and no stability selection; it needs the glasso package.

7.  Returns a `CelltypeNetworkResult` object.

## References

Albanese A, Kohlen W, Behrouzi P (2026). Multi-omics network
reconstruction with collaborative graphical lasso. *Bioinformatics*
42(7):btag477.
[doi:10.1093/bioinformatics/btag477](https://doi.org/10.1093/bioinformatics/btag477)

Friedman J, Hastie T, Tibshirani R (2008). Sparse inverse covariance
estimation with the graphical lasso. *Biostatistics* 9(3):432–441.
[doi:10.1093/biostatistics/kxm045](https://doi.org/10.1093/biostatistics/kxm045)

Liu H, Roeder K, Wasserman L (2010). Stability approach to
regularization selection (StARS) for high dimensional graphical models.
*Advances in Neural Information Processing Systems*, 23.

## Examples

``` r
donor_ids <- paste0("D", seq_len(8L))
cell_donor <- rep(donor_ids, each = 10L)
cell_ids <- paste0("cell", seq_along(cell_donor))
counts <- outer(
    seq_len(4L),
    seq_along(cell_ids),
    function(gene, cell) 2L + ((5L * gene + 3L * cell) %% 9L)
)
storage.mode(counts) <- "integer"
dimnames(counts) <- list(paste0("G", seq_len(4L)), cell_ids)
sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = counts),
    colData = S4Vectors::DataFrame(
        cell_id = cell_ids,
        donor_id = cell_donor,
        cell_type = "Monocyte"
    )
)
exposure <- matrix(seq_len(8L), ncol = 1L,
    dimnames = list(donor_ids, "exposure"))
scee <- build_scee(sce, exposure, sample_col = "donor_id")
metabolites <- outer(
    seq_len(8L),
    seq_len(2L),
    function(donor, feature) (donor + feature)^2 / 10
)
dimnames(metabolites) <- list(donor_ids, c("M1", "M2"))
if (requireNamespace("glasso", quietly = TRUE)) {
    network <- run_celltype_network(
        scee,
        metabolites,
        celltype = "Monocyte",
        method = "block_glasso",
        stability = FALSE,
        min_cells = 5L,
        top_var_genes = 3L,
        transform = "log_cpm"
    )
    network
}
#> [exposomeSC] Building block_glasso network for 'Monocyte': n=8 donors, p=5 features (3 transcripts + 2 metabolites)
#> CelltypeNetworkResult for 'Monocyte'
#>   Method: block_glasso
#>   Nodes: 5 (3 transcripts, 2 metabolites)
#>   Edges: 4 (0 cross-omic)
#>   Donors: 8
```
