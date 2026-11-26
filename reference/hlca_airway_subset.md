# Airway epithelial cells of smokers and never smokers from the Human Lung Cell Atlas

A
[`SingleCellExposomeExperiment`](SingleCellExposomeExperiment-class.md)
with the raw counts of 352 tracheal epithelial cells from 12 donors of
the Human Lung Cell Atlas (HLCA) core: six heavy smokers and six never
smokers, with up to 10 basal, 10 secretory and 10 multiciliated-lineage
cells per donor. It is a small real data set for the vignettes and the
tests. It shows how the functions are used on real counts; with 12
donors and at most 10 cells per donor and cell type it is far too small
to support conclusions about smoking.

## Format

An RDS file, `inst/extdata/hlca_airway_subset.rds`, read with
`readRDS(system.file("extdata", "hlca_airway_subset.rds", package = "exposomeSC"))`.

## Source

The CELLxGENE H5AD of the HLCA core, dataset version
`688185ad-11c2-4172-a53a-f4f1f4076860`
(<https://datasets.cellxgene.cziscience.com/688185ad-11c2-4172-a53a-f4f1f4076860.h5ad>,
5,873,612,847 bytes, SHA-256
`1cbdde1e513a31bdb6a3f70296e2d0612760d2af4032e239c94596ea87b44007`),
distributed under the Creative Commons Attribution 4.0 International
licence (CC BY 4.0, <https://creativecommons.org/licenses/by/4.0/>). The
subset was made by `inst/scripts/make_hlca_airway_subset.R`, which
checks the file and carries out the selection described above. Cells and
genes were selected; the counts themselves are unchanged.

## Details

- Cells:

  Study `"Seibold_2020"` of the HLCA core (Goldfarbmuren et al. 2020):
  the heavy smokers T85, T90, T101, T120, T154 and T167 (smokers with at
  least 15 pack-years) and the never smokers T126, T137, T153, T164,
  T165 and T166, the two groups that the authors compare, as classified
  in their public metadata. The other two donors of the study in the
  atlas, a light smoker and a child who never smoked, are left out.
  Cells annotated `"Basal"`, `"Secretory"` or `"Multiciliated lineage"`
  at HLCA annotation level 3 were eligible, and up to 10 per donor and
  cell type were drawn at random with a fixed seed, without looking at
  expression; donor T154 has only two multiciliated-lineage cells. The
  donors and the cell types were fixed before any result on cells of
  these donors was seen. The number of cells per donor and cell type was
  set later, by the size limit of a Bioconductor source package: before
  this subset was drawn, but after the smoking model had been run on
  larger sets of cells of the same donors.

- Counts:

  `assay(x, "counts")`: the raw counts of the H5AD (`/raw/X`), 20,111
  genes by 352 cells, as a sparse matrix. Genes without a count in these
  cells are left out, so the column sums are those of the complete
  transcriptome. The rows are Ensembl gene identifiers.
  `rowData(x)$symbol` holds the gene names of the atlas
  (`feature_name`): for 2,301 genes that name is the Ensembl identifier
  again, and three names (MATR3, RAET1E-AS1 and SPATA13) are each on two
  rows.

- Cell annotations:

  `colData(x)`: `donor_id` (the donor code of the study), `cell_type`,
  `assay` (10x 3' v2 or v3) and `hlca_donor_id` (the donor identifier of
  the HLCA).

- Donor variables:

  `exposureData(x)`, one row per donor: `smoking` (1 = heavy smoker, 0 =
  never smoker), `chemistry_v3` (1 = 10x 3' v3, 0 = 10x 3' v2), `age`
  (years) and `sex_male` (1 = male, 0 = female).

- Provenance:

  `metadata(x)$source`: the URL, dataset version, size and SHA-256 of
  the H5AD, the licence, the references, the selection rule and, in
  `eligible_cells`, the number of eligible cells of each donor and cell
  type in the atlas (29,217 in all), from which the cells of the subset
  were drawn.

## References

Sikkema L, et al. (2023). An integrated cell atlas of the lung in health
and disease. *Nature Medicine*, 29, 1563-1577.
[doi:10.1038/s41591-023-02327-2](https://doi.org/10.1038/s41591-023-02327-2)

Goldfarbmuren KC, et al. (2020). Dissecting the cellular specificity of
smoking effects and reconstructing lineages in the human airway
epithelium. *Nature Communications*, 11, 2485.
[doi:10.1038/s41467-020-16239-z](https://doi.org/10.1038/s41467-020-16239-z)

## Examples

``` r
x <- readRDS(system.file("extdata", "hlca_airway_subset.rds",
    package = "exposomeSC"))
x
#> class: SingleCellExposomeExperiment 
#> dim: 20111 352 
#> metadata(1): source
#> assays(1): counts
#> rownames(20111): ENSG00000000003 ENSG00000000419 ... ENSG00000283103
#>   ENSG00000283125
#> rowData names(1): symbol
#> colnames(352): TCTTCCTTCAGAGCAG_T101 CTACTATCAACGGGTA_T101 ...
#>   CCTCTAGAGCAAGCCA_T90 CCTCATGGTCTACATG_T90
#> colData names(4): donor_id cell_type assay hlca_donor_id
#> reducedDimNames(0):
#> mainExpName: NULL
#> altExpNames(0):
#> exposureData(12 samples x 4 exposures): smoking, chemistry_v3, age, sex_male
#> sampleColumn: donor_id
table(x$donor_id, x$cell_type)
#>       
#>        Basal Multiciliated lineage Secretory
#>   T101    10                    10        10
#>   T120    10                    10        10
#>   T126    10                    10        10
#>   T137    10                    10        10
#>   T153    10                    10        10
#>   T154    10                     2        10
#>   T164    10                    10        10
#>   T165    10                    10        10
#>   T166    10                    10        10
#>   T167    10                    10        10
#>   T85     10                    10        10
#>   T90     10                    10        10
exposureData(x)
#>      smoking chemistry_v3 age sex_male
#> T101       1            1  55        1
#> T120       1            0  57        0
#> T126       0            0  35        0
#> T137       0            0  27        0
#> T153       0            1  38        1
#> T154       1            1  61        1
#> T164       0            1  66        1
#> T165       0            1  64        1
#> T166       0            1  68        0
#> T167       1            1  66        0
#> T85        1            0  59        1
#> T90        1            1  44        1
S4Vectors::metadata(x)$source$licence
#> [1] "CC BY 4.0"
```
