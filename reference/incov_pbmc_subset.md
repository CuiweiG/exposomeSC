# Blood cells, surface proteins and plasma measurements of patients with COVID-19 from the INCOV cohort

A
[`SingleCellExposomeExperiment`](SingleCellExposomeExperiment-class.md)
with the raw counts of 1,032 blood mononuclear cells from 43 patients of
the INCOV cohort (Su et al. 2020) at their first blood draw: 14
monocytes and 10 T cells per patient, with the counts of 192 antibodies
for the same cells, four of them isotype controls, 13 plasma compounds
of exogenous origin among the donor variables and 182 plasma metabolites
of the same draw. It is a small real data set for the functions that
need more than transcript counts: surface proteins, exposure
measurements with missing values, and a second omic layer measured once
per donor. It shows how these functions are used on real data. With at
most 14 cells per donor and cell type it is far too small to support
conclusions, and all its donors had an acute viral infection.

## Format

An RDS file, `inst/extdata/incov_pbmc_subset.rds`, read with
`readRDS(system.file("extdata", "incov_pbmc_subset.rds", package = "exposomeSC"))`.

## Source

The clinical records and the plasma measurements are Table S1 of
Mendeley Data,
[doi:10.17632/tzydswhhb5.5](https://doi.org/10.17632/tzydswhhb5.5)
(version 5, published on 4 November 2020), distributed under the
Creative Commons Attribution 4.0 International licence (CC BY 4.0,
<https://creativecommons.org/licenses/by/4.0/>). The plasma values were
selected, not changed; smoking and sex were recoded as described above.

The single-cell values of genes and surface proteins are from
ArrayExpress accession E-MTAB-9357
(<https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-9357>),
released on 29 October 2020 and downloaded on 19 November 2026. On 21
November 2026 the authors who deposited the accession gave permission in
these words: "The single-cell values from E-MTAB-9357 included in
exposomeSC may be used and redistributed under CC BY 4.0, with citation
of Su et al. (2020)." The object holds a subset of the values that the
authors released, transformed as described above.

The object as a whole is therefore distributed under CC BY 4.0. Work
that uses it should cite Su et al. (2020). The permission does not imply
that the authors endorse exposomeSC.

The article itself is published under CC BY-NC-ND 4.0. The values in
this object come from the data that the authors deposited in Mendeley
Data and ArrayExpress, not from the article.

The object was made by `inst/scripts/make_incov_pbmc_subset.R`, which
checks the size and the SHA-256 of every input file and carries out the
steps described above.

## Details

- Donors:

  Patients of the INCOV cohort who, at their first blood draw, had
  single cells, plasma metabolomics and plasma proteomics measured, a
  recorded history of cigarette smoking, and were not in hospital (WHO
  ordinal scale 1 or 2). 43 patients meet these conditions: 37 women and
  6 men aged 19 to 81, 16 of them current or former smokers. The rule
  was fixed from the clinical sheet and the list of files before any
  single-cell value of these donors was read.

- Counts:

  `assay(x, "counts")`: 15,406 genes by 1,032 cells, as a sparse matrix;
  the rows are gene symbols. The archive holds the natural logarithm of
  one plus the counts per million of each cell, from which the counts
  were recovered: every value of a cell is a whole multiple of the value
  of one count, and that value times the total of the cell is one
  million, which the script checks for every cell. Genes without a count
  in these cells are left out, so the column sums are those of all
  24,966 genes of the archive files. The authors had kept cells with at
  least 200 genes, fewer than 10,000 counts, fewer than 2,500 genes and
  under 10 per cent mitochondrial counts.

- Surface proteins:

  `altExp(x, "CITE")`: the counts of 192 antibodies, four of them
  isotype controls, for the same cells, recovered in the same way.

- Cell types:

  `x$cell_type`, assigned for this data set; the archive holds no cell
  types. Each sample has two gates, one for the protein CD3 and one for
  the sum of the proteins CD14, CD11c, CD33 and CD64. A gate is the
  midpoint of the two centres that a split of the cells into two groups
  by the logarithm of one plus the count gives. A `"T cell"` has CD3
  above its gate, the sum not above its gate, and fewer than 10 counts
  of the genes *LYZ*, *S100A8* and *S100A9* together. A `"Monocyte"` has
  the sum above its gate, CD3 not above its gate, at least 5 counts of
  those three genes and fewer than 3 counts of *CD3D*, *CD3E* and *TRAC*
  together. Of the 63,865 cells of the 43 samples, 9,889 are monocytes
  and 38,213 T cells by this rule. In one sample, INCOV132, the two
  groups of the sum lie close together (its gate is 8.8 counts, the
  lowest; the others are 12 or more): 121 of its 328 cells lie above
  both gates and are neither kind, against fewer than a tenth of the
  cells of any other sample, and of the 43 cells that the gates alone
  would make monocytes the conditions on the genes keep 16. The rule was
  chosen by how well alternatives agreed with these genes in 42 of the
  samples and was applied to the 43rd as it stood. The five proteins and
  six genes of the rule are constrained by it and must not be read as
  results.

- Cells:

  14 monocytes and 10 T cells per donor, drawn at random with a fixed
  seed from the cells of each kind. The two numbers were set by the size
  limit of a Bioconductor source package.

- Cell annotations:

  `colData(x)`: `donor_id` (the patient identifier of the study),
  `cell_type`, `sample_id` (the name of the sample in the archive, the
  patient identifier followed by `-BL`) and `archive_cell` (the name of
  the cell in the files of the archive). The column names of the object
  are the patient identifier and the barcode of the cell. The rows are
  named by gene symbol, which `rowData(x)$symbol` repeats.

- Donor variables:

  `exposureData(x)`, one row per donor: `ever_smoker` (1 = current or
  former smoker, 0 = never smoker), `age` (years), `sex_male` (1 = male,
  0 = female), and the plasma levels of 13 compounds of exogenous
  origin: `cotinine`; `caffeine`, `paraxanthine` and `theobromine`;
  `pfos` and `pfoa` (perfluorooctanesulfonate and perfluorooctanoate);
  `methylparaben_sulfate` and `propylparaben_sulfate`;
  `hydroxychlorothalonil`; `naphthol_sulfate`; `ethyl_glucuronide`;
  `saccharin` and `acesulfame`. The levels are relative abundances,
  scaled by the authors so that the median of a compound within a batch
  is 1, and `NA` where the table has no value for the sample. Ten of the
  13 compounds lack a value for at least one donor: cotinine has a value
  for only 5 donors, propylparaben sulfate for 17, acesulfame for 19 and
  saccharin for 27. The table does not say why a value is missing. In
  measurements of this kind the usual reason is that the compound was
  not detected in the sample, but it is not the only one: the table also
  lacks a value of glucose for one of these donors. `exposureInfo(x)`
  gives the name of each compound in the source table (`source_name`),
  its family, and `lowest_reported`, the smallest value of the compound
  among the 387 samples of the table. The table reports no limits of
  detection. For the metabolites that they analysed, those with fewer
  than 20 per cent of the values missing, the authors put this smallest
  value in place of a missing one; six of the 13 compounds here have
  more missing values than that. The 13 compounds were chosen by their
  kind, from their names and the numbers of samples with a value.

- Plasma metabolites:

  `metadata(x)$plasma_metabolites`, 43 donors by 182 metabolites: the
  metabolites that the source table assigns to a super-pathway other
  than xenobiotics and reports for all 43 donors, on the same relative
  scale. Their pathways are in `metadata(x)$metabolite_info`.

- Provenance:

  `metadata(x)$source`: the reference, the sources and their licences,
  the permission for the single-cell values (`single_cells$permission`),
  the names, sizes and SHA-256 of the 88 input files, the selection
  rule, the gates of each sample (`gates`) and the numbers of monocytes,
  T cells and other cells of each donor (`gated_cells`).

The rule for the donors, the rule for the cell types, the numbers of
cells and the seed were fixed before any association of a donor variable
with expression, proteins or metabolites was computed on these donors.
The object was first made from 42 of the 43 donors, because the first
transfer of the expression file of one was incomplete, and the vignettes
were first written on that version. It was then made again from all 43
by the script without that exclusion; the cells of the 14 donors that
follow that donor in the order of the draw are a new random draw.

## References

Su Y, Chen D, Yuan D, et al. (2020). Multi-omics resolves a sharp
disease-state shift between mild and moderate COVID-19. *Cell*, 183(6),
1479-1495.e20.
[doi:10.1016/j.cell.2020.10.037](https://doi.org/10.1016/j.cell.2020.10.037)

## Examples

``` r
x <- readRDS(system.file("extdata", "incov_pbmc_subset.rds",
    package = "exposomeSC"))
x
#> class: SingleCellExposomeExperiment 
#> dim: 15406 1032 
#> metadata(3): source plasma_metabolites metabolite_info
#> assays(1): counts
#> rownames(15406): A1BG A1BG-AS1 ... ZYX ZZEF1
#> rowData names(1): symbol
#> colnames(1032): INCOV039_AAGACCTAGGGTCTCC INCOV039_ACTGTCCGTGAGGCTA ...
#>   INCOV135_GGCGTGTGTACCGTAT INCOV135_GTATCTTAGAAACGCC
#> colData names(4): donor_id cell_type sample_id archive_cell
#> reducedDimNames(0):
#> mainExpName: NULL
#> altExpNames(1): CITE
#> exposureData(43 samples x 16 exposures): ever_smoker, age, sex_male, cotinine, caffeine ...
#> sampleColumn: donor_id
table(x$cell_type)
#> 
#> Monocyte   T cell 
#>      602      430 
exposureData(x)[1:3, 1:6]
#>          ever_smoker age sex_male cotinine  caffeine paraxanthine
#> INCOV039           0  38        0       NA 9.8244059           NA
#> INCOV042           0  60        0       NA 0.9791921    0.3850829
#> INCOV056           1  61        0       NA 6.3010545    2.1963925
SingleCellExperiment::altExp(x, "CITE")
#> class: SummarizedExperiment 
#> dim: 192 1032 
#> metadata(0):
#> assays(1): counts
#> rownames(192): B7-H4 CD10 ... XCR1-1 c-Met
#> rowData names(0):
#> colnames(1032): INCOV039_AAGACCTAGGGTCTCC INCOV039_ACTGTCCGTGAGGCTA ...
#>   INCOV135_GGCGTGTGTACCGTAT INCOV135_GTATCTTAGAAACGCC
#> colData names(0):
dim(S4Vectors::metadata(x)$plasma_metabolites)
#> [1]  43 182
S4Vectors::metadata(x)$source$single_cells$licence
#> [1] "CC BY 4.0"
```
