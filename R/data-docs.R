#' Airway epithelial cells of smokers and never smokers from the Human Lung
#' Cell Atlas
#'
#' A \code{\linkS4class{SingleCellExposomeExperiment}} with the raw counts
#' of 352 tracheal epithelial cells from 12 donors of the Human Lung Cell
#' Atlas (HLCA) core: six heavy smokers and six never smokers, with up to 10
#' basal, 10 secretory and 10 multiciliated-lineage cells per donor. It is a
#' small real data set for the vignettes and the tests. It
#' shows how the functions are used on real counts; with 12 donors and at
#' most 10 cells per donor and cell type it is far too small to support
#' conclusions about smoking.
#'
#' @details
#' \describe{
#'   \item{Cells}{Study \code{"Seibold_2020"} of the HLCA core (Goldfarbmuren
#'   et al. 2020): the heavy smokers T85, T90, T101, T120, T154 and T167
#'   (smokers with at least 15 pack-years) and the never smokers T126,
#'   T137, T153, T164, T165 and T166, the two groups that the authors
#'   compare, as classified in their public metadata. The other two donors
#'   of the study in the atlas, a light smoker and a child who never
#'   smoked, are left out. Cells annotated \code{"Basal"},
#'   \code{"Secretory"} or \code{"Multiciliated lineage"} at HLCA
#'   annotation level 3 were eligible, and up to 10 per donor and cell type
#'   were drawn at random with a fixed seed, without looking at expression;
#'   donor T154 has only two multiciliated-lineage cells. The donors and
#'   the cell types were fixed before any result on cells of these donors
#'   was seen. The number of cells per donor and cell type was set later,
#'   by the size limit of a Bioconductor source package: before this subset
#'   was drawn, but after the smoking model had been run on larger sets of
#'   cells of the same donors.}
#'   \item{Counts}{\code{assay(x, "counts")}: the raw counts of the H5AD
#'   (\code{/raw/X}), 20,111 genes by 352 cells, as a sparse matrix. Genes
#'   without a count in these cells are left out, so the column sums are
#'   those of the complete transcriptome. The rows are Ensembl gene
#'   identifiers. \code{rowData(x)$symbol} holds the gene names of the
#'   atlas (\code{feature_name}): for 2,301 genes that name is the Ensembl
#'   identifier again, and three names (MATR3, RAET1E-AS1 and SPATA13) are
#'   each on two rows.}
#'   \item{Cell annotations}{\code{colData(x)}: \code{donor_id} (the donor
#'   code of the study), \code{cell_type}, \code{assay} (10x 3' v2 or v3)
#'   and \code{hlca_donor_id} (the donor identifier of the HLCA).}
#'   \item{Donor variables}{\code{exposureData(x)}, one row per donor:
#'   \code{smoking} (1 = heavy smoker, 0 = never smoker),
#'   \code{chemistry_v3} (1 = 10x 3' v3, 0 = 10x 3' v2), \code{age} (years)
#'   and \code{sex_male} (1 = male, 0 = female).}
#'   \item{Provenance}{\code{metadata(x)$source}: the URL, dataset version,
#'   size and SHA-256 of the H5AD, the licence, the references, the
#'   selection rule and, in \code{eligible_cells}, the number of eligible
#'   cells of each donor and cell type in the atlas (29,217 in all), from
#'   which the cells of the subset were drawn.}
#' }
#'
#' @format An RDS file, \code{inst/extdata/hlca_airway_subset.rds}, read with
#'   \code{readRDS(system.file("extdata", "hlca_airway_subset.rds",
#'   package = "exposomeSC"))}.
#'
#' @source The CELLxGENE H5AD of the HLCA core, dataset version
#'   \code{688185ad-11c2-4172-a53a-f4f1f4076860}
#'   (\url{https://datasets.cellxgene.cziscience.com/688185ad-11c2-4172-a53a-f4f1f4076860.h5ad},
#'   5,873,612,847 bytes, SHA-256
#'   \code{1cbdde1e513a31bdb6a3f70296e2d0612760d2af4032e239c94596ea87b44007}),
#'   distributed under the Creative Commons Attribution 4.0 International
#'   licence (CC BY 4.0, \url{https://creativecommons.org/licenses/by/4.0/}).
#'   The subset was made by \code{inst/scripts/make_hlca_airway_subset.R},
#'   which checks the file and carries out the selection described above.
#'   Cells and genes were selected; the counts themselves are unchanged.
#'
#' @references
#' Sikkema L, et al. (2023). An integrated cell atlas of the lung in health
#' and disease. \emph{Nature Medicine}, 29, 1563-1577.
#' \doi{10.1038/s41591-023-02327-2}
#'
#' Goldfarbmuren KC, et al. (2020). Dissecting the cellular specificity of
#' smoking effects and reconstructing lineages in the human airway
#' epithelium. \emph{Nature Communications}, 11, 2485.
#' \doi{10.1038/s41467-020-16239-z}
#'
#' @examples
#' x <- readRDS(system.file("extdata", "hlca_airway_subset.rds",
#'     package = "exposomeSC"))
#' x
#' table(x$donor_id, x$cell_type)
#' exposureData(x)
#' S4Vectors::metadata(x)$source$licence
#'
#' @name hlca_airway_subset
#' @docType data
#' @keywords datasets
NULL

#' Blood cells, surface proteins and plasma measurements of patients with
#' COVID-19 from the INCOV cohort
#'
#' A \code{\linkS4class{SingleCellExposomeExperiment}} with the raw counts
#' of 1,032 blood mononuclear cells from 43 patients of the INCOV cohort
#' (Su et al. 2020) at their first blood draw: 14 monocytes and 10 T cells
#' per patient, with the counts of 192 antibodies for the same cells, four
#' of them isotype controls, 13 plasma compounds of exogenous origin among
#' the donor variables and 182 plasma metabolites of the same draw. It is a
#' small real data set for
#' the functions that need more than transcript counts: surface proteins,
#' exposure measurements with missing values, and a second omic layer
#' measured once per donor. It shows how these functions are used on real
#' data. With at most 14 cells per donor and cell type it is far too small
#' to support conclusions, and all its donors had an acute viral infection.
#'
#' @details
#' \describe{
#'   \item{Donors}{Patients of the INCOV cohort who, at their first blood
#'   draw, had single cells, plasma metabolomics and plasma proteomics
#'   measured, a recorded history of cigarette smoking, and were not in
#'   hospital (WHO ordinal scale 1 or 2). 43 patients meet these
#'   conditions: 37 women and 6 men aged 19 to 81, 16 of them current or
#'   former smokers. The rule was fixed from the clinical sheet and the
#'   list of files before any single-cell value of these donors was read.}
#'   \item{Counts}{\code{assay(x, "counts")}: 15,406 genes by 1,032 cells,
#'   as a sparse matrix; the rows are gene symbols. The archive holds the
#'   natural logarithm of one plus the counts per million of each cell,
#'   from which the counts were recovered: every value of a cell is a whole
#'   multiple of the value of one count, and that value times the total of
#'   the cell is one million, which the script checks for every cell. Genes
#'   without a count in these cells are left out, so the column sums are
#'   those of all 24,966 genes of the archive files. The authors had kept
#'   cells with at least 200 genes, fewer than 10,000 counts, fewer than
#'   2,500 genes and under 10 per cent mitochondrial counts.}
#'   \item{Surface proteins}{\code{altExp(x, "CITE")}: the counts of 192
#'   antibodies, four of them isotype controls, for the same cells,
#'   recovered in the same way.}
#'   \item{Cell types}{\code{x$cell_type}, assigned for this data set; the
#'   archive holds no cell types. Each sample has two gates, one for the
#'   protein CD3 and one for the sum of the proteins CD14, CD11c, CD33 and
#'   CD64. A gate is the midpoint of the two centres that a split of the
#'   cells into two groups by the logarithm of one plus the count gives. A
#'   \code{"T cell"} has CD3 above its gate, the sum not above its gate,
#'   and fewer than 10 counts of the genes \emph{LYZ}, \emph{S100A8} and
#'   \emph{S100A9} together. A \code{"Monocyte"} has the sum above its
#'   gate, CD3 not above its gate, at least 5 counts of those three genes
#'   and fewer than 3 counts of \emph{CD3D}, \emph{CD3E} and \emph{TRAC}
#'   together. Of the 63,865 cells of the 43 samples, 9,889 are monocytes
#'   and 38,213 T cells by this rule. In one sample, INCOV132, the two
#'   groups of the sum lie close together (its gate is 8.8 counts, the
#'   lowest; the others are 12 or more): 121 of its 328 cells lie above
#'   both gates and are neither kind, against fewer than a tenth of the
#'   cells of any other sample, and of the 43 cells that the gates alone
#'   would make monocytes the conditions on the genes keep 16. The rule was
#'   chosen by how well
#'   alternatives agreed with these genes in 42 of the samples and was
#'   applied to the 43rd as it stood. The five proteins and six genes of
#'   the rule are constrained by it and must not be read as results.}
#'   \item{Cells}{14 monocytes and 10 T cells per donor, drawn at random
#'   with a fixed seed from the cells of each kind. The two numbers were
#'   set by the size limit of a Bioconductor source package.}
#'   \item{Cell annotations}{\code{colData(x)}: \code{donor_id} (the
#'   patient identifier of the study), \code{cell_type}, \code{sample_id}
#'   (the name of the sample in the archive, the patient identifier
#'   followed by \code{-BL}) and \code{archive_cell} (the name of the cell
#'   in the files of the archive). The column names of the object are the
#'   patient identifier and the barcode of the cell. The rows are named by
#'   gene symbol, which \code{rowData(x)$symbol} repeats.}
#'   \item{Donor variables}{\code{exposureData(x)}, one row per donor:
#'   \code{ever_smoker} (1 = current or former smoker, 0 = never smoker),
#'   \code{age} (years), \code{sex_male} (1 = male, 0 = female), and the
#'   plasma levels of 13 compounds of exogenous origin: \code{cotinine};
#'   \code{caffeine}, \code{paraxanthine} and \code{theobromine};
#'   \code{pfos} and \code{pfoa} (perfluorooctanesulfonate and
#'   perfluorooctanoate); \code{methylparaben_sulfate} and
#'   \code{propylparaben_sulfate}; \code{hydroxychlorothalonil};
#'   \code{naphthol_sulfate}; \code{ethyl_glucuronide}; \code{saccharin}
#'   and \code{acesulfame}. The levels are relative abundances, scaled by
#'   the authors so that the median of a compound within a batch is 1, and
#'   \code{NA} where the table has no value for the sample. Ten of the 13
#'   compounds lack a value for at least one donor: cotinine has a value
#'   for only 5 donors, propylparaben sulfate for 17, acesulfame for 19 and
#'   saccharin for 27. The table does not say why a value is missing. In
#'   measurements of this kind the usual reason is that the compound was
#'   not detected in the sample, but it is not the only one: the table
#'   also lacks a value of glucose for one of these donors.
#'   \code{exposureInfo(x)} gives the name of each compound in the source
#'   table (\code{source_name}), its family, and \code{lowest_reported},
#'   the smallest value of the compound among the 387 samples of the table.
#'   The table reports no limits of detection. For the metabolites that
#'   they analysed, those with fewer than 20 per cent of the values
#'   missing, the authors put this smallest value in place of a missing
#'   one; six of the 13 compounds here have more missing values than that.
#'   The 13 compounds were chosen by their kind, from their names and the
#'   numbers of samples with a value.}
#'   \item{Plasma metabolites}{\code{metadata(x)$plasma_metabolites}, 43
#'   donors by 182 metabolites: the metabolites that the source table
#'   assigns to a super-pathway other than xenobiotics and reports for all
#'   43 donors, on the same relative scale. Their pathways are in
#'   \code{metadata(x)$metabolite_info}.}
#'   \item{Provenance}{\code{metadata(x)$source}: the reference, the
#'   sources and their licences, the permission for the single-cell values
#'   (\code{single_cells$permission}), the names, sizes and SHA-256 of the 88
#'   input files, the selection rule, the gates of each sample
#'   (\code{gates}) and the numbers of monocytes, T cells and other cells
#'   of each donor (\code{gated_cells}).}
#' }
#' The rule for the donors, the rule for the cell types, the numbers of
#' cells and the seed were fixed before any association of a donor variable
#' with expression, proteins or metabolites was computed on these donors.
#' The object was first made from 42 of the 43 donors, because the first
#' transfer of the expression file of one was incomplete, and the vignettes
#' were first written on that version. It was then made again from all 43
#' by the script without that exclusion; the cells of the 14 donors that
#' follow that donor in the order of the draw are a new random draw.
#'
#' @format An RDS file, \code{inst/extdata/incov_pbmc_subset.rds}, read with
#'   \code{readRDS(system.file("extdata", "incov_pbmc_subset.rds",
#'   package = "exposomeSC"))}.
#'
#' @source The clinical records and the plasma measurements are Table S1 of
#'   Mendeley Data, \doi{10.17632/tzydswhhb5.5} (version 5, published on 4
#'   November 2020), distributed under the Creative Commons Attribution 4.0
#'   International licence (CC BY 4.0,
#'   \url{https://creativecommons.org/licenses/by/4.0/}). The plasma values
#'   were selected, not changed; smoking and sex were recoded as described
#'   above.
#'
#'   The single-cell values of genes and surface proteins are from
#'   ArrayExpress accession E-MTAB-9357
#'   (\url{https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-9357}),
#'   released on 29 October 2020 and downloaded on 19 November 2026. On 21
#'   November 2026 the authors who deposited the accession gave permission
#'   in these words: "The single-cell values from E-MTAB-9357 included in
#'   exposomeSC may be used and redistributed under CC BY 4.0, with citation
#'   of Su et al. (2020)." The object holds a subset of the values that the
#'   authors released, transformed as described above.
#'
#'   The object as a whole is therefore distributed under CC BY 4.0. Work
#'   that uses it should cite Su et al. (2020). The permission does not imply
#'   that the authors endorse exposomeSC.
#'
#'   The article itself is published under CC BY-NC-ND 4.0. The values in
#'   this object come from the data that the authors deposited in Mendeley
#'   Data and ArrayExpress, not from the article.
#'
#'   The object was made by \code{inst/scripts/make_incov_pbmc_subset.R},
#'   which checks the size and the SHA-256 of every input file and carries
#'   out the steps described above.
#'
#' @references
#' Su Y, Chen D, Yuan D, et al. (2020). Multi-omics resolves a sharp
#' disease-state shift between mild and moderate COVID-19. \emph{Cell},
#' 183(6), 1479-1495.e20. \doi{10.1016/j.cell.2020.10.037}
#'
#' @examples
#' x <- readRDS(system.file("extdata", "incov_pbmc_subset.rds",
#'     package = "exposomeSC"))
#' x
#' table(x$cell_type)
#' exposureData(x)[1:3, 1:6]
#' SingleCellExperiment::altExp(x, "CITE")
#' dim(S4Vectors::metadata(x)$plasma_metabolites)
#' S4Vectors::metadata(x)$source$single_cells$licence
#'
#' @name incov_pbmc_subset
#' @docType data
#' @keywords datasets
NULL
