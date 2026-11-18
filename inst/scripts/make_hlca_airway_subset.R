#!/usr/bin/env Rscript
## Build inst/extdata/hlca_airway_subset.rds, the real-data example of
## exposomeSC, from the Human Lung Cell Atlas core.
##
## Usage:
##   Rscript make_hlca_airway_subset.R <hlca_core.h5ad> <output.rds>
##
## Source data
##   Human Lung Cell Atlas (HLCA) core: Sikkema L et al. (2023) An integrated
##   cell atlas of the lung in health and disease. Nature Medicine 29,
##   1563-1577, doi:10.1038/s41591-023-02327-2. Licence: CC BY 4.0.
##   File: the CELLxGENE H5AD of dataset version
##   688185ad-11c2-4172-a53a-f4f1f4076860,
##   https://datasets.cellxgene.cziscience.com/688185ad-11c2-4172-a53a-f4f1f4076860.h5ad
##   (5,873,612,847 bytes; SHA-256
##   1cbdde1e513a31bdb6a3f70296e2d0612760d2af4032e239c94596ea87b44007).
##   The script checks the size and the SHA-256 before it reads the file.
##   Raw counts are read from /raw/X.
##
## Selection (the donors and the cell types were fixed before any result on
## cells of these donors was seen; the number of cells per donor and cell
## type was set later, by the size limit of a source package: before this
## subset was drawn, but after the smoking model had been run on larger sets
## of cells of the same donors)
##   - HLCA study "Seibold_2020", which is Goldfarbmuren KC et al. (2020)
##     Dissecting the cellular specificity of smoking effects and
##     reconstructing lineages in the human airway epithelium. Nature
##     Communications 11, 2485, doi:10.1038/s41467-020-16239-z (GEO
##     GSE134174): tracheal epithelium of smokers and never smokers.
##   - 12 donors: the six heavy smokers T85, T90, T101, T120, T154 and T167
##     (smokers with at least 15 pack-years, as the paper defines them) and
##     the six never smokers T126, T137, T153, T164, T165 and T166, the two
##     groups that the paper compares, as classified in the authors' public
##     metadata (github.com/seiboldlab/SingleCell_smoking, Metadata.txt,
##     columns smoke and smoke_noT89). The other two donors of the study in
##     the atlas are left out: T121, a light smoker, and T89, a child who
##     never smoked and whom the paper leaves out of its comparison.
##     The HLCA writes a donor as a long identifier that ends in the study's
##     code (for example "..._483747T101"); the script takes the code from
##     the end and checks that codes and identifiers correspond one to one.
##     It also checks that the HLCA records the heavy smokers as "active"
##     and the never smokers as "never" smokers.
##   - Cells annotated "Basal", "Secretory" or "Multiciliated lineage" at
##     HLCA annotation level 3 (ann_level_3).
##   - Up to 10 cells per donor and cell type, drawn at random without
##     looking at expression: the strata are taken in the order of the donor
##     codes and cell-type names (sorted as in the C locale), the cells of a
##     stratum in the order of the file, and sample.int() draws from them
##     after
##     set.seed(20261030, kind = "Mersenne-Twister",
##              normal.kind = "Inversion", sample.kind = "Rejection").
##     The number of eligible cells of each donor and cell type in the
##     atlas, 29,217 in all, is kept in metadata(x)$source$eligible_cells.
##   - Genes without a count in the selected cells are left out; the
##     column sums, and so the library sizes of every pseudobulk, are those
##     of the complete transcriptome.
##
## Donor variables (exposureData of the result)
##   smoking       1 = heavy smoker, 0 = never smoker
##   chemistry_v3  1 = 10x Genomics 3' v3 chemistry, 0 = 3' v2 (from 'assay')
##   age           age in years ('age_or_mean_of_age_range')
##   sex_male      1 = male, 0 = female ('sex')
##
## Packages: rhdf5, Matrix, S4Vectors, SingleCellExperiment, exposomeSC.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
    stop("Usage: Rscript make_hlca_airway_subset.R <hlca_core.h5ad> ",
         "<output.rds>")
}
h5ad <- args[[1]]
out_file <- args[[2]]
## The paths are checked before any package is loaded or any data are read
if (!file.exists(h5ad)) stop("File not found: ", h5ad)
if (!dir.exists(dirname(out_file))) {
    stop("The directory of the output file does not exist: ",
         dirname(out_file))
}
expected_bytes <- 5873612847
expected_sha256 <-
    "1cbdde1e513a31bdb6a3f70296e2d0612760d2af4032e239c94596ea87b44007"
if (file.size(h5ad) != expected_bytes) {
    stop("The file has ", file.size(h5ad), " bytes, not the ",
         format(expected_bytes, big.mark = ","), " of the pinned HLCA H5AD.")
}
for (pkg in c("rhdf5", "Matrix", "S4Vectors", "SingleCellExperiment",
              "exposomeSC")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        stop("Package '", pkg, "' is required.")
    }
}
if (unname(tools::sha256sum(h5ad)) != expected_sha256) {
    stop("The SHA-256 of the file is not that of the pinned HLCA H5AD.")
}

## ---- cell annotations ---------------------------------------------------
## A categorical column of /obs is a group with integer codes (0-based, -1 for
## a missing value) and its categories.
obs_nodes <- rhdf5::h5ls(h5ad, recursive = 2, datasetinfo = FALSE)
obs_nodes <- obs_nodes[obs_nodes$group == "/obs", ]
read_obs <- function(name) {
    path <- paste0("obs/", name)
    if (obs_nodes$otype[match(name, obs_nodes$name)] == "H5I_GROUP") {
        categories <- as.character(rhdf5::h5read(h5ad,
            paste0(path, "/categories")))
        codes <- as.integer(rhdf5::h5read(h5ad, paste0(path, "/codes")))
        out <- rep(NA_character_, length(codes))
        out[codes >= 0L] <- categories[codes[codes >= 0L] + 1L]
        return(out)
    }
    as.vector(rhdf5::h5read(h5ad, path))
}
obs <- data.frame(
    cell_id = as.character(rhdf5::h5read(h5ad, "obs/_index")),
    study = read_obs("study"),
    hlca_donor_id = read_obs("donor_id"),
    cell_type = read_obs("ann_level_3"),
    assay = read_obs("assay"),
    smoking_status = read_obs("smoking_status"),
    age = read_obs("age_or_mean_of_age_range"),
    sex = read_obs("sex"),
    stringsAsFactors = FALSE
)
obs$row <- seq_len(nrow(obs))

heavy <- c("T85", "T90", "T101", "T120", "T154", "T167")
never <- c("T126", "T137", "T153", "T164", "T165", "T166")
donors <- sort(c(heavy, never), method = "radix")
celltypes <- sort(c("Basal", "Secretory", "Multiciliated lineage"),
                  method = "radix")

## The donor code of the study is the end of the HLCA donor identifier.
study_cells <- obs[!is.na(obs$study) & obs$study == "Seibold_2020", ]
if (!nrow(study_cells) ||
    !all(grepl("T[0-9]+$", study_cells$hlca_donor_id))) {
    stop("The donor identifiers of study Seibold_2020 do not all end in a ",
         "donor code of the form T<number>.")
}
study_cells$donor_id <- regmatches(study_cells$hlca_donor_id,
    regexpr("T[0-9]+$", study_cells$hlca_donor_id))
pairs <- unique(study_cells[, c("donor_id", "hlca_donor_id")])
if (anyDuplicated(pairs$donor_id) || anyDuplicated(pairs$hlca_donor_id)) {
    stop("Donor codes and HLCA donor identifiers do not correspond one to ",
         "one in study Seibold_2020.")
}

sel <- study_cells[study_cells$donor_id %in% donors &
                   study_cells$cell_type %in% celltypes, ]
status <- unique(sel[, c("donor_id", "smoking_status")])
if (anyDuplicated(status$donor_id) ||
    !setequal(status$donor_id, donors) ||
    !all(status$smoking_status[match(heavy, status$donor_id)] == "active") ||
    !all(status$smoking_status[match(never, status$donor_id)] == "never")) {
    stop("The smoking status of the 12 donors in the file is not the ",
         "expected one.")
}
## The number of eligible cells of each donor and cell type in the atlas,
## kept in the provenance record: the cells are drawn from these.
eligible_cells <- table(factor(sel$donor_id, levels = donors),
                        factor(sel$cell_type, levels = celltypes))
eligible_cells <- matrix(as.integer(eligible_cells), nrow = length(donors),
                         dimnames = list(donors, celltypes))

## ---- random selection of up to 10 cells per donor and cell type ---------
set.seed(20261030, kind = "Mersenne-Twister", normal.kind = "Inversion",
         sample.kind = "Rejection")
picked <- integer()
for (donor in donors) {
    for (ct in celltypes) {
        rows <- sort(sel$row[sel$donor_id == donor & sel$cell_type == ct])
        if (!length(rows)) next
        take <- sort(rows[sample.int(length(rows), min(10L, length(rows)))])
        picked <- c(picked, take)
    }
}
cells <- sel[match(picked, sel$row), ]

## ---- raw counts of the selected cells -----------------------------------
## /raw/X is stored by row (CSR): the counts of cell r (1-based) are
## data[indptr[r] + 1 .. indptr[r + 1]] at the gene positions indices[...].
indptr <- as.numeric(rhdf5::h5read(h5ad, "raw/X/indptr",
    bit64conversion = "double"))
gene_ids <- as.character(rhdf5::h5read(h5ad, "raw/var/_index"))
symbol_codes <- as.integer(rhdf5::h5read(h5ad, "raw/var/feature_name/codes"))
symbol_cats <- as.character(rhdf5::h5read(h5ad,
    "raw/var/feature_name/categories"))
symbols <- symbol_cats[symbol_codes + 1L]
ii <- vector("list", nrow(cells))
xx <- vector("list", nrow(cells))
for (j in seq_len(nrow(cells))) {
    r <- cells$row[[j]]
    first <- indptr[[r]] + 1
    n <- indptr[[r + 1L]] - indptr[[r]]
    if (n == 0) next
    ii[[j]] <- as.integer(rhdf5::h5read(h5ad, "raw/X/indices",
        start = first, count = n)) + 1L
    xx[[j]] <- as.numeric(rhdf5::h5read(h5ad, "raw/X/data",
        start = first, count = n))
}
jj <- rep(seq_along(ii), lengths(ii))
counts <- Matrix::sparseMatrix(
    i = unlist(ii), j = jj, x = unlist(xx),
    dims = c(length(gene_ids), nrow(cells)),
    dimnames = list(gene_ids, cells$cell_id)
)
if (any(counts@x != round(counts@x)) || any(counts@x < 0)) {
    stop("/raw/X does not hold non-negative whole counts.")
}
keep <- Matrix::rowSums(counts) > 0
counts <- counts[keep, ]

## ---- the SingleCellExposomeExperiment -----------------------------------
sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = counts),
    colData = S4Vectors::DataFrame(
        donor_id = cells$donor_id,
        cell_type = cells$cell_type,
        assay = cells$assay,
        hlca_donor_id = cells$hlca_donor_id,
        row.names = cells$cell_id
    ),
    rowData = S4Vectors::DataFrame(
        symbol = symbols[keep],
        row.names = gene_ids[keep]
    )
)
donor_table <- unique(cells[, c("donor_id", "smoking_status", "assay",
                                "age", "sex")])
if (anyDuplicated(donor_table$donor_id)) {
    stop("A donor has more than one value of a donor-level variable.")
}
donor_table <- donor_table[match(donors, donor_table$donor_id), ]
exposures <- cbind(
    smoking = as.numeric(donor_table$smoking_status == "active"),
    chemistry_v3 = as.numeric(donor_table$assay == "10x 3' v3"),
    age = as.numeric(donor_table$age),
    sex_male = as.numeric(donor_table$sex == "male")
)
rownames(exposures) <- donor_table$donor_id
info <- S4Vectors::DataFrame(
    exposure = colnames(exposures),
    family = c("tobacco smoke", "technical", "demographic", "demographic"),
    unit = c("1 = heavy smoker, 0 = never smoker",
             "1 = 10x 3' v3, 0 = 10x 3' v2", "years", "1 = male, 0 = female")
)
scee <- exposomeSC::build_scee(sce, exposures, sample_col = "donor_id",
                               exposure_info = info)
S4Vectors::metadata(scee)$source <- list(
    data = "Human Lung Cell Atlas (HLCA) core, study Seibold_2020",
    url = paste0("https://datasets.cellxgene.cziscience.com/",
                 "688185ad-11c2-4172-a53a-f4f1f4076860.h5ad"),
    dataset_version = "688185ad-11c2-4172-a53a-f4f1f4076860",
    bytes = expected_bytes,
    sha256 = expected_sha256,
    licence = "CC BY 4.0",
    references = c("doi:10.1038/s41591-023-02327-2",
                   "doi:10.1038/s41467-020-16239-z"),
    selection = paste("12 donors (6 heavy, 6 never smokers);",
                      "Basal, Secretory and Multiciliated lineage cells",
                      "(ann_level_3); up to 10 cells per donor and cell",
                      "type, seed 20261030"),
    eligible_cells = eligible_cells,
    script = "inst/scripts/make_hlca_airway_subset.R"
)
saveRDS(scee, out_file, compress = "xz")
message(sprintf("%d cells, %d genes, %d donors; %s bytes written to %s",
                ncol(scee), nrow(scee), nrow(exposures),
                format(file.size(out_file), big.mark = ","), out_file))
