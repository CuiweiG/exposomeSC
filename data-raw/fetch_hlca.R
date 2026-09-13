## data-raw/fetch_hlca.R
## Build the HLCA-core smoking SingleCellExposomeExperiment from the CZ
## CELLxGENE download, PINNED TO THE VERIFIED OBS SCHEMA of the actual file
## (inspected 2026-07-03 from data/hlca_core.h5ad).
##
## Verified file facts:
##   584,944 cells x 27,402 genes; sparse CSR. X = log-normalised;
##   RAW COUNTS live in /raw/X (same 27,402 genes) -> WE USE RAW COUNTS.
##   107 donors (obs$donor_id); all disease == "normal".
##   obs$smoking_status = factor c("active","former","never")
##     (cells: never 301791, former 92329, active 89642, NA 101182).
##   obs$BMI (float64, some NaN), obs$age_or_mean_of_age_range (float64),
##   obs$sex = c("female","male").
##   Cell-type annotations: ann_level_1..5, ann_finest_level, cell_type.
##
## Exposure coding (numeric, because exposureData must be a numeric matrix):
##   smoking_ever : never=0, {former,active}=1   (PRIMARY exposure)
##   smoking_ord  : never=0, former=1, active=2   (ordinal alternative)
## Covariates: age (years), sex (female=0/male=1), BMI (kept, may be NA).
##
## STATUS: executed against the pinned local asset. Re-run after changes to the
## S4 container or provenance fields, then regenerate the input manifest.
## Windows: no fork; large matrix read on-disk via use_hdf5, realised to a
##   sparse dgCMatrix only for the retained (non-missing-smoking) cells.

## ------------------------- configuration --------------------------------
.script_file <- function() {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
    if (length(file_arg))
        return(normalizePath(file_arg[[1]], mustWork = TRUE))
    ofile <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
    if (!is.null(ofile)) return(normalizePath(ofile, mustWork = TRUE))
    stop("Cannot identify the script path; run with Rscript or source().")
}
PROJECT_ROOT <- dirname(dirname(.script_file()))
IN_H5AD   <- file.path(PROJECT_ROOT, "data", "hlca_core.h5ad")
OUT_SCEE  <- file.path(PROJECT_ROOT, "data", "hlca_smoking_scee.rds")
DATASET_URL <- paste0(
    "https://datasets.cellxgene.cziscience.com/",
    "688185ad-11c2-4172-a53a-f4f1f4076860.h5ad"
)
COLLECTION_URL <- paste0(
    "https://cellxgene.cziscience.com/collections/",
    "6f6d381a-7701-4781-935c-db10d30de293"
)
PUBLICATION_DOI <- "10.1038/s41591-023-02327-2"
EXPECTED_H5AD_BYTES <- 5873612847
CELLTYPE_FIELD <- "ann_level_2"     # 11 clean types; use "ann_level_3" for finer
## The HLCA matrix is SPARSE (CSR). reader="R" loads it as an in-memory sparse
## dgCMatrix (a few GB of the ~488M nonzeros), which is correct here. Do NOT set
## use_hdf5=TRUE: zellkonverter writes DENSE HDF5, i.e. ~64 GB for 585k x 27k.
USE_HDF5  <- FALSE                   # keep FALSE (sparse in-memory)
DROP_CT   <- c("None", "Rare", "Unknown", "NA", "nan")  # non-cell-type labels

need <- c("devtools", "zellkonverter", "SingleCellExperiment",
          "SummarizedExperiment", "S4Vectors", "Matrix")
miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) stop("Install first: ",
                       paste(miss, collapse = ", "))
suppressPackageStartupMessages({
    library(SingleCellExperiment); library(SummarizedExperiment); library(S4Vectors)
})
devtools::load_all(PROJECT_ROOT, quiet = TRUE)
stopifnot(file.exists(IN_H5AD))
if (unname(file.info(IN_H5AD)$size) != EXPECTED_H5AD_BYTES) {
    stop("HLCA H5AD byte size differs from the pinned dataset version.")
}
manifest_path <- file.path(
    PROJECT_ROOT,
    "analysis",
    "results",
    "hlca_input_manifest.rds"
)
input_sha256 <- NA_character_
if (file.exists(manifest_path)) {
    manifest <- readRDS(manifest_path)
    manifest_index <- match("data/hlca_core.h5ad", manifest$relative_path)
    if (!is.na(manifest_index) &&
            identical(manifest$bytes[[manifest_index]], EXPECTED_H5AD_BYTES)) {
        input_sha256 <- manifest$sha256[[manifest_index]]
    }
}

## ------------------------- read the .h5ad -------------------------------
## reader = "R" uses rhdf5 (installed) -- no Python/basilisk needed.
message("Reading obs/var metadata from ", IN_H5AD, " (skip_assays) ...")
sce <- tryCatch(
    zellkonverter::readH5AD(IN_H5AD, reader = "R", raw = FALSE,
                            skip_assays = TRUE, verbose = TRUE),
    error = function(e) {
        message("  skip_assays read failed (", conditionMessage(e),
                "); reading X too ...")
        zellkonverter::readH5AD(IN_H5AD, reader = "R", raw = FALSE,
                                verbose = TRUE)
    })
message("Metadata: ", nrow(sce), " genes x ", ncol(sce), " cells")

## ------------------------- RAW COUNTS from /raw/X -----------------------
## CELLxGENE stores RAW counts in /raw/X (CSR sparse); /X is log-normalised.
## zellkonverter's raw import is unreliable across readers (it gave only X),
## so read /raw/X directly with rhdf5 and build a genes x cells dgCMatrix.
.read_raw_counts_h5ad <- function(path) {
    if (!requireNamespace("rhdf5", quietly = TRUE)) stop("rhdf5 required")
    grp <- "raw/X"
    tab <- rhdf5::h5ls(path)
    if (!any(tab$group == "/raw" & tab$name == "X")) {
        stop(
            "/raw/X not found. The HLCA analysis requires raw integer counts; ",
            "normalised /X is not a valid fallback for negative-binomial models."
        )
    }
    at  <- rhdf5::h5readAttributes(path, grp)
    shp <- as.integer(at[["shape"]])                 # c(n_obs, n_var)
    if (length(shp) != 2L) stop("bad 'shape' attribute on ", grp)
    message("  reading /", grp, " (", shp[1], " cells x ", shp[2], " genes) ...")
    data    <- as.numeric(rhdf5::h5read(path, paste0(grp, "/data")))
    indices <- as.integer(rhdf5::h5read(path, paste0(grp, "/indices")))
    indptr  <- rhdf5::h5read(path, paste0(grp, "/indptr"))
    nnz <- length(data)
    message("  nnz = ", nnz)
    if (nnz > .Machine$integer.max)
        stop("nnz exceeds 2^31 -- a chunked reader is needed.")
    ## CSR (cells x genes) -> dgRMatrix -> transpose -> dgCMatrix (genes x cells)
    m <- new("dgRMatrix", j = indices, p = as.integer(indptr),
             x = data, Dim = shp)
    methods::as(Matrix::t(m), "CsparseMatrix")
}
counts_all <- .read_raw_counts_h5ad(IN_H5AD)
stopifnot(ncol(counts_all) == ncol(sce), nrow(counts_all) == nrow(sce))
rownames(counts_all) <- rownames(sce)
colnames(counts_all) <- colnames(sce)

## ------------------------- obs -> donor-level table ---------------------
cd <- as.data.frame(colData(sce))
get1 <- function(col) {
    if (!col %in% colnames(cd)) stop("obs column '", col, "' not found")
    as.character(cd[[col]])
}
donor <- get1("donor_id")
if (!CELLTYPE_FIELD %in% colnames(cd))
    stop("CELLTYPE_FIELD '", CELLTYPE_FIELD, "' not in obs")
ctype <- as.character(cd[[CELLTYPE_FIELD]])

smk   <- get1("smoking_status")                       # active/former/never/NA
bmi   <- suppressWarnings(as.numeric(cd[["BMI"]]))
age   <- suppressWarnings(as.numeric(cd[["age_or_mean_of_age_range"]]))
bmi[is.nan(bmi)] <- NA_real_
age[is.nan(age)] <- NA_real_
sexch <- get1("sex")
sexn  <- ifelse(sexch == "male", 1, ifelse(sexch == "female", 0, NA_real_))
sample_id <- get1("sample")
study <- get1("study")
dataset <- get1("dataset")
assay_type <- get1("assay")
platform <- get1("sequencing_platform")
fresh_frozen <- get1("fresh_or_frozen")
sampling_method <- get1("tissue_sampling_method")
tissue_level_2 <- get1("tissue_level_2")
tissue_level_3 <- get1("tissue_level_3")

## numeric smoking codes (per cell; donor-constant)
smk_ord  <- c(never = 0, former = 1, active = 2)[smk]
smk_ever <- ifelse(is.na(smk_ord), NA_real_, as.numeric(smk_ord > 0))

## one row per donor (validate donor-constancy of exposure/covariates)
cell_df <- data.frame(donor = donor, smoking_ever = smk_ever,
                      smoking_ord = as.numeric(smk_ord), age = age,
                      sex = sexn, BMI = bmi, study = study,
                      dataset = dataset, assay = assay_type,
                      sequencing_platform = platform,
                      fresh_or_frozen = fresh_frozen,
                      stringsAsFactors = FALSE)

.assert_constant <- function(data, id, fields, level) {
    index <- split(seq_len(nrow(data)), data[[id]])
    violations <- character()
    for (field in fields) {
        bad <- vapply(index, function(i) {
            value <- unique(data[[field]][i])
            value <- value[!is.na(value)]
            length(value) > 1L
        }, logical(1))
        if (any(bad)) {
            violations <- c(
                violations,
                paste0(field, " in ", paste(names(bad)[bad], collapse = ", "))
            )
        }
    }
    if (length(violations)) {
        stop(level, " metadata are not constant: ",
             paste(violations, collapse = "; "))
    }
}

.assert_constant(
    cell_df,
    "donor",
    setdiff(colnames(cell_df), "donor"),
    "Donor-level"
)
sample_df <- data.frame(
    sample = sample_id,
    donor = donor,
    study = study,
    dataset = dataset,
    assay = assay_type,
    sequencing_platform = platform,
    fresh_or_frozen = fresh_frozen,
    tissue_sampling_method = sampling_method,
    tissue_level_2 = tissue_level_2,
    tissue_level_3 = tissue_level_3,
    stringsAsFactors = FALSE
)
.assert_constant(
    sample_df,
    "sample",
    setdiff(colnames(sample_df), "sample"),
    "Sample-level"
)

donor_tab <- cell_df[!duplicated(cell_df$donor), , drop = FALSE]
rownames(donor_tab) <- donor_tab$donor
## keep only donors with a defined smoking status
donor_tab <- donor_tab[!is.na(donor_tab$smoking_ever), , drop = FALSE]
message("Donors total: ", length(unique(donor)),
        " | with smoking status: ", nrow(donor_tab))
message("Donor smoking (ever) table: ",
        paste(names(table(donor_tab$smoking_ever)),
              table(donor_tab$smoking_ever), sep = "=", collapse = "  "))

## ------------------------- subset cells ---------------------------------
keep <- donor %in% rownames(donor_tab) &
        !(ctype %in% DROP_CT) & !is.na(ctype)
message("Retaining ", sum(keep), " / ", length(keep), " cells")
counts_sub <- counts_all[, keep, drop = FALSE]
## realise the retained block to an in-memory sparse dgCMatrix (portable rds)
counts_sub <- as(counts_sub, "CsparseMatrix")
donor_sub  <- donor[keep]; ct_sub <- ctype[keep]

## ------------------------- assemble SCE + SCEE --------------------------
sub_sce <- SingleCellExperiment(
    assays = list(counts = counts_sub),
    rowData = rowData(sce),
    colData = DataFrame(cell_id = colnames(counts_sub),
                        donor_id = donor_sub,
                        sample = sample_id[keep],
                        sample_id = sample_id[keep],
                        study = study[keep],
                        dataset = dataset[keep],
                        assay = assay_type[keep],
                        sequencing_platform = platform[keep],
                        fresh_or_frozen = fresh_frozen[keep],
                        tissue_sampling_method = sampling_method[keep],
                        tissue_level_2 = tissue_level_2[keep],
                        tissue_level_3 = tissue_level_3[keep],
                        cell_type = ct_sub))
rownames(sub_sce) <- rownames(sce)

exp_mat <- as.matrix(donor_tab[sort(unique(donor_sub)),
                               c("smoking_ever", "smoking_ord",
                                 "age", "sex", "BMI"), drop = FALSE])
exp_info <- DataFrame(
    exposure = colnames(exp_mat),
    family = c("lifestyle", "lifestyle", "demographic", "demographic",
               "anthropometric"),
    unit = c("binary", "ordinal", "years", "binary", "kg/m2"))

scee <- exposomeSC::build_scee(sub_sce, exp_mat, sample_col = "donor_id",
                               exposure_info = exp_info)
S4Vectors::metadata(scee)$source <- list(
    asset = "HLCA core H5AD from CZ CELLxGENE",
    title = "An integrated cell atlas of the human lung in health and disease (core)",
    publication_doi = PUBLICATION_DOI,
    collection_url = COLLECTION_URL,
    dataset_version_url = DATASET_URL,
    cellxgene_schema_version = "7.1.0",
    licence = "CC BY 4.0",
    input_path = normalizePath(IN_H5AD, winslash = "/", mustWork = TRUE),
    input_bytes = unname(file.info(IN_H5AD)$size),
    input_sha256 = input_sha256,
    raw_counts_path = "/raw/X",
    built = format(Sys.time(), tz = "UTC", usetz = TRUE),
    celltype_field = CELLTYPE_FIELD,
    retained_cells = ncol(scee),
    retained_donors = length(unique(scee$donor_id))
)
if (nrow(scee) != 27402L || ncol(scee) != 483762L ||
        length(unique(scee$donor_id)) != 98L) {
    stop(
        "Pinned HLCA retained dimensions changed: ",
        nrow(scee), " features x ", ncol(scee), " cells from ",
        length(unique(scee$donor_id)), " donors."
    )
}
if (anyDuplicated(rownames(scee)) || anyDuplicated(colnames(scee))) {
    stop("SCEE feature and cell identifiers must be unique.")
}
count_values <- SummarizedExperiment::assay(scee, "counts")@x
validation_chunk <- 5000000L
for (chunk_start in seq.int(1L, length(count_values), by = validation_chunk)) {
    chunk_end <- min(length(count_values), chunk_start + validation_chunk - 1L)
    value <- count_values[chunk_start:chunk_end]
    if (any(!is.finite(value)) || any(value < 0) ||
            any(value != floor(value))) {
        stop(
            "The retained /raw/X matrix is not a non-negative integer ",
            "count matrix."
        )
    }
}
methods::validObject(scee)

## Write to the same directory and swap only after saveRDS completes, so an
## interrupted multi-gigabyte serialisation cannot truncate the canonical file.
temporary_output <- paste0(OUT_SCEE, ".tmp-", Sys.getpid())
backup_output <- paste0(OUT_SCEE, ".previous")
if (file.exists(temporary_output) || file.exists(backup_output)) {
    stop("Refusing to overwrite a stale temporary or backup SCEE file.")
}
on.exit(unlink(temporary_output), add = TRUE)
message("Serialising validated SCEE to temporary file ...")
saveRDS(scee, temporary_output, compress = "gzip")
if (!file.exists(temporary_output) || file.info(temporary_output)$size <= 0) {
    stop("Temporary SCEE serialisation is missing or empty.")
}
had_previous <- file.exists(OUT_SCEE)
if (had_previous && !file.rename(OUT_SCEE, backup_output)) {
    stop("Could not stage the previous canonical SCEE for replacement.")
}
if (!file.rename(temporary_output, OUT_SCEE)) {
    if (had_previous) file.rename(backup_output, OUT_SCEE)
    stop("Could not atomically promote the rebuilt SCEE.")
}
if (had_previous && unlink(backup_output) != 0L && file.exists(backup_output)) {
    warning("Rebuilt SCEE is valid, but the previous-file backup remains: ",
            backup_output, call. = FALSE)
}

## ------------------------- QC report ------------------------------------
message("\n================= HLCA smoking SCEE built =================")
message("cells: ", ncol(scee), " | genes: ", nrow(scee),
        " | donors: ", length(unique(scee$donor_id)))
message("cell types (", CELLTYPE_FIELD, ", n=",
        length(unique(scee$cell_type)), "): ",
        paste(sort(unique(scee$cell_type)), collapse = ", "))
message("BMI missing in ", sum(is.na(exp_mat[, "BMI"])), " / ",
        nrow(exp_mat), " donors")
message("saved: ", OUT_SCEE)
message("NEXT: data-raw/extract_hlca_design_metadata.R, then ",
        "data-raw/extract_hlca_feature_metadata.R")
