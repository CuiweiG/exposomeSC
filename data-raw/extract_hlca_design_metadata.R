## Extract the donor-, sample-, and cell-level design metadata required for
## confounding control and sampling-process sensitivity analyses. This script
## reads only /obs from the canonical HLCA H5AD; it never touches the assay.

.script_file <- function() {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
    if (length(file_arg)) {
        return(normalizePath(file_arg[[1]], mustWork = TRUE))
    }
    ofile <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
    if (!is.null(ofile)) {
        return(normalizePath(ofile, mustWork = TRUE))
    }
    stop("Cannot identify the script path; run with Rscript or source().")
}

project_root <- dirname(dirname(.script_file()))
h5ad_path <- file.path(project_root, "data", "hlca_core.h5ad")
output_path <- file.path(project_root, "data", "hlca_design_metadata.rds")
dataset_url <- paste0(
    "https://datasets.cellxgene.cziscience.com/",
    "688185ad-11c2-4172-a53a-f4f1f4076860.h5ad"
)
collection_url <- paste0(
    "https://cellxgene.cziscience.com/collections/",
    "6f6d381a-7701-4781-935c-db10d30de293"
)
expected_h5ad_bytes <- 5873612847

if (!requireNamespace("rhdf5", quietly = TRUE)) {
    stop("Package 'rhdf5' is required.")
}
if (!file.exists(h5ad_path)) {
    stop("Canonical HLCA H5AD not found: ", h5ad_path)
}
if (unname(file.info(h5ad_path)$size) != expected_h5ad_bytes) {
    stop("HLCA H5AD byte size differs from the pinned dataset version.")
}

h5_nodes <- rhdf5::h5ls(h5ad_path, recursive = TRUE)

.read_obs <- function(field) {
    node <- h5_nodes[
        h5_nodes$group == "/obs" & h5_nodes$name == field,
        ,
        drop = FALSE
    ]
    if (nrow(node) != 1L) {
        stop("Expected exactly one /obs node named '", field, "'.")
    }
    path <- paste0("/obs/", field)
    if (node$otype[[1]] == "H5I_GROUP") {
        categories <- rhdf5::h5read(
            h5ad_path,
            paste0(path, "/categories")
        )
        codes <- as.integer(rhdf5::h5read(
            h5ad_path,
            paste0(path, "/codes")
        ))
        value <- rep(NA_character_, length(codes))
        observed <- codes >= 0L
        if (any(codes[observed] >= length(categories))) {
            stop("Out-of-range categorical code in /obs/", field, ".")
        }
        value[observed] <- categories[codes[observed] + 1L]
        return(value)
    }
    rhdf5::h5read(h5ad_path, path)
}

fields <- c(
    "_index",
    "donor_id",
    "sample",
    "study",
    "dataset",
    "assay",
    "sequencing_platform",
    "fresh_or_frozen",
    "tissue_sampling_method",
    "tissue_level_2",
    "tissue_level_3",
    "smoking_status",
    "age_or_mean_of_age_range",
    "sex",
    "BMI",
    "ann_level_1",
    "ann_level_2",
    "ann_level_3",
    "ann_level_4",
    "ann_level_5",
    "ann_finest_level",
    "cell_type",
    "scanvi_label",
    "reannotation_type",
    "n_genes_detected",
    "log10_total_counts"
)

message("Reading HLCA /obs design fields from ", h5ad_path, " ...")
obs <- as.data.frame(
    setNames(lapply(fields, .read_obs), fields),
    stringsAsFactors = FALSE,
    check.names = FALSE
)
names(obs)[names(obs) == "_index"] <- "cell_id"

if (anyDuplicated(obs$cell_id)) {
    stop("HLCA /obs cell identifiers are not unique.")
}
if (anyNA(obs$donor_id) || anyNA(obs$sample)) {
    stop("Missing donor_id or sample in HLCA /obs.")
}

obs$smoking_ord <- unname(c(
    never = 0,
    former = 1,
    active = 2
)[obs$smoking_status])
obs$smoking_ever <- ifelse(
    is.na(obs$smoking_ord),
    NA_real_,
    as.numeric(obs$smoking_ord > 0)
)
obs$age <- suppressWarnings(as.numeric(obs$age_or_mean_of_age_range))
obs$sex_binary <- ifelse(
    obs$sex == "male",
    1,
    ifelse(obs$sex == "female", 0, NA_real_)
)
obs$BMI <- suppressWarnings(as.numeric(obs$BMI))

.assert_constant <- function(data, id, fields, level) {
    split_index <- split(seq_len(nrow(data)), data[[id]])
    violations <- character()
    for (field in fields) {
        bad <- vapply(split_index, function(index) {
            value <- unique(data[[field]][index])
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
        stop(
            level,
            " metadata are not constant: ",
            paste(violations, collapse = "; ")
        )
    }
    invisible(TRUE)
}

donor_constant <- c(
    "smoking_status",
    "smoking_ord",
    "smoking_ever",
    "age",
    "sex",
    "sex_binary",
    "BMI",
    "study",
    "dataset",
    "assay",
    "sequencing_platform",
    "fresh_or_frozen"
)
.assert_constant(obs, "donor_id", donor_constant, "Donor-level")

sample_constant <- c(
    "donor_id",
    donor_constant,
    "tissue_sampling_method",
    "tissue_level_2",
    "tissue_level_3"
)
.assert_constant(obs, "sample", sample_constant, "Sample-level")

.first_observed <- function(value) {
    value <- value[!is.na(value)]
    if (length(value)) value[[1]] else NA
}

.collapse_level <- function(data, id, fields) {
    ids <- unique(data[[id]])
    index <- split(seq_len(nrow(data)), data[[id]])[ids]
    out <- data.frame(ids, stringsAsFactors = FALSE)
    names(out) <- id
    for (field in fields) {
        template <- if (is.numeric(data[[field]])) NA_real_ else NA_character_
        out[[field]] <- vapply(
            index,
            function(i) {
                value <- data[[field]][i]
                value <- value[!is.na(value)]
                if (length(value)) value[[1]] else template
            },
            template
        )
    }
    out$n_cells <- lengths(index)
    out
}

drop_celltype <- c("None", "Rare", "Unknown", "NA", "nan")
keep <- !is.na(obs$smoking_ever) &
    !is.na(obs$ann_level_2) &
    !(obs$ann_level_2 %in% drop_celltype)
cell_metadata <- obs[keep, c(
    "cell_id",
    "donor_id",
    "sample",
    "study",
    "dataset",
    "assay",
    "sequencing_platform",
    "fresh_or_frozen",
    "tissue_sampling_method",
    "tissue_level_2",
    "tissue_level_3",
    "ann_level_1",
    "ann_level_2",
    "ann_level_3",
    "ann_level_4",
    "ann_level_5",
    "ann_finest_level",
    "cell_type",
    "scanvi_label",
    "reannotation_type",
    "n_genes_detected",
    "log10_total_counts"
)]
names(cell_metadata)[names(cell_metadata) == "cell_type"] <-
    "cell_type_original"
rownames(cell_metadata) <- cell_metadata$cell_id

donor_metadata <- .collapse_level(
    obs[keep, , drop = FALSE],
    "donor_id",
    donor_constant
)
sample_metadata <- .collapse_level(
    obs[keep, , drop = FALSE],
    "sample",
    sample_constant
)

samples_per_donor <- table(sample_metadata$donor_id)
donor_metadata$n_samples <- as.integer(samples_per_donor[donor_metadata$donor_id])

if (nrow(donor_metadata) != 98L) {
    stop("Expected 98 smoking-annotated donors; found ", nrow(donor_metadata), ".")
}
if (nrow(cell_metadata) != 483762L) {
    stop("Expected 483,762 retained cells; found ", nrow(cell_metadata), ".")
}
if (!identical(sort(unique(cell_metadata$donor_id)), sort(donor_metadata$donor_id))) {
    stop("Cell and donor metadata contain different donor sets.")
}

source_info <- list(
    title = "An integrated cell atlas of the human lung in health and disease (core)",
    publication_doi = "10.1038/s41591-023-02327-2",
    collection_url = collection_url,
    dataset_version_url = dataset_url,
    cellxgene_schema_version = "7.1.0",
    licence = "CC BY 4.0",
    path = normalizePath(h5ad_path, winslash = "/", mustWork = TRUE),
    bytes = unname(file.info(h5ad_path)$size),
    modified = format(file.info(h5ad_path)$mtime, tz = "UTC", usetz = TRUE),
    extracted = format(Sys.time(), tz = "UTC", usetz = TRUE),
    obs_rows = nrow(obs),
    retained_cells = nrow(cell_metadata),
    retained_donors = nrow(donor_metadata),
    retained_samples = nrow(sample_metadata),
    fields = fields
)

result <- list(
    donor = donor_metadata,
    sample = sample_metadata,
    cell = cell_metadata,
    source = source_info
)
saveRDS(result, output_path, compress = "gzip")

message(
    "Saved ", output_path, ": ",
    nrow(donor_metadata), " donors, ",
    nrow(sample_metadata), " samples, ",
    nrow(cell_metadata), " cells."
)
