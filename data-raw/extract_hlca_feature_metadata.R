## Extract immutable feature annotations from the canonical HLCA H5AD.
## This script reads only /var and writes a small, derived RDS used to annotate
## Ensembl-keyed statistical results. It does not read or transform counts.

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
output_path <- file.path(project_root, "data", "hlca_feature_metadata.rds")
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
if (!requireNamespace("AnnotationDbi", quietly = TRUE) ||
        !requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    stop("Packages 'AnnotationDbi' and 'org.Hs.eg.db' are required.")
}
if (!file.exists(h5ad_path)) {
    stop("Canonical HLCA H5AD not found: ", h5ad_path)
}
if (unname(file.info(h5ad_path)$size) != expected_h5ad_bytes) {
    stop("HLCA H5AD byte size differs from the pinned dataset version.")
}

h5_nodes <- rhdf5::h5ls(h5ad_path, recursive = TRUE)

.read_var <- function(field) {
    node <- h5_nodes[
        h5_nodes$group == "/var" & h5_nodes$name == field,
        ,
        drop = FALSE
    ]
    if (nrow(node) != 1L) {
        stop("Expected exactly one /var node named '", field, "'.")
    }
    path <- paste0("/var/", field)
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
            stop("Out-of-range categorical code in /var/", field, ".")
        }
        value[observed] <- categories[codes[observed] + 1L]
        return(value)
    }
    rhdf5::h5read(h5ad_path, path)
}

fields <- c(
    "_index",
    "feature_name",
    "feature_biotype",
    "feature_type",
    "feature_length",
    "feature_reference",
    "feature_is_filtered"
)

message("Reading HLCA /var annotations from ", h5ad_path, " ...")
feature <- as.data.frame(
    setNames(lapply(fields, .read_var), fields),
    stringsAsFactors = FALSE,
    check.names = FALSE
)
names(feature)[names(feature) == "_index"] <- "ensembl_id"
names(feature)[names(feature) == "feature_name"] <- "gene_symbol"
feature$feature_length <- suppressWarnings(as.integer(feature$feature_length))
feature$feature_is_filtered <- as.logical(feature$feature_is_filtered)

if (nrow(feature) != 27402L) {
    stop("Expected 27,402 HLCA features; found ", nrow(feature), ".")
}
if (anyNA(feature$ensembl_id) || anyDuplicated(feature$ensembl_id)) {
    stop("HLCA Ensembl feature identifiers must be complete and unique.")
}
if (anyNA(feature$gene_symbol) || any(!nzchar(feature$gene_symbol))) {
    stop("HLCA feature_name contains missing or empty gene symbols.")
}
rownames(feature) <- feature$ensembl_id

.map_orgdb <- function(column) {
    unname(AnnotationDbi::mapIds(
        org.Hs.eg.db::org.Hs.eg.db,
        keys = feature$ensembl_id,
        column = column,
        keytype = "ENSEMBL",
        multiVals = "first"
    ))
}

message("Adding auxiliary org.Hs.eg.db annotations ...")
feature$orgdb_symbol <- .map_orgdb("SYMBOL")
feature$gene_name <- .map_orgdb("GENENAME")
feature$entrez_id <- .map_orgdb("ENTREZID")
feature$chromosome <- .map_orgdb("CHR")
feature$sex_chromosome <- feature$chromosome %in% c("X", "Y")
feature$annotation_mapping_status <- ifelse(
    is.na(feature$orgdb_symbol),
    "unmapped",
    ifelse(
        feature$orgdb_symbol == feature$gene_symbol,
        "symbol-concordant",
        "symbol-discordant"
    )
)

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
    var_rows = nrow(feature),
    fields = fields,
    auxiliary_annotation = list(
        package = "org.Hs.eg.db",
        version = as.character(utils::packageVersion("org.Hs.eg.db")),
        keytype = "ENSEMBL",
        multiVals = "first"
    )
)

saveRDS(
    list(feature = feature, source = source_info),
    output_path,
    compress = "gzip"
)
message("Saved ", output_path, ": ", nrow(feature), " annotated features.")
