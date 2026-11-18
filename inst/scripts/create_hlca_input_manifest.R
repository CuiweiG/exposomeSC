## Create a cryptographic manifest for the canonical HLCA source and derived
## analysis inputs. Hashing is intentionally explicit because the two large
## objects can take several minutes to stream from disk.
##
## Usage:
##   Rscript create_hlca_input_manifest.R <output_dir> <hlca_core.h5ad> \
##       <hlca_smoking_scee.rds> <hlca_design_metadata.rds> \
##       <hlca_feature_metadata.rds>
## The four inputs are the pinned H5AD and the files written by
## fetch_hlca.R, extract_hlca_design_metadata.R and
## extract_hlca_feature_metadata.R. hlca_input_manifest.rds and
## hlca_input_manifest.csv are written to <output_dir>; a manifest that is
## already there is used to skip the hashing of files that have not changed.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) {
    stop("Usage: Rscript create_hlca_input_manifest.R <output_dir> ",
         "<hlca_core.h5ad> <hlca_smoking_scee.rds> ",
         "<hlca_design_metadata.rds> <hlca_feature_metadata.rds>")
}
output_dir <- args[[1]]
paths <- unlist(args[2:5])

## The inputs are checked before the output directory is made
files <- basename(paths)
if (anyDuplicated(files)) {
    stop("The four inputs must have different file names.")
}
if (any(!file.exists(paths))) {
    stop(
        "Missing canonical input(s): ",
        paste(paths[!file.exists(paths)], collapse = ", ")
    )
}
if (!requireNamespace("digest", quietly = TRUE)) {
    stop("Package 'digest' is required.")
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(output_dir)) {
    stop("The output directory cannot be used: ", output_dir)
}

source_url <- c(
    paste0(
        "https://datasets.cellxgene.cziscience.com/",
        "688185ad-11c2-4172-a53a-f4f1f4076860.h5ad"
    ),
    rep(NA_character_, 3L)
)
role <- c(
    "immutable upstream source",
    "derived raw-count SCEE",
    "derived design metadata",
    "derived feature metadata"
)

previous_path <- file.path(output_dir, "hlca_input_manifest.rds")
previous <- if (file.exists(previous_path)) readRDS(previous_path) else NULL
current_bytes <- unname(file.info(paths)$size)
current_modified <- format(file.info(paths)$mtime, tz = "UTC", usetz = TRUE)

sha256 <- vapply(seq_along(paths), function(i) {
    previous_index <- if (!is.null(previous)) {
        match(files[[i]], previous$file)
    } else {
        NA_integer_
    }
    reusable <- !is.na(previous_index) &&
        identical(previous$bytes[[previous_index]], current_bytes[[i]]) &&
        identical(previous$modified_utc[[previous_index]], current_modified[[i]]) &&
        grepl("^[0-9a-f]{64}$", previous$sha256[[previous_index]])
    if (reusable) {
        message("Reusing verified hash for ", files[[i]], ".")
        return(previous$sha256[[previous_index]])
    }
    message(
        "Hashing ", files[[i]], " (",
        format(current_bytes[[i]], big.mark = ","), " bytes) ..."
    )
    digest::digest(
        paths[[i]],
        algo = "sha256",
        serialize = FALSE,
        file = TRUE
    )
}, character(1))

manifest <- data.frame(
    role = role,
    file = files,
    bytes = current_bytes,
    modified_utc = current_modified,
    sha256 = unname(sha256),
    source_url = source_url,
    publication_doi = c("10.1038/s41591-023-02327-2", rep(NA_character_, 3L)),
    licence = c("CC BY 4.0", rep("derived from CC BY 4.0 source", 3L)),
    generated_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    stringsAsFactors = FALSE
)

saveRDS(
    manifest,
    file.path(output_dir, "hlca_input_manifest.rds"),
    compress = "xz"
)
utils::write.csv(
    manifest,
    file.path(output_dir, "hlca_input_manifest.csv"),
    row.names = FALSE
)
message("Saved SHA-256 manifest for ", nrow(manifest), " canonical inputs.")
