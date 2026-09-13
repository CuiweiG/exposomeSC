## Create a cryptographic manifest for the canonical HLCA source and derived
## analysis inputs. Hashing is intentionally explicit because the two large
## objects can take several minutes to stream from disk.

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
summary_dir <- file.path(project_root, "analysis", "summaries")
result_dir <- file.path(project_root, "analysis", "results")
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

if (!requireNamespace("digest", quietly = TRUE)) {
    stop("Package 'digest' is required.")
}

relative_paths <- c(
    "data/hlca_core.h5ad",
    "data/hlca_smoking_scee.rds",
    "data/hlca_design_metadata.rds",
    "data/hlca_feature_metadata.rds"
)
paths <- file.path(project_root, relative_paths)
if (any(!file.exists(paths))) {
    stop(
        "Missing canonical input(s): ",
        paste(relative_paths[!file.exists(paths)], collapse = ", ")
    )
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

previous_path <- file.path(result_dir, "hlca_input_manifest.rds")
previous <- if (file.exists(previous_path)) readRDS(previous_path) else NULL
current_bytes <- unname(file.info(paths)$size)
current_modified <- format(file.info(paths)$mtime, tz = "UTC", usetz = TRUE)

sha256 <- vapply(seq_along(paths), function(i) {
    previous_index <- if (!is.null(previous)) {
        match(relative_paths[[i]], previous$relative_path)
    } else {
        NA_integer_
    }
    reusable <- !is.na(previous_index) &&
        identical(previous$bytes[[previous_index]], current_bytes[[i]]) &&
        identical(previous$modified_utc[[previous_index]], current_modified[[i]]) &&
        grepl("^[0-9a-f]{64}$", previous$sha256[[previous_index]])
    if (reusable) {
        message("Reusing verified hash for ", relative_paths[[i]], ".")
        return(previous$sha256[[previous_index]])
    }
    message(
        "Hashing ", relative_paths[[i]], " (",
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
    relative_path = gsub("\\\\", "/", relative_paths),
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
    file.path(result_dir, "hlca_input_manifest.rds"),
    compress = "xz"
)
utils::write.csv(
    manifest,
    file.path(summary_dir, "hlca_input_manifest.csv"),
    row.names = FALSE
)
message("Saved SHA-256 manifest for ", nrow(manifest), " canonical inputs.")
