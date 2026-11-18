.hlca_subset <- function() {
    readRDS(system.file("extdata", "hlca_airway_subset.rds",
        package = "exposomeSC"))
}

test_that("the HLCA subset is a valid object with the documented content", {
    x <- .hlca_subset()
    expect_s4_class(x, "SingleCellExposomeExperiment")
    expect_true(methods::validObject(x))
    expect_identical(dim(x), c(20111L, 352L))
    expect_identical(SummarizedExperiment::assayNames(x), "counts")
    counts <- SummarizedExperiment::assay(x, "counts")
    expect_s4_class(counts, "dgCMatrix")
    ## raw counts: whole numbers, and no gene without a count
    expect_true(all(counts@x == round(counts@x)))
    expect_true(all(counts@x > 0))
    expect_identical(sum(counts@x), 11337267)
    expect_true(all(Matrix::rowSums(counts) > 0))
    expect_true(all(startsWith(rownames(x), "ENSG")))
    ## the gene names of the atlas: for 2,301 genes the identifier again,
    ## and three names on two rows each
    symbol <- SummarizedExperiment::rowData(x)$symbol
    expect_false(anyNA(symbol))
    expect_identical(symbol[match("ENSG00000138061", rownames(x))], "CYP1B1")
    expect_identical(sum(symbol == rownames(x)), 2301L)
    expect_identical(sort(symbol[duplicated(symbol)]),
        c("MATR3", "RAET1E-AS1", "SPATA13"))

    ## up to 10 cells per donor and cell type
    expect_identical(sampleColumn(x), "donor_id")
    donors <- c("T101", "T120", "T126", "T137", "T153", "T154", "T164",
        "T165", "T166", "T167", "T85", "T90")
    cells <- table(x$donor_id, x$cell_type)
    expect_identical(rownames(cells), donors)
    expect_identical(colnames(cells),
        c("Basal", "Multiciliated lineage", "Secretory"))
    expect_identical(cells[["T154", "Multiciliated lineage"]], 2L)
    expect_identical(sum(cells == 10L), 35L)
    expect_true(all(endsWith(x$hlca_donor_id, x$donor_id)))
    expect_setequal(x$assay, c("10x 3' v2", "10x 3' v3"))

    ## donor variables
    exposures <- exposureData(x)
    expect_identical(rownames(exposures), donors)
    expect_identical(colnames(exposures),
        c("smoking", "chemistry_v3", "age", "sex_male"))
    expect_identical(unname(exposures[, "smoking"]),
        c(1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1))
    expect_identical(unname(exposures[, "chemistry_v3"]),
        c(1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 1))
    expect_identical(unname(exposures[, "age"]),
        c(55, 57, 35, 27, 38, 61, 66, 64, 68, 66, 59, 44))
    expect_identical(unname(exposures[, "sex_male"]),
        c(1, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1))
    expect_identical(exposureInfo(x)$exposure, colnames(exposures))
    ## the chemistry of a donor is that of each of its cells
    expect_identical(unname(exposures[x$donor_id, "chemistry_v3"]),
        as.numeric(x$assay == "10x 3' v3"))

    ## provenance
    source <- S4Vectors::metadata(x)$source
    expect_identical(source$url, paste0(
        "https://datasets.cellxgene.cziscience.com/",
        "688185ad-11c2-4172-a53a-f4f1f4076860.h5ad"))
    expect_match(source$selection,
        "up to 10 cells per donor and cell type, seed 20261030",
        fixed = TRUE)
    ## the eligible cells of the atlas, of which up to 10 were drawn
    eligible <- source$eligible_cells
    expect_identical(dimnames(eligible), list(donors, colnames(cells)))
    expect_identical(sum(eligible), 29217L)
    expect_identical(eligible["T154", ],
        c(Basal = 1413L, "Multiciliated lineage" = 2L, Secretory = 210L))
    expect_identical(as.vector(cells), as.vector(pmin(eligible, 10L)))
    expect_identical(source$dataset_version,
        "688185ad-11c2-4172-a53a-f4f1f4076860")
    expect_identical(source$bytes, 5873612847)
    expect_identical(source$sha256, paste0("1cbdde1e513a31bdb6a3f70296e2d061",
        "2760d2af4032e239c94596ea87b44007"))
    expect_identical(source$licence, "CC BY 4.0")
    expect_identical(source$references,
        c("doi:10.1038/s41591-023-02327-2", "doi:10.1038/s41467-020-16239-z"))
    expect_identical(source$script, "inst/scripts/make_hlca_airway_subset.R")
})

test_that("the smoking model on the HLCA subset gives the recorded result", {
    ## the result of the default backend when the data were made
    ## (edgeR 4.99.13), with room for differences between platforms
    x <- .hlca_subset()
    fit <- run_sc_exwas(x, exposure = "smoking", covariates = "chemistry_v3")
    expect_identical(nrow(fit), 32890L)
    by_type <- split(as.data.frame(fit), fit$celltype)
    expect_identical(names(by_type),
        c("Basal", "Multiciliated lineage", "Secretory"))
    expect_identical(unname(vapply(by_type, nrow, 0L)),
        c(10708L, 11139L, 11043L))
    ## donor T154 has too few multiciliated-lineage cells
    expect_identical(unname(vapply(by_type, function(d) d$n_donors[1], 0L)),
        c(12L, 11L, 12L))
    expect_identical(unname(vapply(by_type, function(d) d$n_exposed[1], 0L)),
        c(6L, 5L, 6L))
    expect_equal(fit$padj_global, stats::p.adjust(fit$pvalue, "BH"))
    cyp1b1 <- by_type$Basal[by_type$Basal$gene == "ENSG00000138061", ]
    expect_identical(nrow(cyp1b1), 1L)
    expect_equal(cyp1b1$log2FC, 5.974974, tolerance = 1e-4)
    ## on the log scale: a tolerance is absolute for a value as small as
    ## this p-value
    expect_equal(log10(cyp1b1$pvalue), log10(2.862663e-06), tolerance = 1e-4)
    expect_equal(cyp1b1$padj_global, 0.02353825, tolerance = 1e-3)
    design <- S4Vectors::metadata(fit)$diagnostics$Basal$design_columns
    expect_identical(design, c("(Intercept)", "exposure", "chemistry_v3"))
})


## Runs an installed script as Rscript would with these arguments, in an
## environment of its own in which commandArgs() returns them
.run_hlca_script <- function(script, args) {
    env <- new.env(parent = globalenv())
    env$commandArgs <- function(...) args
    source(system.file("scripts", script, package = "exposomeSC"),
        local = env)
    invisible(env)
}

test_that("the HLCA scripts are installed and check their arguments", {
    scripts <- c("make_hlca_airway_subset.R", "fetch_hlca.R",
        "extract_hlca_design_metadata.R", "extract_hlca_feature_metadata.R",
        "create_hlca_input_manifest.R")
    ## .Rbuildignore keeps an allow-list of the files in inst/scripts
    found <- system.file("scripts", scripts, package = "exposomeSC")
    expect_identical(basename(found), scripts)
    installed <- list.files(system.file("scripts", package = "exposomeSC"))
    expect_setequal(installed, c(scripts, "make_incov_pbmc_subset.R"))
    ## and one of the files in inst/extdata: the two data sets and no more
    expect_setequal(
        list.files(system.file("extdata", package = "exposomeSC")),
        c("hlca_airway_subset.rds", "incov_pbmc_subset.rds"))
    for (script in scripts) {
        code <- readLines(system.file("scripts", script,
            package = "exposomeSC"), warn = FALSE)
        ## no path derived from the place of the script in a clone, and no
        ## clone loaded in place of the installed package
        expect_identical(grep("data-raw|project_root|PROJECT_ROOT|load_all",
            code, value = TRUE), character(), info = script)
        ## without arguments a script says how it is called
        expect_error(.run_hlca_script(script, character()),
            paste0("Usage: Rscript ", script, " <"), fixed = TRUE,
            info = script)
    }

    ## the four scripts that read the H5AD name a path that does not exist,
    ## then an output directory that does not exist, then refuse a file of
    ## another size, before they need a package
    folder <- tempfile("hlca_scripts_")
    dir.create(folder)
    on.exit(unlink(folder, recursive = TRUE), add = TRUE)
    absent <- file.path(folder, "absent.h5ad")
    small <- file.path(folder, "small.h5ad")
    writeLines("not the atlas", small)
    nowhere <- file.path(folder, "absent", "out.rds")
    for (script in scripts[1:4]) {
        expect_error(
            .run_hlca_script(script, c(absent, file.path(folder, "out.rds"))),
            paste0("not found: ", absent), fixed = TRUE, info = script)
        expect_error(
            .run_hlca_script(script, c(small, nowhere)),
            paste0("The directory of the output file does not exist: ",
                dirname(nowhere)), fixed = TRUE, info = script)
        expect_error(
            .run_hlca_script(script, c(small, file.path(folder, "out.rds"))),
            "pinned", fixed = TRUE, info = script)
    }
    expect_identical(list.files(folder), "small.h5ad")
})

test_that("the HLCA manifest script lists its inputs by file name", {
    skip_if_not_installed("digest")
    folder <- tempfile("hlca_manifest_")
    dir.create(folder)
    on.exit(unlink(folder, recursive = TRUE), add = TRUE)
    inputs <- file.path(folder, c("a.h5ad", "b.rds", "c.rds", "d.rds"))
    for (k in seq_along(inputs)) {
        writeBin(as.raw(seq_len(10 * k)), inputs[k])
    }
    out <- file.path(folder, "manifest")
    ## the inputs are checked before the output directory is made
    expect_error(
        .run_hlca_script("create_hlca_input_manifest.R",
            c(out, inputs[1:3], file.path(folder, "absent.rds"))),
        paste0("Missing canonical input(s): ",
            file.path(folder, "absent.rds")), fixed = TRUE)
    expect_error(
        .run_hlca_script("create_hlca_input_manifest.R",
            c(out, inputs[1:3], file.path(folder, "other", "a.h5ad"))),
        "The four inputs must have different file names.", fixed = TRUE)
    expect_false(dir.exists(out))

    expect_message(
        .run_hlca_script("create_hlca_input_manifest.R", c(out, inputs)),
        "Saved SHA-256 manifest for 4 canonical inputs.", fixed = TRUE)
    manifest <- readRDS(file.path(out, "hlca_input_manifest.rds"))
    expect_identical(manifest$file, basename(inputs))
    expect_false("relative_path" %in% names(manifest))
    expect_identical(manifest$bytes, c(10, 20, 30, 40))
    expect_identical(manifest$sha256, unname(tools::sha256sum(inputs)))
    expect_identical(
        utils::read.csv(file.path(out, "hlca_input_manifest.csv"))$file,
        basename(inputs))

    ## fetch_hlca.R stops at a manifest that does not exist and at one
    ## without an entry for its H5AD, here a.h5ad with another size; neither
    ## is passed over in silence
    fetch <- function(manifest_path) {
        .run_hlca_script("fetch_hlca.R",
            c(inputs[1], file.path(folder, "scee.rds"), manifest_path))
    }
    absent <- file.path(folder, "absent.rds")
    expect_error(fetch(absent),
        paste0("Input manifest not found: ", absent), fixed = TRUE)
    written <- file.path(out, "hlca_input_manifest.rds")
    no_entry <- paste0("The manifest ", written,
        " has no entry for a.h5ad with the pinned size.")
    expect_error(fetch(written), no_entry, fixed = TRUE)
    ## a manifest of 0.99.1 named its files in the column relative_path
    names(manifest)[names(manifest) == "file"] <- "relative_path"
    saveRDS(manifest, written)
    expect_error(fetch(written), no_entry, fixed = TRUE)
    ## with an entry of the pinned size the manifest is accepted, and the
    ## script goes on to the size of the file itself
    names(manifest)[names(manifest) == "relative_path"] <- "file"
    manifest$bytes[1] <- 5873612847
    saveRDS(manifest, written)
    expect_error(fetch(written),
        "HLCA H5AD byte size differs from the pinned dataset version.",
        fixed = TRUE)
})
