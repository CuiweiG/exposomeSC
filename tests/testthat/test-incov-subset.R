.incov_subset <- function() {
    readRDS(system.file("extdata", "incov_pbmc_subset.rds",
        package = "exposomeSC"))
}

test_that("the INCOV subset is a valid object with the documented content", {
    x <- .incov_subset()
    expect_s4_class(x, "SingleCellExposomeExperiment")
    expect_true(methods::validObject(x))
    expect_identical(dim(x), c(15406L, 1032L))
    expect_identical(SummarizedExperiment::assayNames(x), "counts")
    counts <- SummarizedExperiment::assay(x, "counts")
    expect_s4_class(counts, "dgCMatrix")
    ## raw counts: whole numbers, and no gene without a count
    expect_true(all(counts@x == round(counts@x)))
    expect_true(all(counts@x > 0))
    expect_identical(sum(counts@x), 2874868)
    expect_true(all(Matrix::rowSums(counts) > 0))
    expect_identical(SummarizedExperiment::rowData(x)$symbol, rownames(x))

    ## 14 monocytes and 10 T cells of each of 43 donors, first blood draw
    expect_identical(sampleColumn(x), "donor_id")
    cells <- table(x$donor_id, x$cell_type)
    expect_identical(dim(cells), c(43L, 2L))
    expect_identical(colnames(cells), c("Monocyte", "T cell"))
    expect_true(all(cells[, "Monocyte"] == 14L))
    expect_true(all(cells[, "T cell"] == 10L))
    expect_identical(x$sample_id, paste0(x$donor_id, "-BL"))
    ## a cell is named by its donor and its barcode, and carries its name
    ## in the files of the archive
    expect_false(anyDuplicated(x$archive_cell) > 0L)
    expect_identical(colnames(x),
        paste(x$donor_id, sub("-.*$", "", x$archive_cell), sep = "_"))

    ## the surface proteins of the same cells
    expect_identical(SingleCellExperiment::altExpNames(x), "CITE")
    protein <- SummarizedExperiment::assay(
        SingleCellExperiment::altExp(x, "CITE"), "counts")
    expect_identical(dim(protein), c(192L, 1032L))
    expect_identical(colnames(protein), colnames(x))
    expect_true(all(protein@x == round(protein@x)))
    expect_identical(sum(protein@x), 1652054)
    expect_identical(sum(grepl("Isotype-Control", rownames(protein))), 4L)

    ## the cells follow the rule that assigned their cell types: the gates
    ## of their sample for CD3 and for four monocyte proteins, and counts
    ## of three genes of each kind
    gates <- S4Vectors::metadata(x)$source$gates
    expect_identical(rownames(gates), rownames(exposureData(x)))
    is_monocyte <- x$cell_type == "Monocyte"
    cd3 <- protein["CD3", ]
    myeloid <- Matrix::colSums(
        protein[c("CD14-1", "CD11C", "CD33-1", "CD64"), ])
    expect_true(all((cd3 > gates[x$donor_id, 1]) == !is_monocyte))
    expect_true(all((myeloid > gates[x$donor_id, 2]) == is_monocyte))
    monocyte_rna <- Matrix::colSums(counts[c("LYZ", "S100A8", "S100A9"), ])
    t_cell_rna <- Matrix::colSums(counts[c("CD3D", "CD3E", "TRAC"), ])
    expect_true(all(monocyte_rna[is_monocyte] >= 5))
    expect_true(all(t_cell_rna[is_monocyte] < 3))
    expect_true(all(monocyte_rna[!is_monocyte] < 10))

    ## donor variables: smoking, age, sex and 13 plasma compounds
    exposures <- exposureData(x)
    expect_identical(colnames(exposures),
        c("ever_smoker", "age", "sex_male", "cotinine", "caffeine",
            "paraxanthine", "theobromine", "pfos", "pfoa",
            "methylparaben_sulfate", "propylparaben_sulfate",
            "hydroxychlorothalonil", "naphthol_sulfate", "ethyl_glucuronide",
            "saccharin", "acesulfame"))
    expect_identical(nrow(exposures), 43L)
    expect_identical(rownames(exposures)[c(1, 43)], c("INCOV039", "INCOV135"))
    expect_false(is.unsorted(rownames(exposures)))
    expect_identical(sum(exposures[, "ever_smoker"]), 16)
    expect_identical(sum(exposures[, "sex_male"]), 6)
    expect_identical(range(exposures[, "age"]), c(19, 81))
    ## compounds without a value for a donor
    expect_identical(unname(colSums(is.na(exposures))),
        c(0, 0, 0, 38, 3, 3, 0, 1, 1, 0, 26, 0, 3, 8, 16, 24))
    expect_identical(rownames(exposures)[!is.na(exposures[, "cotinine"])],
        c("INCOV081", "INCOV108", "INCOV120", "INCOV122", "INCOV127"))
    expect_equal(unname(exposures[1:3, "pfos"]),
        c(0.530135011, 2.543051293, 0.479798764), tolerance = 1e-9)
    info <- exposureInfo(x)
    expect_identical(info$exposure, colnames(exposures))
    expect_identical(info$source_name[c(4, 8, 10)], c("cotinine",
        "perfluorooctanesulfonate (PFOS)", "methyl-4-hydroxybenzoate sulfate"))
    expect_identical(info$family, c("tobacco smoke", "demographic",
        "demographic", "tobacco", "caffeine", "caffeine", "caffeine",
        "perfluoroalkyl", "perfluoroalkyl", "paraben", "paraben",
        "pesticide", "naphthalene", "alcohol", "sweetener", "sweetener"))
    expect_identical(unique(info$unit[4:16]),
        "relative abundance, median of a batch = 1")
    ## no reported value is below the smallest one of its compound among
    ## all samples of the study
    lowest <- info$lowest_reported
    expect_true(all(is.na(lowest[1:3])))
    for (k in 4:16) {
        expect_true(all(exposures[, k] >= lowest[k], na.rm = TRUE),
            info = colnames(exposures)[k])
    }
    expect_equal(lowest[4:5], c(0.016246075, 0.018811849), tolerance = 1e-8)

    ## plasma metabolites with a value for every donor, and their pathways
    plasma <- S4Vectors::metadata(x)$plasma_metabolites
    expect_identical(dim(plasma), c(43L, 182L))
    expect_identical(rownames(plasma), rownames(exposures))
    expect_false(anyNA(plasma))
    expect_true(all(plasma > 0))
    pathways <- S4Vectors::metadata(x)$metabolite_info
    expect_identical(pathways$metabolite, colnames(plasma))
    expect_identical(as.vector(table(pathways$super_pathway)),
        c(35L, 6L, 9L, 1L, 109L, 12L, 6L, 4L))
    expect_false("Xenobiotics" %in% pathways$super_pathway)
    expect_false(anyNA(pathways$sub_pathway))
    expect_equal(plasma[["INCOV039", "kynurenine"]], 1.628158693,
        tolerance = 1e-9)

    ## provenance
    source <- S4Vectors::metadata(x)$source
    expect_identical(source$references, "doi:10.1016/j.cell.2020.10.037")
    expect_identical(source$tables$doi, "10.17632/tzydswhhb5.5")
    expect_identical(source$tables$licence, "CC BY 4.0")
    expect_identical(source$single_cells$accession, "E-MTAB-9357")
    expect_match(source$single_cells$licence, "none stated by the authors",
        fixed = TRUE)
    expect_match(source$selection, "seed 20261119", fixed = TRUE)
    expect_identical(colSums(source$gated_cells),
        c(Monocyte = 9889, `T cell` = 38213, other = 15763))
    expect_identical(source$script, "inst/scripts/make_incov_pbmc_subset.R")

    ## the input files of the record are those that the installed script
    ## pins: the two tables, then a file of genes and one of proteins for
    ## each donor, with the same sizes and SHA-256
    files <- source$files
    expect_identical(nrow(files), 88L)
    expect_identical(sum(files$bytes), 173831460)
    expect_identical(files$file[1:2],
        c("Su2020_TableS1.xlsx", "E-MTAB-9357.sdrf.txt"))
    expect_identical(files$bytes[1:2], c(7296770, 351598))
    expect_identical(files$sha256[1:2], c(
        "1082603fe92429d2d115c6ee9e2ed8223f8ce0699a95669347fca8307dbc7744",
        "a3ed7ce93c21f708fb22cd09519d11af1fe473a74e5e3c07368f34cffc959c5a"))
    script <- readLines(system.file("scripts", "make_incov_pbmc_subset.R",
        package = "exposomeSC"))
    pinned <- utils::read.table(
        text = grep("^[gp] [0-9]+ [0-9]+ [0-9a-f]{64}$", script, value = TRUE),
        col.names = c("kind", "number", "bytes", "sha256"),
        colClasses = c("character", "integer", "numeric", "character"))
    expect_identical(as.vector(table(pinned$kind)), c(43L, 43L))
    expect_identical(files$file[-(1:2)],
        sprintf("heathlab_dc_9_17_pbmc_%s_library_%d_1.txt.gz",
            ifelse(pinned$kind == "g", "gex", "pro"), pinned$number))
    expect_identical(files$bytes[-(1:2)], pinned$bytes)
    expect_identical(files$sha256[-(1:2)], pinned$sha256)
    ## the two tables are pinned in the script as well, each SHA-256 in two
    ## halves
    expect_true(any(grepl("bytes = c(7296770, 351598)", script,
        fixed = TRUE)))
    halves <- c(substr(files$sha256[1:2], 1, 32),
        substr(files$sha256[1:2], 33, 64))
    for (half in halves) {
        expect_true(any(grepl(half, script, fixed = TRUE)), info = half)
    }
    ## a file of genes and one of proteins for each donor, in the order of
    ## the donors
    expect_identical(sprintf("INCOV%03d", pinned$number[pinned$kind == "g"]),
        rownames(exposures))
    expect_identical(sprintf("INCOV%03d", pinned$number[pinned$kind == "p"]),
        rownames(exposures))

    ## the description of the data and the release of the cells
    expect_match(source$data, "43 patients", fixed = TRUE)
    expect_match(source$single_cells$released, "2020-10-29", fixed = TRUE)
    ## INCOV132 has the lowest gate of the monocyte proteins, the others
    ## one of 12 or more
    expect_identical(names(which.min(gates[, 2])), "INCOV132")
    expect_true(sort(gates[, 2])[[2]] >= 12)
    expect_identical(source$gated_cells["INCOV132", ],
        c(Monocyte = 16L, `T cell` = 127L, other = 185L))

    ## the units of the donor variables and the pathways of three
    ## metabolites
    expect_identical(info$unit[1:3], c(
        "1 = current or former smoker, 0 = never smoker", "years",
        "1 = male, 0 = female"))
    expect_identical(
        pathways$sub_pathway[match(c("kynurenine", "arginine", "methionine"),
            pathways$metabolite)],
        c("Tryptophan Metabolism",
            "Urea cycle; Arginine and Proline Metabolism",
            "Methionine, Cysteine, SAM and Taurine Metabolism"))
})

test_that("analyses of the INCOV subset give the recorded results", {
    x <- .incov_subset()
    pseudobulk <- donor_pseudobulk(x)
    expect_identical(names(pseudobulk), c("Monocyte", "T cell"))
    totals <- vapply(c("Monocyte", "T cell"), function(celltype) {
        sum(SummarizedExperiment::assay(pseudobulk[[celltype]], "counts"))
    }, numeric(1))
    expect_equal(unname(totals), c(1626291, 1248577))
    ## the association with age, adjusted for sex, as obtained when the
    ## data were made (edgeR 4.99.13), with room for differences between
    ## platforms
    fit <- run_sc_exwas(x, exposure = "age", covariates = "sex_male")
    expect_identical(as.vector(table(fit$celltype)), c(7120L, 6574L))
    expect_identical(unique(fit$n_donors), 43L)
    expect_identical(fit$gene[1:2], c("FOLR3", "IL32"))
    expect_identical(fit$celltype[1:2], c("Monocyte", "T cell"))
    expect_equal(fit$log2FC[1], 0.06437527, tolerance = 1e-4)
    ## on the log scale: a tolerance is absolute for a value as small as
    ## this p-value
    expect_equal(log10(fit$pvalue[1]), log10(1.334647e-12),
        tolerance = 1e-3)
})

## Runs the installed script as Rscript would with these arguments, in an
## environment of its own in which commandArgs() returns them
.run_incov_script <- function(args) {
    env <- new.env(parent = globalenv())
    env$commandArgs <- function(...) args
    source(system.file("scripts", "make_incov_pbmc_subset.R",
        package = "exposomeSC"), local = env)
    invisible(env)
}

test_that("the INCOV script is installed and checks its arguments", {
    installed <- list.files(system.file("scripts", package = "exposomeSC"))
    expect_true("make_incov_pbmc_subset.R" %in% installed)
    expect_error(.run_incov_script(character()),
        "Usage: Rscript make_incov_pbmc_subset.R <folder> <output.rds>",
        fixed = TRUE)

    folder <- tempfile("incov_script_")
    dir.create(folder)
    on.exit(unlink(folder, recursive = TRUE), add = TRUE)
    absent <- file.path(folder, "absent")
    expect_error(
        .run_incov_script(c(absent, file.path(folder, "out.rds"))),
        paste0("Folder not found: ", absent), fixed = TRUE)
    expect_error(
        .run_incov_script(c(folder, file.path(absent, "out.rds"))),
        "The directory of the output file does not exist: ", fixed = TRUE)
    ## an empty folder: the 88 input files are missed, the table first
    expect_error(
        .run_incov_script(c(folder, file.path(folder, "out.rds"))),
        paste0("88 input file(s) not found in ", folder,
            ", the first: Su2020_TableS1.xlsx"), fixed = TRUE)
    ## the names of the input files are those of the provenance record of
    ## the data; files of another size are refused before any is read
    files <- S4Vectors::metadata(.incov_subset())$source$files$file
    for (file in files) writeLines("not the data", file.path(folder, file))
    expect_error(
        .run_incov_script(c(folder, file.path(folder, "out.rds"))),
        paste0("File(s) of another size than the pinned one: ",
            "Su2020_TableS1.xlsx, E-MTAB-9357.sdrf.txt"), fixed = TRUE)
    expect_setequal(list.files(folder), files)
})
