#!/usr/bin/env Rscript
## Build inst/extdata/incov_pbmc_subset.rds, the real multi-omic example of
## exposomeSC: blood cells with surface proteins, plasma metabolites and
## plasma exposure biomarkers of 43 patients of the INCOV cohort.
##
## Usage:
##   Rscript make_incov_pbmc_subset.R <folder> <output.rds>
##
## <folder> holds the 88 input files, under the names of their sources:
##   - Su2020_TableS1.xlsx: "Table S1. Human subject details, plasma
##     proteomic and metabolomic datasets and analysis, and CITE-seq
##     antibodies" of Mendeley Data, doi:10.17632/tzydswhhb5.5 (version 5,
##     published 2020-11-04, CC BY 4.0; 7,296,770 bytes), saved under this
##     name;
##   - E-MTAB-9357.sdrf.txt and, for each of the 43 donors, the files
##     heathlab_dc_9_17_pbmc_gex_library_<n>_1.txt.gz (gene expression) and
##     heathlab_dc_9_17_pbmc_pro_library_<n>_1.txt.gz (surface proteins) of
##     ArrayExpress E-MTAB-9357, from
##     https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/357/E-MTAB-9357/Files/
##     (released in October 2020; the files carry the date 2024-01-18 and
##     were downloaded on 2026-11-19). On 2026-11-21 the authors who
##     deposited the accession gave permission in these words: "The
##     single-cell values from E-MTAB-9357 included in exposomeSC may be
##     used and redistributed under CC BY 4.0, with citation of Su et al.
##     (2020)."
## The script checks the size and the SHA-256 of every file before it reads
## any of them.
##
## Reference: Su Y, Chen D, Yuan D, et al. (2020). Multi-omics resolves a
## sharp disease-state shift between mild and moderate COVID-19. Cell
## 183(6), 1479-1495.e20, doi:10.1016/j.cell.2020.10.037.
##
## Donors (fixed from the clinical sheet and the list of files, before any
## single-cell value of these donors was read)
##   - the first blood draw of a patient (sample "<subject>-BL"; "T1" in the
##     clinical sheet S1.1);
##   - single cells, plasma metabolomics and plasma proteomics all available
##     for that draw;
##   - cigarette smoking recorded as never, former or current;
##   - not in hospital at that draw (WHO ordinal scale "1", "1 or 2" or "2").
##   43 patients meet these conditions.
##
## Counts
##   The archive holds ln(1 + counts per million) of each cell (STAR Methods
##   of the paper), so expm1() of a value is the count times 1e6 / total of
##   the cell, and the smallest non-zero value of a cell is that of its
##   smallest count. This count is 1 in nearly every cell; where it is not,
##   it is the smallest whole number that turns every value of the cell into
##   a whole number. The script stops unless every value is a whole multiple
##   of the unit of its cell (to 1 part in 10,000) and the unit times the
##   total of the cell is 1e6 (to 1 part in 1,000).
##
## Cell types
##   Two gates per sample, each the mean of the two centres that two-means
##   clustering finds in log(1 + count), started at the 10 % and 90 %
##   quantiles: one for the protein CD3 and one for the sum of the proteins
##   CD14, CD11c, CD33 and CD64.
##   - "T cell": CD3 above its gate, the sum at or below its gate, and fewer
##     than 10 counts of the genes LYZ, S100A8 and S100A9 together.
##   - "Monocyte": the sum above its gate, CD3 at or below its gate, at
##     least 5 counts of LYZ, S100A8 and S100A9 and fewer than 3 counts of
##     CD3D, CD3E and TRAC together.
##   All other cells are counted as "other". The numbers of cells of the
##   three kinds are kept in metadata(x)$source$gated_cells.
##
## Cells
##   Up to 14 monocytes and 10 T cells per donor, drawn at random without
##   looking at any other value: donors and cell types are taken in the
##   order of their names (sorted as in the C locale), the cells of a
##   stratum in the order of the file, and sample.int() draws from them after
##   set.seed(20261119, kind = "Mersenne-Twister",
##            normal.kind = "Inversion", sample.kind = "Rejection").
##   The numbers 14 and 10 were set by the size limit of a source package.
##   Genes without a count in the selected cells are left out; the column
##   sums are those of all genes of the archive file.
##
## Donor variables (exposureData of the result)
##   ever_smoker  1 = current or former smoker, 0 = never ("Cigarette
##                Smoking" of sheet S1.1)
##   age          years; sex_male  1 = male, 0 = female
##   and 13 plasma compounds of exogenous origin from sheet S1.4 (see
##   `biomarkers` below): relative abundances, scaled by the authors to a
##   median of 1 within each batch, NA where the sheet has no value for the
##   sample. The sheet does not say why a value is missing.
##   exposureInfo(x)$lowest_reported holds the smallest value of each
##   compound among the 387 samples of the sheet. For the metabolites that
##   they analysed, those with fewer than 20 % of the values missing, the
##   authors put this smallest value in place of a missing one (STAR
##   Methods); six of the 13 compounds have more missing values than that.
##
## Plasma metabolites (metadata(x)$plasma_metabolites)
##   The metabolites of sheet S1.4 that sheet S1.6 assigns to a super-pathway
##   other than "Xenobiotics" and that have a value for all 43 donors, with
##   their pathways in metadata(x)$metabolite_info.
##
## Packages: data.table with R.utils, which it needs to read the compressed
## files, Matrix, S4Vectors, SummarizedExperiment, SingleCellExperiment,
## exposomeSC. The xlsx file is read with base R.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
    stop("Usage: Rscript make_incov_pbmc_subset.R <folder> <output.rds>")
}
folder <- args[[1]]
out_file <- args[[2]]
## The paths are checked before any package is loaded or any data are read
if (!dir.exists(folder)) stop("Folder not found: ", folder)
if (!dir.exists(dirname(out_file))) {
    stop("The directory of the output file does not exist: ",
         dirname(out_file))
}

## ---- the input files: name, bytes, SHA-256 ------------------------------
## g = gene expression, p = surface proteins; the number is that of the file
## name, heathlab_dc_9_17_pbmc_<gex|pro>_library_<number>_1.txt.gz
sample_files <- utils::read.table(text = "
g 39 4500481 9cb6ff3834df49856b3f503f083630254398a7f7ed61eeb7fbe8fc1262ef9f5d
g 42 4178718 c34e0f335315fa541ccc6108b620bcee3902d3f2588ab2b110b863c55ddd5770
g 56 8791207 5674314d4ec3148f59555d7c43069fb8d737f59a0908f6ee953ba8d77112fb63
g 57 3320977 04acff6e3f4beaf00534633fbb5a5e9602e95434f94744cff354d254f33f775e
g 58 4239019 23258e0346749e7a8e1bcd539b2299253d0b4e0eac37affcbc240f149ec787fb
g 60 3037071 1e5f0bddbd862e30198bac72db2a94c69686748f3c8a1784556037a3782b537f
g 61 2796982 1659400a9ca443ccbcbf7486474463d8c88bbc64a7a872171cf78288719ca813
g 63 3034929 63e28e4f923936569826b95f64c2fc60afdbc75fad434272cdf3d88ab0916e5a
g 66 8990879 64df4cd2027e59f866bb497b3e0581062c3d0ad9eb52750260d60fee6289d27d
g 67 3169511 0fbbdb9bb8cf365a736b11b03fc746b334ef5a8ece837711a304addf8035c76c
g 68 4314965 f836f35f12228b4451a0415ffda09ecf94544aff6eea78fd42ec668753b02291
g 74 1803489 b7a57ec6c91711e085021b2c107ccbbc0cc3d75cad04b294bc0921cb8bdc0d8c
g 76 2399004 324cf5bf57d2def7d7bb15dff3c875ced043d078141c41020773b268bbcaa0c9
g 77 8947020 c752d90fa392eb2d3d990da0aa726d7873b88540c1846cc8d1acb3dbd79faa33
g 78 2595346 24d98515e49794d73dd43d3f27ddb67b283ff5c1fb88abb5ba0fabf4c4f6dbd7
g 81 2539939 657e576c69efbb20393f70f07e36053c6ed618a51285210635aa905d887e9c90
g 82 9967503 07820ad55dffc4ffd88edaafaee6b99e96e4b7a1e60b32a5581f3b4573c89044
g 83 7893558 39b6c930ff1c0b6e6cd0b0021f514afd028a884fee63c931429f94ef9f0e23f4
g 88 1649742 65530ec064a7ef78d9772b86b9de377295c7f7fe0f5e44e0b1daa5ad0ec6979a
g 92 2765480 6b3ea635bdb1e47d531e652c14ed44932a8170d11b995376b886f66633be15b9
g 94 2601115 b68a2d90d7a511fad08ef89dd354af621b501008cdb706131494f9181dfb74aa
g 97 1001358 9939c1bc6d8f52a75a645f356d57b101d314939708f7ca3dbb601cc32694b3d2
g 100 2195721 7e9651daed2de5d8b391ca48e88504708c996d2fbc0998510ff460ef0495f95c
g 101 2600139 037d92924cf3306b0eef914f9cbd0a7c579c613e7ca238467f555c6c5c8c29d3
g 102 2398496 8eafd8fff1949516d6a96ec74a9d3f5227a06a1687810b0aa3f8b1637581ef14
g 103 1356407 477e3195829838aff2f04c5160016e2d45eae930b8ff7ea27cd825401047cb62
g 107 1545794 4493fd2e731c4175472a398fbda95a1ae3aab01b7701212b09692fa9e20ae6a2
g 108 1164799 f96cb0ad762f71687573f4cf321ef75adcc4b2da74c491543468759f0a9fa873
g 109 593737 63e62ac619f939a1bf3cfe57f7bb165dd880ebc8e35d7c2e52f4102d4f097b4c
g 110 2153666 39024675f5726a6b84631e3ff8e2e9ed0b36f2cd2d40496361aae6e6f374aed4
g 113 3237197 9e71f94af7d0697c45123cc97d882fa392ff063e1abb7de863be768fdcf8af15
g 114 1582398 c3c8f03331b508edd3e9414687b31923df69ae0279ddf177f9821a50323749aa
g 116 4958563 7b270ae952d1851cc78b40d99a7e89a65299b0add5daa47ce2e6cd5f2bd70229
g 120 5470940 378edf6090a6cf60e76663b096774f2ec246c0a458596d10ed299bd05dbb128a
g 122 4050479 c5d64bd7f66c2a110958ef1de733321332a1cf42e83d5be9b6141bd2d1ad41b6
g 123 5128713 f1f740ecd174a5c7a913a4f698ce4616ccc25d0c3c0844a48fe65395a2ce68b1
g 126 4869947 96ce4cf545aa2b3d01ea9f75926d0190c085efd7ec7a2f86c3ef0c1e8e925c4d
g 127 1270927 49b35c6e2b76642f1034d0c54d201432690816dfc0a4067e9ba7f362deffb720
g 128 1599072 a635af087eda843292840deeacff3e669870cb01fb391c1327ebf79ee268224f
g 131 1842095 271d3dcbc6223ac0d997a7e360efba55a0ab272274022011ad37b42f476a5de2
g 132 675547 4ebfce25ef6c001eab9e96a18ce9288285861b34dd098cb155f10a48f2e5d157
g 134 1132575 ff8ef6f86dcae19621c1aaaa9fdf80647281ce7b5a4a636d255b7e621392489d
g 135 1227061 4ae7018d7cf86ba4c8fcca4a31828a476ea85ba95896d13132195e9d9c858497
p 39 582330 75a9f21fb5dc8697b51798eb46e6bbf90c50581e5f2c3768bf4901805fb88364
p 42 459301 e14bda35a090e2e35c74c12622768b11f14ad8e9f4a619a0657e664121b21a2d
p 56 946956 9428a5eeddf4548c93a3718d3a5c4b854fe6f1e42891910ad7a3d55b4bbd5e45
p 57 350345 3dc0ce612a4ba46b5dd9837d2ef7cd505e82f2f8c1cf1307c1906cf625e87efe
p 58 418129 1c3a814ec374b1a4e910ccde02c2debf0a9c63b9fa928be8b8dee010459f5afb
p 60 311455 a7e23728ca6060351d2c75c539d8a2097383096094ddbe51575046d872ba942c
p 61 275211 1faeb6edaa32743acd2ac3e06a380739a3f0aef5b63f82e5b781dd1dc04b71a5
p 63 316697 d40aa2d03183be5c043fe9864cbd26af2134cf27bd4181437037077e07c1bc24
p 66 1096690 e2637159379f82ee1d0f3e162e4d863e2d39e075861f0c4417dd0145df697bf4
p 67 376187 f06888edfdaf538178932ffe8f47047e20dfcfffa87b30068c06ce9ea2f6e925
p 68 518894 6c271a43f15633377d280c766f343f99680e913b7081062edc4e78c9ccdf4900
p 74 310479 c7e7b1cd754e5254b8906415b8e5bcecfdf10c5e829f7817a82d489e935ff454
p 76 316132 13c73aca8b30156f9df72dbb89718ee963456f136307aeb9ace79b3683437bb7
p 77 1148393 d6d7710122f0c5b74269a8b09efb2094b53fb6e44821fd3a34f07bd514a553e9
p 78 425430 abc844e216696cfb6d41e6debdad36863d2726c294cc45528dcf59c887de1001
p 81 310115 070ce630d19bd9d3c891fe15122a0b4e0594b52763e277c955a02cab4775c26d
p 82 1212029 e63553db15959abfdaa36d6e3ecf3c227a8b426c670a9596f201aaffb5593051
p 83 1023076 c6d14bbbf4239a15c3a471c8180ccc300c6e2c3ad8a5d12117b4a72d8326cb4c
p 88 260591 5d09a63a831f7c89dfe2eb1c2046b5f239f08899b5d7233ce7b279423b03fa7d
p 92 450765 20a7148696796f3d6ed323b7d965d8023ce2869a92e1ff13808ce9ed2f50a150
p 94 439431 1c67f1018427639a20249eac1ad834480febe318d5792824f478ef11b57f2dd1
p 97 161361 0363839d198d91f73938fe3c9dd6e479bd25fa4dfbfaba7ee7c72784db1c4d67
p 100 369578 0dfdcb80c7694f043fcf6e30e2cb80dac414b84f83282e9956fe2f3f5e324094
p 101 423875 74ddc8b1feb3141f733c0cbcda6a3a52d3a6848d7ea03ca07b03ad0999b243e0
p 102 388293 0404b970073eeba79b1360ee368d12559ba4e371f02ebede3c05f962f2351f0a
p 103 216033 698f28413055e755d395075dd3448bf89cf5ff67c21f598ea3b40d085c1d1c7f
p 107 249626 f8c28251a87ec074576c685cfcbb7d473626a5e66eacbb63cf37e857c7e4321d
p 108 174123 0a24de37f40e267c9813141fd5fe7025ee4dc277fe14505acc99273e24e8a474
p 109 86235 85db513f8306cf3215678a6f5f2f96da82b82529dc97eadf3b25e5bc87c6d4f8
p 110 348022 16e17ffcf0ab2c0682a9be32d8e3aa0c8c2848f45850d817970929e4ce9bce0d
p 113 591875 4d9c603c745683c325c047d031559497395c5346e889e7c09f2a91ed43819b8a
p 114 266101 65a4e6986df4fa301a5116e94e9b79f224dae81b337bac072fe59fcfa700a335
p 116 912827 00db40bbd6e0282ab8dd43ca2999f20eb47ad8d3243a57874d9db5ab78b09396
p 120 955333 1006b862c813005302b6d766d1285d796e64a711c2b63a8c9779029cd36fedae
p 122 756826 43e55edeecdd70029057fb74ccff7304267d9607342aa11c083b8c46c28cc355
p 123 950459 dc66161d92fef4e7822b3618c49a8d3c6bf1e03525a47a0505d3d6f50504f23e
p 126 890320 c8ca95160f322cd3efc16856b20e0a030f04bac1aaa8f3963a6d8547b8ba131c
p 127 211260 4a9cc6a929e6a0c094bc9c586a3862b2708391676538b659993428e23f64c0a7
p 128 268845 e8e39110ac6823b0c25194e7716e672c345c529b82759347d2b19cc4356ba275
p 131 312937 fad0a2c0d76446edf3091e05e9b3aa627c4c28843d2ab02fb7ad5189769c06b0
p 132 111106 bfbee20533ed19202e5d87fce8f8256aa9077400e14be6c77df4657cf5a05f41
p 134 192568 931b363c392e04871dfe27c7c5d8d0c7d7d4380e07a943ea7ccf3522f08ce372
p 135 204287 00a33842bb795202a75d84673d2ea56948f2c3bc09bfa9e65494f89a1fafc109
", col.names = c("kind", "number", "bytes", "sha256"),
    colClasses = c("character", "integer", "numeric", "character"))
sample_files$file <- sprintf("heathlab_dc_9_17_pbmc_%s_library_%d_1.txt.gz",
    ifelse(sample_files$kind == "g", "gex", "pro"), sample_files$number)
inputs <- rbind(
    data.frame(
        file = c("Su2020_TableS1.xlsx", "E-MTAB-9357.sdrf.txt"),
        bytes = c(7296770, 351598),
        sha256 = c(
            paste0("1082603fe92429d2d115c6ee9e2ed822",
                   "3f8ce0699a95669347fca8307dbc7744"),
            paste0("a3ed7ce93c21f708fb22cd09519d11af",
                   "1fe473a74e5e3c07368f34cffc959c5a")),
        stringsAsFactors = FALSE),
    sample_files[, c("file", "bytes", "sha256")])
paths <- file.path(folder, inputs$file)
absent <- inputs$file[!file.exists(paths)]
if (length(absent)) {
    stop(length(absent), " input file(s) not found in ", folder, ", the ",
         "first: ", absent[[1]])
}
wrong_size <- inputs$file[file.size(paths) != inputs$bytes]
if (length(wrong_size)) {
    stop("File(s) of another size than the pinned one: ",
         paste(wrong_size, collapse = ", "))
}
## data.table::fread() reads the compressed files through R.utils
for (pkg in c("data.table", "R.utils", "Matrix", "S4Vectors",
              "SummarizedExperiment", "SingleCellExperiment",
              "exposomeSC")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        stop("Package '", pkg, "' is required.")
    }
}
wrong_hash <- inputs$file[unname(tools::sha256sum(paths)) != inputs$sha256]
if (length(wrong_hash)) {
    stop("File(s) whose SHA-256 is not the pinned one: ",
         paste(wrong_hash, collapse = ", "))
}

## ---- reading a sheet of an xlsx file with base R ------------------------
## An xlsx file is a zip archive of XML files. The cells of a sheet are
## returned as a character matrix.
xml_unescape <- function(x) {
    x <- gsub("&lt;", "<", x, fixed = TRUE)
    x <- gsub("&gt;", ">", x, fixed = TRUE)
    x <- gsub("&quot;", "\"", x, fixed = TRUE)
    x <- gsub("&apos;", "'", x, fixed = TRUE)
    gsub("&amp;", "&", x, fixed = TRUE)
}
read_whole <- function(path) readChar(path, file.size(path), useBytes = TRUE)
read_xlsx_sheet <- function(file, sheet) {
    work <- tempfile("xlsx_")
    dir.create(work)
    on.exit(unlink(work, recursive = TRUE), add = TRUE)
    utils::unzip(file, files = c("xl/workbook.xml",
        "xl/_rels/workbook.xml.rels", "xl/sharedStrings.xml"), exdir = work)
    workbook <- read_whole(file.path(work, "xl/workbook.xml"))
    sheets <- regmatches(workbook, gregexpr("<sheet [^>]*>", workbook))[[1]]
    relations <- read_whole(file.path(work, "xl/_rels/workbook.xml.rels"))
    relations <- regmatches(relations,
        gregexpr("<Relationship [^>]*>", relations))[[1]]
    target <- sub('.*Target="([^"]*)".*', "\\1", relations)
    names(target) <- sub('.*Id="([^"]*)".*', "\\1", relations)
    sheet_names <- xml_unescape(sub('.*name="([^"]*)".*', "\\1", sheets))
    wanted <- match(sheet, sheet_names)
    if (is.na(wanted)) stop("No sheet '", sheet, "' in ", file)
    part <- target[[sub('.*r:id="([^"]*)".*', "\\1", sheets[[wanted]])]]
    utils::unzip(file, files = file.path("xl", part), exdir = work)
    strings <- read_whole(file.path(work, "xl/sharedStrings.xml"))
    items <- regmatches(strings,
        gregexpr("<si>.*?</si>", strings, perl = TRUE))[[1]]
    shared <- vapply(regmatches(items, gregexpr("<t[^>]*>[^<]*</t>", items)),
        function(pieces) paste(gsub("<[^>]*>", "", pieces), collapse = ""),
        character(1))
    shared <- xml_unescape(shared)
    xml <- read_whole(file.path(work, "xl", part))
    cells <- regmatches(xml, gregexpr(
        "<c r=\"[A-Z]+[0-9]+\"[^>]*?(/>|>.*?</c>)", xml, perl = TRUE))[[1]]
    reference <- sub("^<c r=\"([A-Z]+[0-9]+)\".*", "\\1", cells)
    type <- ifelse(grepl("^<c [^>]* t=\"", cells),
        sub("^<c [^>]* t=\"([^\"]*)\".*", "\\1", cells), "n")
    value <- ifelse(grepl("<v>", cells, fixed = TRUE),
        sub(".*<v>([^<]*)</v>.*", "\\1", cells), NA_character_)
    inline <- type == "inlineStr"
    value[inline] <- gsub("<[^>]*>", "",
        sub(".*<is>(.*)</is>.*", "\\1", cells[inline]))
    is_shared <- type == "s" & !is.na(value)
    value[is_shared] <- shared[as.integer(value[is_shared]) + 1L]
    text <- type %in% c("str", "inlineStr")
    value[text] <- xml_unescape(value[text])
    row <- as.integer(sub("^[A-Z]+", "", reference))
    column <- vapply(strsplit(sub("[0-9]+$", "", reference), ""),
        function(chars) {
            sum(match(chars, LETTERS) * 26^(rev(seq_along(chars)) - 1L))
        }, numeric(1))
    out <- matrix(NA_character_, max(row), max(column))
    out[cbind(row, column)] <- value
    out
}

## ---- the donors ---------------------------------------------------------
table_s1 <- file.path(folder, "Su2020_TableS1.xlsx")
clinical <- read_xlsx_sheet(table_s1, "S1.1 Patient Clinical Data")
colnames(clinical) <- clinical[1, ]
clinical <- as.data.frame(clinical[-1, ], stringsAsFactors = FALSE,
    check.names = FALSE)
clinical[] <- lapply(clinical, trimws)
first <- clinical[clinical[["Blood draw time point"]] == "T1", ]
if (anyDuplicated(first[["Study Subject ID"]]) ||
    !all(first[["Sample ID"]] == paste0(first[["Study Subject ID"]], "-1"))) {
    stop("Sheet S1.1 does not hold one first draw per subject, named ",
         "<subject>-1.")
}
first$sample <- paste0(first[["Study Subject ID"]], "-BL")

sdrf <- utils::read.delim(file.path(folder, "E-MTAB-9357.sdrf.txt"),
    check.names = FALSE, stringsAsFactors = FALSE)
with_files <- function(pattern) {
    column <- which(vapply(sdrf, function(v) any(grepl(pattern, v)),
        logical(1)))[[1]]
    sdrf[[column]]
}
archive <- data.frame(sample = sdrf[["Source Name"]],
    gex = with_files("_gex_library_"), pro = with_files("_pro_library_"),
    stringsAsFactors = FALSE)

metabolite_sheet <- read_xlsx_sheet(table_s1, "S1.4 Plasma Metabolomic Data")
plasma <- metabolite_sheet[-1, -(1:2)]
storage.mode(plasma) <- "numeric"
dimnames(plasma) <- list(trimws(metabolite_sheet[-1, 2]),
    metabolite_sheet[1, -(1:2)])
if (anyDuplicated(rownames(plasma)) || anyDuplicated(colnames(plasma))) {
    stop("Sheet S1.4 has duplicated sample or metabolite names.")
}

eligible <- first[["10X?"]] == "yes" & first[["Metabolomics"]] == "yes" &
    first[["Proteomics?"]] == "yes" &
    first$sample %in% archive$sample & first$sample %in% rownames(plasma) &
    first[["Cigarette Smoking"]] %in% c("Never", "Former", "Current") &
    first[["Who Ordinal Scale"]] %in% c("1", "1 or 2", "2")
if (sum(eligible) != 43L) {
    stop("The rule selects ", sum(eligible), " donors, not 43.")
}
donors_table <- first[eligible, ]
donors_table <- donors_table[order(donors_table[["Study Subject ID"]],
    method = "radix"), ]
donors <- donors_table[["Study Subject ID"]]
files <- archive[match(donors_table$sample, archive$sample), ]
if (!setequal(c(files$gex, files$pro), sample_files$file)) {
    stop("The files of the selected donors are not the pinned ones.")
}

## ---- counts of one archive file -----------------------------------------
read_counts <- function(file, tolerance = 1e-4) {
    values <- data.table::fread(file, header = TRUE, showProgress = FALSE)
    cells <- values[[1]]
    features <- names(values)[-1]
    m <- as.matrix(values[, -1])
    rm(values)
    if (anyDuplicated(cells) || anyDuplicated(features)) {
        stop("Duplicated cell or feature names in ", file)
    }
    position <- which(m != 0, arr.ind = TRUE)
    value <- expm1(m[position])
    rm(m)
    cell <- position[, 1L]
    feature <- position[, 2L]
    by_cell <- factor(cell, levels = seq_along(cells))
    smallest <- as.vector(tapply(value, by_cell, min))
    if (anyNA(smallest)) stop("A cell without a count in ", file)
    off <- function(ratio) abs(ratio - round(ratio)) / pmax(ratio, 1)
    ## the smallest count of a cell: 1, or the smallest whole number that
    ## turns every value of the cell into a whole number
    smallest_count <- rep(1L, length(cells))
    ratio <- value / smallest[cell]
    doubtful <- unique(cell[off(ratio) > tolerance])
    if (length(doubtful)) {
        members <- split(seq_along(cell), by_cell)[doubtful]
        for (k in seq_along(doubtful)) {
            scaled <- value[members[[k]]] / smallest[doubtful[[k]]]
            for (multiple in 2:50) {
                if (max(off(scaled * multiple)) <= tolerance) {
                    smallest_count[doubtful[[k]]] <- multiple
                    break
                }
            }
            if (smallest_count[doubtful[[k]]] == 1L) {
                stop("The counts of a cell cannot be recovered in ", file)
            }
        }
        ratio <- value * smallest_count[cell] / smallest[cell]
    }
    count <- round(ratio)
    total <- as.vector(tapply(count, by_cell, sum))
    scale <- smallest / smallest_count * total
    if (max(off(ratio)) > tolerance || max(abs(scale / 1e6 - 1)) > 1e-3) {
        stop("The values of ", file, " are not ln(1 + counts per million).")
    }
    Matrix::sparseMatrix(i = feature, j = cell, x = count,
        dims = c(length(features), length(cells)),
        dimnames = list(features, cells))
}

## ---- gates --------------------------------------------------------------
## The gate between the two modes of log(1 + count), on the count scale
gate <- function(count) {
    z <- log1p(count)
    centres <- as.vector(stats::quantile(z, c(0.1, 0.9)))
    for (step in seq_len(200L)) {
        high <- abs(z - centres[[2]]) < abs(z - centres[[1]])
        if (!any(high) || all(high)) stop("A gate found one mode only.")
        new <- c(mean(z[!high]), mean(z[high]))
        if (all(abs(new - centres) < 1e-10)) break
        centres <- new
    }
    expm1(mean(centres))
}
monocyte_proteins <- c("CD14-1", "CD11C", "CD33-1", "CD64")
monocyte_genes <- c("LYZ", "S100A8", "S100A9")
t_cell_genes <- c("CD3D", "CD3E", "TRAC")
celltypes <- c("Monocyte", "T cell")
cells_per_donor <- c(Monocyte = 14L, `T cell` = 10L)

## ---- the cells ----------------------------------------------------------
set.seed(20261119, kind = "Mersenne-Twister", normal.kind = "Inversion",
         sample.kind = "Rejection")
rna_parts <- list()
protein_parts <- list()
cell_parts <- list()
gated_cells <- matrix(0L, length(donors), 3L,
    dimnames = list(donors, c(celltypes, "other")))
gates <- matrix(NA_real_, length(donors), 2L,
    dimnames = list(donors, c("CD3", "CD14_CD11c_CD33_CD64")))
genes <- NULL
proteins <- NULL
for (k in seq_along(donors)) {
    rna <- read_counts(file.path(folder, files$gex[[k]]))
    protein <- read_counts(file.path(folder, files$pro[[k]]))
    if (is.null(genes)) {
        genes <- rownames(rna)
        proteins <- rownames(protein)
    }
    if (!identical(colnames(rna), colnames(protein)) ||
        !identical(rownames(rna), genes) ||
        !identical(rownames(protein), proteins)) {
        stop("The files of ", donors[[k]], " do not hold the same cells, ",
             "genes and proteins as the others.")
    }
    cd3 <- protein["CD3", ]
    myeloid <- Matrix::colSums(protein[monocyte_proteins, ])
    gates[k, ] <- c(gate(cd3), gate(myeloid))
    monocyte_rna <- Matrix::colSums(rna[monocyte_genes, ])
    t_cell_rna <- Matrix::colSums(rna[t_cell_genes, ])
    is_t <- cd3 > gates[k, 1] & myeloid <= gates[k, 2] & monocyte_rna < 10
    is_monocyte <- myeloid > gates[k, 2] & cd3 <= gates[k, 1] &
        monocyte_rna >= 5 & t_cell_rna < 3
    label <- ifelse(is_monocyte, "Monocyte", ifelse(is_t, "T cell", "other"))
    gated_cells[k, ] <- as.integer(table(factor(label,
        levels = colnames(gated_cells))))
    for (celltype in celltypes) {
        rows <- which(label == celltype)
        if (!length(rows)) next
        take <- sort(rows[sample.int(length(rows),
            min(cells_per_donor[[celltype]], length(rows)))])
        key <- paste(donors[[k]], celltype)
        rna_parts[[key]] <- rna[, take, drop = FALSE]
        protein_parts[[key]] <- protein[, take, drop = FALSE]
        barcode <- sub("-.*$", "", colnames(rna)[take])
        cell_parts[[key]] <- data.frame(
            cell = paste(donors[[k]], barcode, sep = "_"),
            donor_id = donors[[k]], cell_type = celltype,
            sample_id = donors_table$sample[[k]],
            archive_cell = colnames(rna)[take], stringsAsFactors = FALSE)
    }
    rm(rna, protein)
    invisible(gc())
}
counts <- do.call(cbind, rna_parts)
protein_counts <- do.call(cbind, protein_parts)
cells <- do.call(rbind, cell_parts)
if (anyDuplicated(cells$cell)) stop("Duplicated cell names.")
colnames(counts) <- colnames(protein_counts) <- cells$cell
keep <- Matrix::rowSums(counts) > 0
counts <- counts[keep, ]

## ---- donor variables ----------------------------------------------------
## Plasma compounds of exogenous origin: column of exposureData, name in
## sheet S1.4, family
biomarkers <- utils::read.table(sep = "|", strip.white = TRUE,
    col.names = c("exposure", "metabolite", "family"), text = "
cotinine              | cotinine                         | tobacco
caffeine              | caffeine                         | caffeine
paraxanthine          | paraxanthine                     | caffeine
theobromine           | theobromine                      | caffeine
pfos                  | perfluorooctanesulfonate (PFOS)  | perfluoroalkyl
pfoa                  | perfluorooctanoate (PFOA)        | perfluoroalkyl
methylparaben_sulfate | methyl-4-hydroxybenzoate sulfate | paraben
propylparaben_sulfate | propyl 4-hydroxybenzoate sulfate | paraben
hydroxychlorothalonil | 4-hydroxychlorothalonil          | pesticide
naphthol_sulfate      | 2-naphthol sulfate               | naphthalene
ethyl_glucuronide     | ethyl glucuronide                | alcohol
saccharin             | saccharin                        | sweetener
acesulfame            | acesulfame                       | sweetener
", stringsAsFactors = FALSE)
if (!all(biomarkers$metabolite %in% colnames(plasma))) {
    stop("A plasma compound of the list is not in sheet S1.4.")
}
levels <- plasma[donors_table$sample, biomarkers$metabolite, drop = FALSE]
colnames(levels) <- biomarkers$exposure
exposures <- cbind(
    ever_smoker = as.numeric(
        donors_table[["Cigarette Smoking"]] != "Never"),
    age = as.numeric(donors_table[["Age"]]),
    sex_male = as.numeric(donors_table[["Sex"]] == "Male"),
    levels)
rownames(exposures) <- donors
if (anyNA(exposures[, c("ever_smoker", "age", "sex_male")])) {
    stop("A donor without smoking status, age or sex.")
}
relative <- "relative abundance, median of a batch = 1"
info <- S4Vectors::DataFrame(
    exposure = colnames(exposures),
    family = c("tobacco smoke", "demographic", "demographic",
               biomarkers$family),
    unit = c("1 = current or former smoker, 0 = never smoker", "years",
             "1 = male, 0 = female", rep(relative, nrow(biomarkers))),
    source_name = c("Cigarette Smoking", "Age", "Sex",
                    biomarkers$metabolite),
    lowest_reported = c(rep(NA_real_, 3L), unname(apply(
        plasma[, biomarkers$metabolite, drop = FALSE], 2L, min,
        na.rm = TRUE)))
)

## ---- plasma metabolites -------------------------------------------------
## Sheet S1.6 names the pathways of the metabolites of three comparisons,
## side by side in blocks of 11 columns
changes <- read_xlsx_sheet(table_s1, "S1.6 Pro_Met changes vs WHO")
header <- which(apply(changes, 1L, function(row) {
    any(row == "metabolite_id", na.rm = TRUE)
}))[[1]]
blocks <- which(changes[header, ] == "metabolite_id")
annotation <- unique(do.call(rbind, lapply(blocks, function(start) {
    block <- changes[(header + 1L):nrow(changes), start:(start + 10L),
        drop = FALSE]
    colnames(block) <- changes[header, start:(start + 10L)]
    block <- block[!is.na(block[, "metabolite"]), , drop = FALSE]
    data.frame(metabolite = block[, "metabolite"],
        super_pathway = block[, "superpathway"],
        sub_pathway = block[, "subpathway"], stringsAsFactors = FALSE)
})))
if (anyDuplicated(annotation$metabolite)) {
    stop("A metabolite with two pathways in sheet S1.6.")
}
annotation <- annotation[annotation$super_pathway != "nan" &
    annotation$super_pathway != "Xenobiotics" &
    annotation$metabolite %in% colnames(plasma), ]
donor_plasma <- plasma[donors_table$sample, annotation$metabolite,
    drop = FALSE]
complete <- colSums(is.na(donor_plasma)) == 0L
plasma_metabolites <- donor_plasma[, complete, drop = FALSE]
rownames(plasma_metabolites) <- donors
annotation <- annotation[complete, ]
metabolite_order <- order(annotation$super_pathway, annotation$metabolite,
    method = "radix")
plasma_metabolites <- plasma_metabolites[, metabolite_order, drop = FALSE]
metabolite_info <- S4Vectors::DataFrame(annotation[metabolite_order, ],
    row.names = NULL)

## ---- the SingleCellExposomeExperiment -----------------------------------
sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = counts),
    colData = S4Vectors::DataFrame(
        donor_id = cells$donor_id,
        cell_type = cells$cell_type,
        sample_id = cells$sample_id,
        archive_cell = cells$archive_cell,
        row.names = cells$cell
    ),
    rowData = S4Vectors::DataFrame(symbol = rownames(counts),
        row.names = rownames(counts))
)
SingleCellExperiment::altExp(sce, "CITE") <-
    SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = protein_counts))
scee <- exposomeSC::build_scee(sce, exposures, sample_col = "donor_id",
                               exposure_info = info)
manifest <- inputs
rownames(manifest) <- NULL
S4Vectors::metadata(scee)$source <- list(
    data = paste("INCOV cohort (Su et al. 2020): blood mononuclear cells,",
                 "plasma metabolites and clinical records of 43 patients",
                 "with COVID-19 at their first blood draw"),
    references = "doi:10.1016/j.cell.2020.10.037",
    tables = list(
        doi = "10.17632/tzydswhhb5.5",
        file = "Table S1 (Mendeley Data, version 5, 2020-11-04)",
        licence = "CC BY 4.0"),
    single_cells = list(
        accession = "E-MTAB-9357",
        url = paste0("https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/357/",
                     "E-MTAB-9357/Files/"),
        released = paste("2020-10-17 (public release date of the IDF",
                         "file), 2020-10-29 (BioStudies record)"),
        downloaded = "2026-11-19",
        licence = "CC BY 4.0",
        permission = paste("2026-11-21, from the authors who deposited the",
                           "accession: \"The single-cell values from",
                           "E-MTAB-9357 included in exposomeSC may be used",
                           "and redistributed under CC BY 4.0, with",
                           "citation of Su et al. (2020).\"")),
    files = manifest,
    selection = paste("first blood draw; single cells, plasma metabolomics",
                      "and proteomics available; smoking status recorded;",
                      "WHO ordinal scale 1 or 2; up to 14 monocytes and",
                      "10 T cells per donor, seed 20261119"),
    gated_cells = gated_cells,
    gates = gates,
    script = "inst/scripts/make_incov_pbmc_subset.R"
)
S4Vectors::metadata(scee)$plasma_metabolites <- plasma_metabolites
S4Vectors::metadata(scee)$metabolite_info <- metabolite_info
saveRDS(scee, out_file, compress = "xz")
message(sprintf(
    "%d cells, %d genes, %d proteins, %d donors, %d metabolites; %s bytes %s",
    ncol(scee), nrow(scee), nrow(protein_counts), length(donors),
    ncol(plasma_metabolites), format(file.size(out_file), big.mark = ","),
    paste("written to", out_file)))
