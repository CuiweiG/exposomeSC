# tests/testthat/test-bridge-spatial.R

test_that("as_scee reads the ExposomeSet exposure matrix", {
    skip_if_not_installed("Biobase")
    setClass("ExposomeSetForTest", contains = "eSet",
             where = environment())
    samples <- paste0("D", 1:5)
    exposures <- matrix(c(1:5, 11:15), nrow = 2, byrow = TRUE,
        dimnames = list(c("PM25", "NO2"), samples))
    es <- methods::new("ExposomeSetForTest",
        assayData = Biobase::assayDataNew("environment", exp = exposures),
        phenoData = Biobase::AnnotatedDataFrame(
            data.frame(row.names = samples, sex = rep(1:2, length.out = 5))),
        featureData = Biobase::AnnotatedDataFrame(
            data.frame(row.names = c("PM25", "NO2"), family = c("air", "air"))))
    counts <- matrix(1L, nrow = 3, ncol = 10,
        dimnames = list(paste0("G", 1:3), paste0("c", 1:10)))
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = counts),
        colData = S4Vectors::DataFrame(donor_id = rep(samples, each = 2)))
    scee <- suppressWarnings(as_scee(es, sce, sample_col = "donor_id",
                                     exposures = "NO2"))
    expect_equal(unname(exposureData(scee)[samples, "NO2"]), 11:15)
})

test_that("seurat_to_exposure warns on within-donor variation", {
    meta <- data.frame(donor = rep(c("D1", "D2"), each = 3),
                       pm = c(1, 1, 1, 2, 2, 3),
                       sex = rep(c("m", "f"), each = 3))
    expect_warning(m <- seurat_to_exposure(meta, "donor", "pm"),
                   "vary within a donor")
    expect_equal(unname(m[, "pm"]), c(1, 7 / 3))
    expect_error(seurat_to_exposure(meta, "donor", "sex"), "Non-numeric")
})

test_that("run_spatial_exwas refuses unregistered pooled k-means regions", {
    skip_if_not_installed("SpatialExperiment")
    set.seed(1)
    donors <- paste0("D", 1:6)
    spot_donor <- rep(donors, each = 40)
    region <- rep(rep(c("epithelium", "stroma"), each = 20), times = 6)
    counts <- matrix(stats::rpois(200 * length(spot_donor), 5), nrow = 200,
        dimnames = list(paste0("G", 1:200), paste0("s", seq_along(spot_donor))))
    spe <- SpatialExperiment::SpatialExperiment(
        assays = list(counts = counts),
        colData = S4Vectors::DataFrame(sample_id = spot_donor, region = region),
        spatialCoords = cbind(x = stats::runif(240), y = stats::runif(240)))
    exposure <- matrix(stats::rnorm(6), ncol = 1,
        dimnames = list(donors, "PM2.5"))
    expect_error(run_spatial_exwas(spe, exposure, "PM2.5"), "registered")
    res <- run_spatial_exwas(spe, exposure, "PM2.5", region_col = "region",
                             target_genes = paste0("G", 1:10))
    expect_setequal(unique(res$region), c("epithelium", "stroma"))
    expect_true(all(res$gene %in% paste0("G", 1:10)))
    expect_true("padj_global" %in% colnames(res))
})
