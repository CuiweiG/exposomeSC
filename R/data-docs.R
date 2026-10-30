#' Demo SingleCellExperiment for exposomeSC
#'
#' A simulated SingleCellExperiment with 200 genes and 2000 cells
#' from 20 donors. Counts are independent Poisson draws with mean 8,
#' with no donor-level variation and no overdispersion, and each cell
#' is given one of three labels (Monocyte, NK, T_CD4) at random.
#' Genes 1-10 receive an additive count increase of
#' \code{round(PM2.5 / 10)} in cells labelled Monocyte. The data are
#' for quick illustrative runs; they lack the sparsity and donor
#' effects of real single-cell data, so results on them say nothing
#' about performance on real data.
#'
#' @format An RDS file at \code{inst/extdata/demo_sce.rds}.
#'   Load with \code{readRDS(system.file("extdata",
#'   "demo_sce.rds", package = "exposomeSC"))}.
#'
#' @source Simulated using \code{inst/scripts/create_demo_data.R}.
#'
#' @name demo_data
#' @docType data
#' @keywords datasets
NULL

#' Demo exposure matrix for exposomeSC
#'
#' A CSV file with 20 donors (rows), five simulated exposure
#' variables (PM2.5, NO2, Pb, Cd, BPA) and the covariates age and
#' sex. PM2.5 and NO2 are drawn from normal distributions, so NO2 can
#' be negative, and the metals and BPA from log-normal distributions.
#' The values have no units and are not measurements.
#'
#' @format A CSV file at \code{inst/extdata/demo_exposures.csv}.
#'   Load with \code{read.csv(system.file("extdata",
#'   "demo_exposures.csv", package = "exposomeSC"),
#'   row.names = 1)}.
#'
#' @source Simulated using \code{inst/scripts/create_demo_data.R}.
#'
#' @name demo_data
#' @docType data
#' @keywords datasets
NULL
