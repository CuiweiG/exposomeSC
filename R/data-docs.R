#' Demo SingleCellExperiment for exposomeSC
#'
#' A simulated SingleCellExperiment with 200 genes, 2000 cells
#' from 20 donors and 3 cell types (Monocyte, NK, T_CD4).
#' Genes 1-10 have a simulated PM2.5-dependent expression
#' increase in Monocytes.
#'
#' @format An RDS file at \code{inst/extdata/demo_sce.rds}.
#'   Load with \code{readRDS(system.file("extdata",
#'   "demo_sce.rds", package = "exposomeSC"))}.
#'
#' @source Simulated using \code{inst/scripts/create_demo_data.R}.
#'   Parameters based on typical PBMC scRNA-seq studies.
#'
#' @name demo_data
#' @docType data
#' @keywords datasets
NULL

#' Demo exposure matrix for exposomeSC
#'
#' A CSV file with 20 donors (rows) and 5 environmental
#' exposures (PM2.5, NO2, Pb, Cd, BPA). Values simulate
#' realistic urban air pollution and biomonitoring data.
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
