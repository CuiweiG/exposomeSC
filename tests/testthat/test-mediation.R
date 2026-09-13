# tests/testthat/test-mediation.R

test_that("run_mediation works", {
    skip_if_not_installed("mediation")

    set.seed(42)
    counts <- matrix(rpois(10000, 10), nrow = 100,
        dimnames = list(paste0("G", 1:100), paste0("c", 1:100)))
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = counts),
        colData = S4Vectors::DataFrame(
            cell_id = paste0("c", 1:100),
            donor_id = rep(paste0("D", 1:10), each = 10),
            cell_type = sample(c("Mono", "NK", "T"), 100,
                replace = TRUE,
                prob = c(0.4, 0.3, 0.3))))
    exp_mat <- matrix(rnorm(20), nrow = 10,
        dimnames = list(paste0("D", 1:10), c("E1", "E2")))

    for (f in list.files("../../R", full.names = TRUE))
        source(f, local = TRUE)

    scee <- build_scee(sce, exp_mat, sample_col = "donor_id")

    med <- run_mediation(scee, exposure = "E1",
        celltype_col = "cell_type", sample_col = "donor_id",
        mediator_celltype = "Mono",
        outcome_celltype = "NK",
        target_genes = paste0("G", 1:5),
        min_cells = 2L, n_sims = 100L)

    expect_s4_class(med, "DataFrame")
    expect_true(nrow(med) > 0)
    expect_true("ACME" %in% colnames(med))
    expect_true("ADE" %in% colnames(med))
    expect_true("prop_mediated" %in% colnames(med))
    expect_true("total_effect" %in% colnames(med))
})
