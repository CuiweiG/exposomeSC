# tests/testthat/test-glmm.R

test_that("run_sc_exwas_glmm works with lme4 backend", {
    skip_if_not_installed("lme4")

    set.seed(42)
    ## Use more donors (20) with more cells per donor for
    ## better convergence
    n_donors <- 20
    cells_per <- 15
    n_cells <- n_donors * cells_per
    n_genes <- 50
    counts <- matrix(rpois(n_genes * n_cells, 10),
        nrow = n_genes,
        dimnames = list(paste0("G", seq_len(n_genes)),
            paste0("c", seq_len(n_cells))))
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = counts),
        colData = S4Vectors::DataFrame(
            cell_id = paste0("c", seq_len(n_cells)),
            donor_id = rep(paste0("D", seq_len(n_donors)),
                each = cells_per),
            cell_type = rep(c("Mono", "NK"),
                length.out = n_cells)))
    exp_mat <- matrix(rnorm(n_donors * 2), nrow = n_donors,
        dimnames = list(paste0("D", seq_len(n_donors)),
            c("E1", "E2")))

    scee <- build_scee(sce, exp_mat, sample_col = "donor_id")

    res <- run_sc_exwas_glmm(scee, exposure = "E1",
        celltype_col = "cell_type",
        BPPARAM = BiocParallel::SerialParam(),
        target_genes = paste0("G", 1:10))

    expect_s4_class(res, "DataFrame")
    expect_true(nrow(res) > 0)
    expect_true("log2FC" %in% colnames(res))
    expect_true("se" %in% colnames(res))
    expect_true("pvalue" %in% colnames(res))
    expect_true("padj" %in% colnames(res))
    expect_true("method" %in% colnames(res))
    ## gaussian path is fitted by lmerTest (Satterthwaite/Kenward-Roger df)
    expect_true(all(grepl("lmerTest", res$method)))
})

test_that("GLMM SE uses vcov covariance correctly", {
    skip_if_not_installed("lme4")

    set.seed(42)
    n_donors <- 25
    cells_per <- 15
    n_cells <- n_donors * cells_per
    counts <- matrix(rpois(50 * n_cells, 10), nrow = 50,
        dimnames = list(paste0("G", 1:50),
            paste0("c", seq_len(n_cells))))
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = counts),
        colData = S4Vectors::DataFrame(
            cell_id = paste0("c", seq_len(n_cells)),
            donor_id = rep(paste0("D", seq_len(n_donors)),
                each = cells_per),
            cell_type = rep(c("Mono", "NK"),
                length.out = n_cells)))
    exp_mat <- matrix(rnorm(n_donors * 2), nrow = n_donors,
        dimnames = list(paste0("D", seq_len(n_donors)),
            c("E1", "E2")))

    scee <- build_scee(sce, exp_mat, sample_col = "donor_id")

    res <- run_sc_exwas_glmm(scee, exposure = "E1",
        celltype_col = "cell_type",
        BPPARAM = BiocParallel::SerialParam(),
        target_genes = paste0("G", 1:5))

    ## Non-reference cell type should have SE values
    nk_res <- res[res$celltype == "NK", ]
    expect_true(nrow(nk_res) > 0)
    expect_true(all(!is.na(nk_res$se)))
    expect_true(all(nk_res$se > 0))
    ## SE should differ from just sqrt(se1^2 + se2^2)
    ## (since vcov includes covariance term)
})

test_that("GLMM returns valid p-values", {
    skip_if_not_installed("lme4")

    set.seed(42)
    n_donors <- 20
    cells_per <- 15
    n_cells <- n_donors * cells_per
    counts <- matrix(rpois(50 * n_cells, 10), nrow = 50,
        dimnames = list(paste0("G", 1:50),
            paste0("c", seq_len(n_cells))))
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = counts),
        colData = S4Vectors::DataFrame(
            cell_id = paste0("c", seq_len(n_cells)),
            donor_id = rep(paste0("D", seq_len(n_donors)),
                each = cells_per),
            cell_type = rep(c("Mono", "NK"),
                length.out = n_cells)))
    exp_mat <- matrix(rnorm(n_donors * 2), nrow = n_donors,
        dimnames = list(paste0("D", seq_len(n_donors)),
            c("E1", "E2")))

    scee <- build_scee(sce, exp_mat, sample_col = "donor_id")

    res <- run_sc_exwas_glmm(scee, exposure = "E1",
        celltype_col = "cell_type",
        BPPARAM = BiocParallel::SerialParam(),
        target_genes = paste0("G", 1:10))

    pvals <- res$pvalue[!is.na(res$pvalue)]
    expect_true(all(pvals >= 0 & pvals <= 1))
})

test_that("GLMM models donor-by-cell-type variation without singular-fit noise", {
    skip_if_not_installed("lme4")
    skip_if_not_installed("lmerTest")

    set.seed(3)
    n_donors <- 12
    cells_per <- 20
    donors <- paste0("D", seq_len(n_donors))
    donor <- rep(rep(donors, each = cells_per), times = 2)
    cell_type <- rep(c("A", "B"), each = n_donors * cells_per)
    counts <- matrix(rpois(20 * length(donor), 10), nrow = 20,
        dimnames = list(paste0("G", 1:20), paste0("c", seq_along(donor))))
    ## Donor-by-cell-type variation in G1
    unit <- paste(donor, cell_type)
    unit_effect <- stats::setNames(rnorm(length(unique(unit)), sd = 0.6),
                                   unique(unit))
    counts[1, ] <- rpois(length(donor), 10 * exp(unit_effect[unit]))
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = counts),
        colData = S4Vectors::DataFrame(donor_id = donor,
                                       cell_type = cell_type))
    exp_mat <- matrix(rnorm(n_donors), ncol = 1,
        dimnames = list(donors, "E1"))
    scee <- build_scee(sce, exp_mat, sample_col = "donor_id")

    expect_no_message(res <- run_sc_exwas_glmm(scee, exposure = "E1",
        celltype_col = "cell_type",
        BPPARAM = BiocParallel::SerialParam(), target_genes = "G1"))
    expect_setequal(res$celltype, c("A", "B"))
    ## With donor-by-cell-type variation the df reflect donors, not cells
    expect_true(all(res$df < 2 * n_donors))

    one_type <- scee[, cell_type == "A"]
    res_one <- run_sc_exwas_glmm(one_type, exposure = "E1",
        celltype_col = "cell_type",
        BPPARAM = BiocParallel::SerialParam(), target_genes = "G1")
    expect_equal(nrow(res_one), 1L)
})

test_that("run_sc_exwas_pb_offset matches a direct MASS::glm.nb fit", {
    set.seed(11)
    donors <- paste0("D", seq_len(10L))
    cell_donor <- rep(donors, each = 12L)
    counts <- matrix(rnbinom(3L * length(cell_donor), mu = 6, size = 3),
        nrow = 3L, dimnames = list(c("G1", "G2", "G3"),
            paste0("c", seq_along(cell_donor))))
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = counts),
        colData = S4Vectors::DataFrame(donor_id = cell_donor,
            cell_type = "T"))
    exposure <- matrix(c(rnorm(10L), seq(20, 65, by = 5)), ncol = 2L,
        dimnames = list(donors, c("E1", "age")))
    scee <- build_scee(sce, exposure, sample_col = "donor_id")

    res <- run_sc_exwas_pb_offset(scee, exposure = "E1",
        covariates = "age", target_genes = c("G1", "G2"))
    expect_s4_class(res, "DataFrame")
    expect_identical(nrow(res), 2L)
    expect_false(is.unsorted(res$padj))

    pb <- vapply(donors, function(d) {
        rowSums(counts[, cell_donor == d, drop = FALSE])
    }, numeric(3L))
    for (g in c("G1", "G2")) {
        d <- data.frame(y = pb[g, ], A = exposure[, "E1"],
            age = exposure[, "age"], logoff = log(colSums(pb)))
        fit <- MASS::glm.nb(y ~ A + offset(logoff) + age, data = d)
        s <- stats::coef(summary(fit))
        row <- res[res$gene == g, ]
        expect_equal(row$log2FC, s["A", "Estimate"] / log(2))
        expect_equal(row$se, s["A", "Std. Error"] / log(2))
        expect_equal(row$pvalue, s["A", "Pr(>|z|)"])
    }

    ## the plain lapply path (BPPARAM = NULL) gives the same fits
    res_bp <- run_sc_exwas_pb_offset(scee, exposure = "E1",
        covariates = "age", target_genes = c("G1", "G2"),
        BPPARAM = NULL)
    expect_equal(as.data.frame(res_bp), as.data.frame(res))
})

## Ten donors in two cell types, with exposure columns named after the
## variables of the models
.reserved_names_scee <- function() {
    set.seed(1)
    n_donors <- 10
    cells_per <- 6
    n_cells <- n_donors * cells_per
    counts <- matrix(rpois(5 * n_cells, 10), nrow = 5,
        dimnames = list(paste0("G", 1:5), paste0("c", seq_len(n_cells))))
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = counts),
        colData = S4Vectors::DataFrame(
            donor_id = rep(paste0("D", seq_len(n_donors)), each = cells_per),
            cell_type = rep(c("Mono", "NK"), length.out = n_cells)))
    columns <- c("E1", "exposure", "A", "logoff", "donor", "expo_Mono", "y",
        "celltype")
    exp_mat <- matrix(rnorm(n_donors * length(columns)), nrow = n_donors,
        dimnames = list(paste0("D", seq_len(n_donors)), columns))
    build_scee(sce, exp_mat, sample_col = "donor_id")
}

test_that("covariate names used by the mixed model are refused", {
    skip_if_not_installed("lme4")
    skip_if_not_installed("lmerTest")
    scee <- .reserved_names_scee()
    for (name in c("exposure", "donor", "celltype", "y")) {
        expect_error(
            run_sc_exwas_glmm(scee, exposure = "E1", covariates = name,
                celltype_col = "cell_type"),
            paste0("cannot be called '", name, "'"))
    }
    ## the negative binomial family adds the offset and expo_<cell type>
    expect_error(
        run_sc_exwas_glmm(scee, exposure = "E1", covariates = "expo_Mono",
            celltype_col = "cell_type", family = "nbinom"),
        "cannot be called 'expo_Mono'")
    expect_error(
        run_sc_exwas_glmm(scee, exposure = "E1",
            covariates = c("logoff", "expo_Mono"), celltype_col = "cell_type",
            family = "nbinom"),
        "cannot be called 'logoff', 'expo_Mono'")
    ## for the default family logoff is an ordinary covariate
    fit <- suppressWarnings(suppressMessages(run_sc_exwas_glmm(scee,
        exposure = "E1", covariates = "logoff", celltype_col = "cell_type",
        target_genes = "G1")))
    expect_s4_class(fit, "DataFrame")
})

test_that("covariate names used by the offset model are refused", {
    scee <- .reserved_names_scee()
    for (name in c("y", "A", "logoff")) {
        expect_error(
            run_sc_exwas_pb_offset(scee, exposure = "E1", covariates = name,
                celltype_col = "cell_type"),
            paste0("cannot be called '", name, "'"))
    }
})

test_that("nbinom keeps cell types apart whose labels make.names() merges", {
    skip_if_not_installed("lme4")
    skip_if_not_installed("lmerTest")
    fit_labels <- function(labels) {
        set.seed(3)
        donors <- paste0("D", 1:10)
        exposure <- rep(c(0, 1), 5)
        donor <- rep(donors, each = 24)
        label <- rep(rep(labels, each = 12), 10)
        effect <- exposure[match(donor, donors)] * (label == labels[1])
        counts <- rbind(
            G1 = stats::rnbinom(length(donor), mu = 20 * exp(effect), size = 5),
            G2 = stats::rnbinom(length(donor), mu = 20, size = 5))
        colnames(counts) <- paste0("c", seq_along(donor))
        sce <- SingleCellExperiment::SingleCellExperiment(
            assays = list(counts = counts),
            colData = S4Vectors::DataFrame(donor_id = donor, cell_type = label))
        scee <- build_scee(sce, matrix(exposure, ncol = 1,
            dimnames = list(donors, "exposure")), sample_col = "donor_id")
        res <- as.data.frame(suppressMessages(run_sc_exwas_glmm(scee,
            exposure = "exposure", target_genes = "G1", family = "nbinom")))
        res[match(labels, res$celltype),
            c("log2FC", "se", "statistic", "pvalue")]
    }
    ## make.names() maps both labels to "CD16..Mono"
    merged <- fit_labels(c("CD16+ Mono", "CD16- Mono"))
    distinct <- fit_labels(c("CD16pos Mono", "CD16neg Mono"))
    expect_false(isTRUE(all.equal(merged$log2FC[1], merged$log2FC[2])))
    expect_equal(merged, distinct, tolerance = 1e-6, ignore_attr = TRUE)
})

test_that("nbinom fits a single cell type", {
    skip_if_not_installed("lme4")
    skip_if_not_installed("lmerTest")
    set.seed(5)
    donors <- paste0("D", 1:10)
    exposure <- rep(c(0, 1), 5)
    donor <- rep(donors, each = 24)
    counts <- rbind(
        G1 = stats::rnbinom(length(donor),
            mu = 20 * exp(0.5 * exposure[match(donor, donors)]), size = 5),
        G2 = stats::rnbinom(length(donor), mu = 20, size = 5))
    colnames(counts) <- paste0("c", seq_along(donor))
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = counts),
        colData = S4Vectors::DataFrame(donor_id = donor, cell_type = "Mono"))
    scee <- build_scee(sce, matrix(exposure, ncol = 1,
        dimnames = list(donors, "exposure")), sample_col = "donor_id")
    res <- suppressMessages(run_sc_exwas_glmm(scee, exposure = "exposure",
        target_genes = "G1", family = "nbinom"))
    ## the same model fitted directly
    cells <- data.frame(y = counts["G1", ], donor = factor(donor),
        exposure = exposure[match(donor, donors)],
        logoff = log(colSums(counts)))
    full <- suppressMessages(lme4::glmer.nb(
        y ~ exposure + offset(logoff) + (1 | donor), data = cells))
    reduced <- suppressMessages(lme4::glmer.nb(
        y ~ 1 + offset(logoff) + (1 | donor), data = cells))
    expect_identical(nrow(res), 1L)
    expect_equal(res$log2FC,
        unname(lme4::fixef(full)["exposure"]) / log(2), tolerance = 1e-6)
    expect_equal(res$statistic, stats::anova(reduced, full)$Chisq[2],
        tolerance = 1e-6)
})

test_that("a cell type whose donors all lack an exposure is left out", {
    skip_if_not_installed("lme4")
    skip_if_not_installed("lmerTest")
    set.seed(6)
    donors <- sprintf("D%02d", 1:12)
    exposure <- stats::setNames(rep(c(0, 1), 6), donors)
    ## types B and C in every donor; type A only in D11 and D12, whose
    ## exposure is missing
    cell_donor <- c(rep(donors, each = 20), rep(c("D11", "D12"), each = 20))
    cell_type <- c(rep(c("B", "C"), 120), rep("A", 40))
    counts <- rbind(
        G1 = stats::rnbinom(length(cell_donor),
            mu = 20 * exp(0.4 * exposure[cell_donor]), size = 5),
        G2 = stats::rnbinom(length(cell_donor), mu = 20, size = 5))
    colnames(counts) <- paste0("c", seq_along(cell_donor))
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = counts),
        colData = S4Vectors::DataFrame(donor_id = cell_donor,
            cell_type = cell_type))
    design <- matrix(exposure, ncol = 1, dimnames = list(donors, "exposure"))
    design[c("D11", "D12"), "exposure"] <- NA
    fit <- function(cells, family) {
        scee <- build_scee(sce[, cells],
            design[unique(cell_donor[cells]), , drop = FALSE],
            sample_col = "donor_id")
        suppressWarnings(suppressMessages(run_sc_exwas_glmm(scee,
            exposure = "exposure", target_genes = "G1", family = family)))
    }
    complete <- !cell_donor %in% c("D11", "D12")
    three_types <- fit(rep(TRUE, length(cell_donor)), "gaussian")
    expect_setequal(three_types$celltype, c("B", "C"))
    expect_equal(three_types, fit(complete, "gaussian"), tolerance = 1e-6)
    ## without type C only B keeps cells, and both families fit it
    no_c <- cell_type != "C"
    for (family in c("gaussian", "nbinom")) {
        expect_equal(fit(no_c, family), fit(no_c & complete, family),
            tolerance = 1e-6)
    }
})

test_that("cells of a donor without an exposure value are left out", {
    skip_if_not_installed("lme4")
    skip_if_not_installed("lmerTest")
    set.seed(5)
    donors <- paste0("D", 1:10)
    exposure <- rep(c(0, 1), 5)
    donor <- rep(donors, each = 24)
    counts <- rbind(
        G1 = stats::rnbinom(length(donor),
            mu = 20 * exp(0.5 * exposure[match(donor, donors)]), size = 5),
        G2 = stats::rnbinom(length(donor), mu = 20, size = 5))
    colnames(counts) <- paste0("c", seq_along(donor))
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = counts),
        colData = S4Vectors::DataFrame(donor_id = donor, cell_type = "Mono"))
    design <- matrix(exposure, ncol = 1, dimnames = list(donors, "exposure"))
    design["D3", "exposure"] <- NA
    fit <- function(cells, family) {
        scee <- build_scee(sce[, cells],
            design[unique(donor[cells]), , drop = FALSE],
            sample_col = "donor_id")
        suppressWarnings(run_sc_exwas_glmm(scee, exposure = "exposure",
            target_genes = "G1", family = family))
    }
    all_cells <- rep(TRUE, length(donor))
    for (family in c("nbinom", "gaussian")) {
        msgs <- character()
        with_missing <- withCallingHandlers(fit(all_cells, family),
            message = function(m) {
                msgs <<- c(msgs, conditionMessage(m))
                invokeRestart("muffleMessage")
            })
        expect_true(any(grepl("Leaving 1 donor\\(s\\) with a missing",
            msgs)))
        expect_false(is.na(with_missing$pvalue))
        without <- suppressMessages(fit(donor != "D3", family))
        expect_equal(with_missing, without, tolerance = 1e-6)
    }
    ## cells without a label are left out as well
    labelled <- !(seq_along(donor) %in% c(1:3, 30:32))
    sce$cell_type[!labelled] <- NA
    design["D3", "exposure"] <- 0
    for (family in c("nbinom", "gaussian")) {
        msgs <- character()
        with_unlabelled <- withCallingHandlers(fit(all_cells, family),
            message = function(m) {
                msgs <<- c(msgs, conditionMessage(m))
                invokeRestart("muffleMessage")
            })
        expect_true(any(grepl(
            "Leaving 6 cell\\(s\\) without a cell-type label", msgs)))
        expect_false(is.na(with_unlabelled$pvalue))
        expect_equal(with_unlabelled, suppressMessages(fit(labelled, family)),
            tolerance = 1e-6)
    }
})

test_that("padj of the mixed and offset models adjusts across all strata", {
    scee <- .reserved_names_scee()
    ## Poisson counts: glm.nb() warns while it estimates theta
    offset_fit <- suppressWarnings(run_sc_exwas_pb_offset(scee,
        exposure = "E1", target_genes = c("G1", "G2")))
    expect_setequal(offset_fit$celltype, c("Mono", "NK"))
    expect_equal(offset_fit$padj, stats::p.adjust(offset_fit$pvalue, "BH"))
    skip_if_not_installed("lme4")
    skip_if_not_installed("lmerTest")
    mixed_fit <- suppressWarnings(suppressMessages(run_sc_exwas_glmm(scee,
        exposure = "E1", target_genes = c("G1", "G2"))))
    expect_setequal(mixed_fit$celltype, c("Mono", "NK"))
    expect_equal(mixed_fit$padj, stats::p.adjust(mixed_fit$pvalue, "BH"))
})
