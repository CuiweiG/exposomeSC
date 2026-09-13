# tests/testthat/test-network-simulate.R

test_that("simulate_crossomic_network produces valid output", {
    sim <- simulate_crossomic_network(
        n_donors = 30,
        n_celltypes = 3,
        n_transcripts = 10,
        n_metabolites = 5,
        seed = 42)

    expect_type(sim, "list")
    expect_named(sim, c("pseudobulk", "metabolites", "exposure",
                        "composition", "n_cells", "true_networks",
                        "true_adjacency", "celltype_names",
                        "node_info", "params"))

    # Check dimensions
    expect_equal(length(sim$exposure), 30)
    expect_equal(nrow(sim$metabolites), 30)
    expect_equal(ncol(sim$metabolites), 5)
    expect_equal(length(sim$celltype_names), 3)

    # Check pseudobulk per cell type
    expect_equal(length(sim$pseudobulk), 3)
    for (ct in sim$celltype_names) {
        pb <- sim$pseudobulk[[ct]]
        expect_equal(nrow(pb), 30)
        expect_equal(ncol(pb), 10)
        expect_true(all(pb >= 0))  # counts
    }

    # Check true networks
    for (ct in sim$celltype_names) {
        prec <- sim$true_networks[[ct]]
        expect_equal(nrow(prec), 15)  # 10 tx + 5 met
        expect_equal(ncol(prec), 15)
        # Symmetric
        expect_equal(prec, t(prec))
    }

    # Node info
    expect_equal(nrow(sim$node_info), 15)
    expect_true("omic_layer" %in% colnames(sim$node_info))
    expect_equal(
        sum(sim$node_info$omic_layer == "transcript"), 10)
    expect_equal(
        sum(sim$node_info$omic_layer == "metabolite"), 5)

    # Composition with confounding
    expect_equal(nrow(sim$composition), 30)
    expect_equal(ncol(sim$composition), 3)
    # Rows sum to ~1
    expect_true(all(abs(rowSums(sim$composition) - 1) < 0.01))
})

test_that("simulate without composition confounding", {
    sim <- simulate_crossomic_network(
        n_donors = 20,
        n_celltypes = 2,
        n_transcripts = 5,
        n_metabolites = 3,
        composition_confounding = FALSE,
        seed = 123)

    # All compositions should be equal
    expect_true(
        var(sim$composition[, 1]) < 0.001)
})

test_that("true adjacency has expected structure", {
    sim <- simulate_crossomic_network(
        n_donors = 50,
        n_celltypes = 4,
        n_transcripts = 15,
        n_metabolites = 5,
        edge_density = 0.15,
        seed = 99)

    # Cell type 4 has no transcript-transcript or cross-omic edges; its
    # metabolite block is the shared metabolite precision
    adj4 <- sim$true_adjacency[[sim$celltype_names[4]]]
    expect_equal(sum(adj4[1:15, ]), 0)

    # Cell types 1-2 should share edges
    adj1 <- sim$true_adjacency[[sim$celltype_names[1]]]
    adj2 <- sim$true_adjacency[[sim$celltype_names[2]]]
    # Both should have some edges
    expect_true(sum(adj1) > 0)
    expect_true(sum(adj2) > 0)
})

test_that("simulate_crossomic_network runs with the default seed", {
    expect_type(simulate_crossomic_network(n_donors = 10, n_transcripts = 4,
                                           n_metabolites = 2), "list")
})

test_that("simulated data follow the returned precision matrix", {
    sim <- simulate_crossomic_network(n_donors = 4000, n_celltypes = 1,
        n_transcripts = 6, n_metabolites = 3, edge_density = 0.4,
        cross_omic_frac = 0.5, exposure_effect = 0,
        composition_confounding = FALSE, seed = 7)
    ct <- sim$celltype_names[1]
    latent <- log(sim$pseudobulk[[ct]] / sim$n_cells[, ct]) - 1
    ## Poisson noise is small at these depths; compare partial correlations
    joint <- cbind(latent, sim$metabolites)
    estimated <- -stats::cov2cor(solve(stats::cov(joint)))
    truth <- -stats::cov2cor(sim$true_networks[[ct]])
    diag(estimated) <- diag(truth) <- 0
    cross <- truth[1:6, 7:9]
    expect_true(any(cross != 0))
    expect_lt(max(abs(estimated[1:6, 7:9] - cross)), 0.1)
})
