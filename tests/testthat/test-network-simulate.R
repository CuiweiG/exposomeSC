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

# The precision matrices are filled without loops over pairs of nodes. The
# functions with the loops are kept here as the reference: from the same
# state of the generator they must give the same matrix and leave the same
# state behind.

.loop_generate_sparse_precision <- function(p, density, n_tx, n_met,
                                            cross_frac) {
    prec <- diag(1, p)
    n_possible <- p * (p - 1) / 2
    n_edges <- max(1L, floor(n_possible * density))
    n_cross <- max(0L, floor(n_edges * cross_frac))
    n_within <- n_edges - n_cross
    if (n_within > 0) {
        within_pairs <- list()
        idx <- 0L
        for (i in seq_len(max(n_tx - 1L, 0L))) {
            for (j in (i + 1L):n_tx) {
                idx <- idx + 1L
                within_pairs[[idx]] <- c(i, j)
            }
        }
        if (n_met >= 2L) {
            for (i in (n_tx + 1L):(p - 1L)) {
                for (j in (i + 1L):p) {
                    idx <- idx + 1L
                    within_pairs[[idx]] <- c(i, j)
                }
            }
        }
        if (length(within_pairs) > 0) {
            sel <- sample(length(within_pairs),
                min(n_within, length(within_pairs)))
            for (s in sel) {
                pair <- within_pairs[[s]]
                val <- stats::runif(1, 0.2, 0.5) * sample(c(-1, 1), 1)
                prec[pair[1], pair[2]] <- val
                prec[pair[2], pair[1]] <- val
            }
        }
    }
    if (n_cross > 0 && n_tx > 0 && n_met > 0) {
        cross_pairs <- list()
        idx <- 0L
        for (i in seq_len(n_tx)) {
            for (j in (n_tx + 1L):p) {
                idx <- idx + 1L
                cross_pairs[[idx]] <- c(i, j)
            }
        }
        sel <- sample(length(cross_pairs), min(n_cross, length(cross_pairs)))
        for (s in sel) {
            pair <- cross_pairs[[s]]
            val <- stats::runif(1, 0.15, 0.4) * sample(c(-1, 1), 1)
            prec[pair[1], pair[2]] <- val
            prec[pair[2], pair[1]] <- val
        }
    }
    eig <- eigen(prec, symmetric = TRUE, only.values = TRUE)
    if (min(eig$values) < 0.1) {
        prec <- prec + diag(abs(min(eig$values)) + 0.1, p)
    }
    prec
}

.loop_perturb_precision <- function(prec, frac_change) {
    ut <- upper.tri(prec)
    edges <- which(ut & prec != 0, arr.ind = TRUE)
    non_edges <- which(ut & prec == 0, arr.ind = TRUE)
    n_change <- max(1L, floor(nrow(edges) * frac_change))
    if (nrow(edges) > 0 && n_change > 0) {
        to_remove <- sample(nrow(edges), min(n_change, nrow(edges)))
        for (r in to_remove) {
            prec[edges[r, 1], edges[r, 2]] <- 0
            prec[edges[r, 2], edges[r, 1]] <- 0
        }
    }
    if (nrow(non_edges) > 0 && n_change > 0) {
        to_add <- sample(nrow(non_edges), min(n_change, nrow(non_edges)))
        for (a in to_add) {
            val <- stats::runif(1, 0.15, 0.4) * sample(c(-1, 1), 1)
            prec[non_edges[a, 1], non_edges[a, 2]] <- val
            prec[non_edges[a, 2], non_edges[a, 1]] <- val
        }
    }
    eig <- eigen(prec, symmetric = TRUE, only.values = TRUE)
    if (min(eig$values) < 0.1) {
        prec <- prec + diag(abs(min(eig$values)) + 0.1, nrow(prec))
    }
    prec
}

## The value of 'call', which is evaluated only after the seed is set, and
## the next uniform draw after it
.with_next_draw <- function(seed, call) {
    set.seed(seed)
    list(value = call, next_draw = stats::runif(1))
}

test_that("the precision matrices are those of the loops over pairs", {
    sizes <- list(c(1, 1), c(2, 1), c(1, 3), c(2, 2), c(6, 4), c(12, 5),
        c(30, 10))
    settings <- expand.grid(density = c(0.1, 0.5, 1),
        cross_frac = c(0, 0.2, 1))
    for (size in sizes) {
        for (k in seq_len(nrow(settings))) {
            arguments <- list(p = sum(size), density = settings$density[k],
                n_tx = size[1], n_met = size[2],
                cross_frac = settings$cross_frac[k])
            generated <- .with_next_draw(100 + k,
                do.call(.generate_sparse_precision, arguments))
            expect_identical(generated, .with_next_draw(100 + k,
                do.call(.loop_generate_sparse_precision, arguments)))
            for (change in c(0.4, 1)) {
                expect_identical(
                    .with_next_draw(200 + k,
                        .perturb_precision(generated$value, change)),
                    .with_next_draw(200 + k,
                        .loop_perturb_precision(generated$value, change)))
            }
        }
    }
})

test_that("the pairs of nodes come in the order of the loops", {
    expect_identical(.node_pairs(0L), matrix(0L, 0L, 2L))
    expect_identical(.node_pairs(1L), matrix(0L, 0L, 2L))
    expect_identical(.node_pairs(2L), matrix(1:2, 1L, 2L))
    expect_identical(.node_pairs(4L),
        cbind(c(1L, 1L, 1L, 2L, 2L, 3L), c(2L, 3L, 4L, 3L, 4L, 4L)))
    expect_identical(.node_pairs(4), .node_pairs(4L))
})

test_that("a seed gives the networks it gave with the loops over pairs", {
    ## the numbers are those of the version with the loops
    sim <- simulate_crossomic_network(n_donors = 20L, n_celltypes = 4L,
        n_transcripts = 12L, n_metabolites = 5L, seed = 7L)
    n_edges <- vapply(sim$true_adjacency,
        function(a) sum(a[upper.tri(a)] != 0), numeric(1))
    expect_identical(sum(n_edges), 48)
    expect_equal(sum(abs(unlist(sim$true_networks))), 102.25702849949,
        tolerance = 1e-10)
})
