# tests/testthat/test-network-compare.R

test_that("run_comparative_network works with fisher_z", {
    library(S4Vectors)

    # Create two mock CelltypeNetworkResult objects
    p <- 6
    feat <- c("g1", "g2", "g3", "m1", "m2", "m3")
    ni <- DataFrame(
        feature = feat,
        omic_layer = c(rep("transcript", 3),
                       rep("metabolite", 3)),
        block = c(rep(1L, 3), rep(2L, 3)))

    # Network 1: 3 edges
    prec1 <- diag(1, p)
    prec1[1, 4] <- prec1[4, 1] <- 0.3  # cross-omic
    prec1[1, 2] <- prec1[2, 1] <- 0.4  # within tx
    prec1[5, 6] <- prec1[6, 5] <- 0.2  # within met
    adj1 <- (abs(prec1) > 0.1) * 1L
    diag(adj1) <- 0

    net1 <- new("CelltypeNetworkResult",
        precision_matrix = prec1,
        adjacency_matrix = adj1,
        stability_scores = matrix(nrow = 0, ncol = 0),
        node_info = ni,
        celltype = "Mono",
        method = "block_glasso",
        metadata = list(n_donors = 50L))

    # Network 2: different edges
    prec2 <- diag(1, p)
    prec2[2, 5] <- prec2[5, 2] <- 0.5  # different cross-omic
    prec2[1, 2] <- prec2[2, 1] <- 0.4  # shared within tx
    adj2 <- (abs(prec2) > 0.1) * 1L
    diag(adj2) <- 0

    net2 <- new("CelltypeNetworkResult",
        precision_matrix = prec2,
        adjacency_matrix = adj2,
        stability_scores = matrix(nrow = 0, ncol = 0),
        node_info = ni,
        celltype = "NK",
        method = "block_glasso",
        metadata = list(n_donors = 50L))

    comp <- run_comparative_network(net1, net2,
                                     method = "fisher_z")

    expect_s4_class(comp, "NetworkComparison")
    expect_true(nrow(comp@shared_edges) > 0 ||
                nrow(comp@diff_edges) > 0)
    expect_equal(comp@summary$celltypes, c("Mono", "NK"))

    # Show method
    expect_output(show(comp), "NetworkComparison")
})

test_that("run_comparative_network rejects < 2 networks", {
    ni <- S4Vectors::DataFrame(
        feature = "g1", omic_layer = "transcript",
        block = 1L)
    net1 <- new("CelltypeNetworkResult",
        precision_matrix = matrix(1, 1, 1),
        adjacency_matrix = matrix(0, 1, 1),
        stability_scores = matrix(nrow = 0, ncol = 0),
        node_info = ni, celltype = "A",
        method = "test", metadata = list())

    expect_error(run_comparative_network(net1),
                 "At least 2")
})

test_that("NetworkComparison S4 class works", {
    nc <- new("NetworkComparison")
    expect_s4_class(nc, "NetworkComparison")
})
