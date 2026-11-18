# tests/testthat/test-network-compare.R

test_that("run_comparative_network classifies the edges of two networks", {
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

    expect_no_warning(comp <- run_comparative_network(net1, net2))
    ## Edge g1-m1 has precision 0.3 in net1 only: partial correlation -0.3
    edge <- rbind(as.data.frame(diffEdges(comp)),
                  as.data.frame(sharedEdges(comp)))
    edge <- edge[edge$node_i == "g1" & edge$node_j == "m1", ]
    expect_equal(edge$max_diff, 0.3, tolerance = 1e-12)

    expect_s4_class(comp, "NetworkComparison")
    ## one edge in both networks, two in Mono only, one in NK only
    columns <- c("node_i", "node_j", "layer_i", "layer_j", "cross_omic",
        "category", "max_diff")
    shared <- sharedEdges(comp)
    expect_identical(colnames(shared), columns)
    expect_identical(paste(shared$node_i, shared$node_j), "g1 g2")
    expect_equal(shared$max_diff, 0, tolerance = 1e-12)
    differential <- diffEdges(comp)
    expect_identical(colnames(differential), columns)
    expect_identical(paste(differential$node_i, differential$node_j),
        c("g1 m1", "g2 m2", "m2 m3"))
    expect_identical(differential$category,
        c("unique_Mono", "unique_NK", "unique_Mono"))
    expect_identical(differential$cross_omic, c(TRUE, TRUE, FALSE))
    expect_identical(names(celltypeSpecificEdges(comp)), c("Mono", "NK"))
    expect_equal(comparisonSummary(comp)$celltypes, c("Mono", "NK"))
    expect_identical(comparisonSummary(comp)$n_edges, c(Mono = 3L, NK = 2L))
    expect_equal(comparisonSummary(comp)$jaccard$Mono_vs_NK, 1 / 4)

    ## what show() prints: every line has another count
    expect_identical(utils::capture.output(show(comp)),
        c("NetworkComparison", "  Cell types compared: Mono, NK",
            "  Shared edges: 1", "  Differential edges: 3",
            "  Unique to Mono: 2", "  Unique to NK: 1"))
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

test_that("run_comparative_network describes and does not test", {
    feat <- c("g1", "g2", "m1")
    make <- function(ct, w) {
        prec <- diag(3)
        prec[1, 3] <- prec[3, 1] <- w
        adj <- (prec != 0) * 1L
        diag(adj) <- 0L
        new("CelltypeNetworkResult",
            precision_matrix = prec, adjacency_matrix = adj,
            stability_scores = matrix(nrow = 0, ncol = 0),
            node_info = S4Vectors::DataFrame(feature = feat,
                omic_layer = c("transcript", "transcript", "metabolite")),
            celltype = ct, method = "precomputed",
            metadata = list(n_donors = 20L))
    }
    expect_no_warning(
        comp <- run_comparative_network(make("A", -0.4), make("B", -0.1)))
    ## no test: no p-value columns and no method to choose
    expect_identical(names(formals(run_comparative_network)), "...")
    shared <- sharedEdges(comp)
    expect_identical(nrow(shared), 1L)
    expect_false(any(c("p_value", "p_adjusted") %in% colnames(shared)))
    expect_equal(shared$max_diff, 0.3, tolerance = 1e-12)
    expect_identical(names(comparisonSummary(comp)),
        c("celltypes", "n_edges", "jaccard", "n_total_edges", "n_shared",
            "n_differential"))
    ## an argument of the removed test is taken for a network and named
    expect_error(
        run_comparative_network(make("A", -0.4), make("B", -0.1),
            method = "fisher_z"),
        "Argument 3 ('method') is not a CelltypeNetworkResult",
        fixed = TRUE)
    expect_error(
        run_comparative_network(make("A", -0.4), make("B", -0.1),
            adjust = "BH"),
        "Argument 3 ('adjust') is not a CelltypeNetworkResult",
        fixed = TRUE)
    expect_error(run_comparative_network(make("A", -0.4), "B"),
        "Argument 2 is not a CelltypeNetworkResult", fixed = TRUE)
    ## an argument without a name among named ones is not given one
    expect_error(
        run_comparative_network(first = make("A", -0.4), "B",
            third = make("C", -0.1)),
        "Argument 2 is not a CelltypeNetworkResult", fixed = TRUE)
    ## networks without an edge: empty tables with the same columns
    none <- run_comparative_network(make("A", 0), make("B", 0))
    expect_identical(colnames(sharedEdges(none)), colnames(shared))
    expect_identical(nrow(diffEdges(none)), 0L)
    expect_identical(names(comparisonSummary(none)),
        c("celltypes", "n_edges", "jaccard"))
})

test_that("run_comparative_network classifies the edges of three networks", {
    feat <- c("g1", "g2", "m1", "m2")
    make <- function(ct, edges, weights, diagonal = rep(1, 4)) {
        prec <- diag(diagonal)
        adj <- matrix(0L, 4, 4)
        for (k in seq_along(weights)) {
            i <- edges[[k]][1]
            j <- edges[[k]][2]
            prec[i, j] <- prec[j, i] <- weights[k]
            adj[i, j] <- adj[j, i] <- 1L
        }
        dimnames(prec) <- dimnames(adj) <- list(feat, feat)
        CelltypeNetworkResult(prec, adj,
            node_info = S4Vectors::DataFrame(feature = feat,
                omic_layer = c("transcript", "transcript", "metabolite",
                    "metabolite")),
            celltype = ct)
    }
    ## g1-g2 in all three, g1-m1 in A and B only, m1-m2 in C only; the
    ## precision of g1 in A is 4, so that its partial correlations are the
    ## entries divided by 2
    a <- make("A", list(c(1, 2), c(1, 3)), c(-0.8, -0.6),
        diagonal = c(4, 1, 1, 1))
    b <- make("B", list(c(1, 2), c(1, 3)), c(-0.1, 0.2))
    c3 <- make("C", list(c(1, 2), c(3, 4)), c(0.3, -0.5))
    comp <- run_comparative_network(a, b, c3)

    expect_identical(sharedEdges(comp)$category, "shared")
    ## partial correlations of g1-g2: 0.4, 0.1 and -0.3
    expect_equal(sharedEdges(comp)$max_diff, 0.7, tolerance = 1e-12)
    differential <- diffEdges(comp)
    expect_identical(paste(differential$node_i, differential$node_j),
        c("g1 m1", "m1 m2"))
    ## an edge in two of the three networks is differential
    expect_identical(differential$category, c("differential", "unique_C"))
    ## g1-m1: 0.3 in A, -0.2 in B and 0 in C; m1-m2: 0.5 in C only
    expect_equal(differential$max_diff, c(0.5, 0.5), tolerance = 1e-12)
    expect_identical(names(celltypeSpecificEdges(comp)), "C")
    summary <- comparisonSummary(comp)
    expect_identical(summary$celltypes, c("A", "B", "C"))
    expect_identical(summary$n_edges, c(A = 2L, B = 2L, C = 2L))
    expect_equal(summary$jaccard,
        list(A_vs_B = 1, A_vs_C = 1 / 3, B_vs_C = 1 / 3))
    expect_identical(summary$n_shared, 1L)
    expect_identical(summary$n_differential, 2L)
})

test_that("the comparison says how the features of two networks differ", {
    make <- function(ct, feat) {
        p <- length(feat)
        CelltypeNetworkResult(diag(p), matrix(0L, p, p),
            node_info = S4Vectors::DataFrame(feature = feat,
                omic_layer = ifelse(startsWith(feat, "g"), "transcript",
                    "metabolite")),
            celltype = ct)
    }
    first <- make("A", c("g1", "g2", "m1"))
    different <- paste("All networks must have identical features.",
        "Network for 'B' has different features from 'A'.")
    ## a feature that only the second network has is named
    expect_error(run_comparative_network(first,
            make("B", c("g1", "g2", "m1", "m2"))),
        paste(different, "Only in 'B': m2."), fixed = TRUE)
    expect_error(run_comparative_network(first, make("B", c("g1", "m1"))),
        paste(different, "Only in 'A': g2."), fixed = TRUE)
    ## the same features in another order
    in_order <- paste(different,
        "The features are the same, in another order.")
    err <- expect_error(run_comparative_network(first,
        make("B", c("g2", "g1", "m1"))))
    expect_identical(conditionMessage(err), in_order)
    ## the same names, one of them twice: not the same features
    err <- expect_error(run_comparative_network(first,
        make("B", c("g1", "g2", "g2", "m1"))))
    expect_identical(conditionMessage(err), different)
})
