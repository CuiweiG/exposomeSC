## Constructors, accessors and validity of the network result classes

.toy_network <- function(celltype = "Mono", weight = -0.4,
                         stability = NULL) {
    features <- c("G1", "G2", "M1")
    precision <- diag(3)
    precision[1, 3] <- precision[3, 1] <- weight
    adjacency <- (precision != 0) * 1L
    diag(adjacency) <- 0L
    dimnames(precision) <- dimnames(adjacency) <- list(features, features)
    CelltypeNetworkResult(
        precision = precision,
        adjacency = adjacency,
        node_info = data.frame(
            feature = features,
            omic_layer = c("transcript", "transcript", "metabolite")),
        celltype = celltype,
        stability = stability,
        metadata = list(n_donors = 20L))
}

test_that("CelltypeNetworkResult() builds a valid object with accessors", {
    net <- .toy_network()
    expect_s4_class(net, "CelltypeNetworkResult")
    expect_true(methods::validObject(net))
    expect_identical(networkCelltype(net), "Mono")
    expect_identical(networkMethod(net), "precomputed")
    expect_identical(networkPrecision(net)["G1", "M1"], -0.4)
    expect_identical(networkAdjacency(net)["G1", "M1"], 1L)
    expect_identical(networkAdjacency(net)["G1", "G2"], 0L)
    expect_identical(dim(networkStability(net)), c(0L, 0L))
    expect_s4_class(networkNodes(net), "DataFrame")
    expect_identical(networkNodes(net)$feature, c("G1", "G2", "M1"))
    expect_identical(metadata(net), list(n_donors = 20L))

    stab <- networkAdjacency(net) * 0.9
    with_stability <- .toy_network(stability = stab)
    expect_identical(networkStability(with_stability), stab)
})

test_that("CelltypeNetworkResult validity rejects inconsistent parts", {
    good <- .toy_network()
    expect_error(
        CelltypeNetworkResult(networkPrecision(good), networkAdjacency(good),
            node_info = networkNodes(good), celltype = "Mono",
            stability = matrix(0, 2, 2)),
        "stability_scores must be empty")
    expect_error(
        CelltypeNetworkResult(networkPrecision(good), networkAdjacency(good),
            node_info = data.frame(name = c("G1", "G2", "M1")),
            celltype = "Mono"),
        "feature column")
    expect_error(
        CelltypeNetworkResult(networkPrecision(good), networkAdjacency(good),
            node_info = networkNodes(good), celltype = c("Mono", "NK")),
        "single character string")
    expect_error(
        CelltypeNetworkResult(networkPrecision(good), networkAdjacency(good),
            node_info = networkNodes(good)[1:2, ], celltype = "Mono"),
        "must equal nrow")
})

test_that("CelltypeNetworkResult validity ties node_info to the matrices", {
    good <- .toy_network()
    precision <- networkPrecision(good)
    adjacency <- networkAdjacency(good)
    nodes <- networkNodes(good)

    ## a node table in another order would relabel the G1-M1 edge
    expect_error(
        CelltypeNetworkResult(precision, adjacency,
            node_info = nodes[c(2, 3, 1), ], celltype = "Mono"),
        "names of precision_matrix must match node_info\\$feature")
    expect_error(
        CelltypeNetworkResult(precision, unname(adjacency),
            node_info = nodes, celltype = "Mono",
            stability = adjacency[c(2, 3, 1), c(2, 3, 1)] * 0.9),
        "names of stability_scores must match")
    adjacency_named_wrongly <- adjacency
    dimnames(adjacency_named_wrongly) <- list(c("A", "B", "C"), NULL)
    expect_error(
        CelltypeNetworkResult(unname(precision), adjacency_named_wrongly,
            node_info = nodes, celltype = "Mono"),
        "names of adjacency_matrix must match")
    ## unnamed matrices are matched by position
    expect_true(methods::validObject(CelltypeNetworkResult(
        unname(precision), unname(adjacency), node_info = nodes,
        celltype = "Mono")))

    expect_error(
        CelltypeNetworkResult(precision, adjacency,
            node_info = S4Vectors::DataFrame(), celltype = "Mono"),
        "nrow\\(node_info\\) \\[0\\] must equal")
    expect_error(
        CelltypeNetworkResult(precision, adjacency,
            node_info = data.frame(feature = c("G1", NA, "M1")),
            celltype = "Mono"),
        "must not contain NA")
    expect_error(
        CelltypeNetworkResult(matrix(numeric(), 0L, 0L), adjacency,
            node_info = nodes, celltype = "Mono"),
        "adjacency_matrix dimensions must match")
})

test_that("CelltypeNetworkResult validity needs a symmetric 0/1 adjacency", {
    good <- .toy_network()
    precision <- networkPrecision(good)
    adjacency <- networkAdjacency(good)
    nodes <- networkNodes(good)

    weighted <- adjacency
    weighted["G1", "M1"] <- weighted["M1", "G1"] <- 5
    expect_error(
        CelltypeNetworkResult(precision, weighted, node_info = nodes,
            celltype = "Mono"),
        "only 0 and 1")
    missing_edge <- adjacency
    missing_edge["G1", "G2"] <- NA
    expect_error(
        CelltypeNetworkResult(precision, missing_edge, node_info = nodes,
            celltype = "Mono"),
        "only 0 and 1")
    one_way <- adjacency
    one_way["G1", "G2"] <- 1L
    expect_error(
        CelltypeNetworkResult(precision, one_way, node_info = nodes,
            celltype = "Mono"),
        "must be symmetric")
    expect_true(methods::validObject(CelltypeNetworkResult(precision,
        adjacency != 0, node_info = nodes, celltype = "Mono")))
})

test_that("comparison accessors agree with the summary", {
    comp <- run_comparative_network(.toy_network("Mono", -0.4),
                                    .toy_network("NK", 0))
    expect_true(methods::validObject(comp))
    summary <- comparisonSummary(comp)
    expect_identical(summary$celltypes, c("Mono", "NK"))
    expect_identical(nrow(sharedEdges(comp)), summary$n_shared)
    expect_identical(nrow(diffEdges(comp)), summary$n_differential)
    specific <- celltypeSpecificEdges(comp)
    expect_identical(names(specific), "Mono")
    expect_identical(as.character(specific$Mono$node_i), "G1")
    expect_identical(as.character(specific$Mono$node_j), "M1")
    expect_identical(as.character(diffEdges(comp)$category), "unique_Mono")
})

test_that("NetworkComparison() and its validity", {
    shared <- S4Vectors::DataFrame(node_i = "G1", node_j = "M1",
        category = "shared", max_diff = 0.1)
    differential <- S4Vectors::DataFrame(node_i = "G2", node_j = "M1",
        category = "differential", max_diff = 0.5)
    comp <- NetworkComparison(diff_edges = differential,
        shared_edges = shared,
        celltype_specific = list(Mono = differential[0, ]),
        summary = list(celltypes = c("Mono", "NK")))
    expect_true(methods::validObject(comp))
    expect_identical(sharedEdges(comp), shared)
    expect_identical(diffEdges(comp), differential)
    expect_identical(names(celltypeSpecificEdges(comp)), "Mono")
    expect_identical(comparisonSummary(comp)$celltypes, c("Mono", "NK"))

    empty <- NetworkComparison()
    expect_true(methods::validObject(empty))
    expect_identical(nrow(diffEdges(empty)), 0L)

    expect_error(NetworkComparison(
        shared_edges = S4Vectors::DataFrame(from = "G1", to = "M1")),
        "node_i and node_j")
    expect_error(NetworkComparison(shared_edges = differential),
        "category 'shared'")
    expect_error(NetworkComparison(diff_edges = shared),
        "must not contain shared edges")
    expect_error(NetworkComparison(celltype_specific = list(differential)),
        "named by cell type")
    expect_error(NetworkComparison(
        celltype_specific = list(B = differential),
        summary = list(celltypes = c("Mono", "NK"))),
        "among summary")
    expect_error(NetworkComparison(
        celltype_specific = list(Mono = as.data.frame(differential))),
        "list of DataFrames")
})
