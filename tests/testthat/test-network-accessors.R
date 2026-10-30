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
