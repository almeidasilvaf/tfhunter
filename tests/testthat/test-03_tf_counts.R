#----Load data------------------------------------------------------------------
data(gsu_annotation)
gsu_families <- classify_tfs(gsu_annotation)

# Simulate TF classifications for 4 species by sampling 100 genes
set.seed(123)
families <- list(
    Gsu1 = gsu_families[sample(nrow(gsu_families), 100), ],
    Gsu2 = gsu_families[sample(nrow(gsu_families), 100), ],
    Gsu3 = gsu_families[sample(nrow(gsu_families), 100), ],
    Gsu4 = gsu_families[sample(nrow(gsu_families), 100), ]
)

# Create species metadata (in an order different from `families`)
species_metadata <- data.frame(
    row.names = rev(names(families)),
    Division = "Rhodophyta",
    Origin = c("Brazil", "China", "Belgium", "US")
)

# Small example with known counts
toy <- data.frame(
    Gene = c("g1", "g2", "g3", "g4", "g5"),
    PlantTFDB_family = c("GARP;MYB superfamily", "bHLH", NA, NA, "C2C2"),
    PlantTFDB_subfamily = c("G2-like;MYB-related", "bHLH", NA, NA, "Dof"),
    TAPscan_family = c("MYB", "bHLH", "RWP-RK", "FHA", NA),
    TAPscan_subfamily = c("MYB-related", "bHLH", "RKD", "FHA", NA),
    TAP_class = c("TF", "TF", "TF", "TR", NA),
    PlantTFClass_superclass = c(
        "Helix-turn-helix domains", "Basic domains",
        "Helix-turn-helix domains", NA, NA
    ),
    PlantTFClass_class = c(
        "Tryptophan cluster factors", "Basic helix-loop-helix factors (bHLH)",
        "RWP-RK", NA, NA
    ),
    PlantTFClass_family = c("MYB-related", "bHLH", "RWP-RK", NA, NA),
    CISBP_family = c("Myb/SANT", "bHLH", "RWP-RK", NA, "Dof,Myb/SANT"),
    TFClass_superclass = c(
        "Helix-turn-helix domains", "Basic domains", NA, NA,
        "Helix-turn-helix domains"
    ),
    TFClass_class = c(
        "Tryptophan cluster factors", "Basic helix-loop-helix factors (bHLH)",
        NA, NA, "Tryptophan cluster factors"
    ),
    TFClass_family = c("Myb-SANT", NA, NA, NA, "Myb-SANT")
)
# Helper function to get unique values of a variable (split by ';')
level_values <- function(x) {
    x <- unlist(strsplit(x[!is.na(x)], ";", fixed = TRUE))
    return(unique(x))
}

# Helper function to check counts for a scheme and level in the toy example
check_counts <- function(scheme, level) {
    se <- get_tf_counts(list(sp1 = toy), scheme = scheme, level = level)
    counts <- SummarizedExperiment::assay(se)
    col <- .scheme_levels(scheme)[[level]]

    expect_true(is(se, "SummarizedExperiment"))
    expect_true(is.integer(counts))
    expect_true(all(counts >= 0))
    expect_true(setequal(rownames(se), level_values(toy[[col]])))

    # Each gene is counted once for each of its families
    n_values <- lengths(strsplit(toy[[col]][!is.na(toy[[col]])], ";"))
    expect_equal(sum(counts), sum(n_values))
}

#----Start tests----------------------------------------------------------------
test_that("get_tf_counts() returns a SummarizedExperiment object", {
    se <- get_tf_counts(families, species_metadata)
    expect_true(is(se, "SummarizedExperiment"))
    expect_equal(colnames(se), names(families))
    expect_equal(se$Origin, species_metadata[colnames(se), "Origin"])
    expect_true(is.integer(SummarizedExperiment::assay(se)))
    expect_true(nrow(se) > 0)

    # CIS-BP is the default scheme
    expect_equal(
        rownames(se), rownames(get_tf_counts(families, scheme = "CISBP"))
    )

    # Without species metadata
    se2 <- get_tf_counts(families)
    expect_equal(ncol(SummarizedExperiment::colData(se2)), 0)
    expect_equal(
        SummarizedExperiment::assay(se2), SummarizedExperiment::assay(se)
    )
})


test_that("get_tf_counts() counts TFs at each level of each scheme", {
    lapply(list_schemes()$Scheme, function(scheme) {
        lapply(names(.scheme_levels(scheme)), function(level) {
            check_counts(scheme, level)
        })
    })
})


test_that("get_tf_counts() adds higher levels to rowData", {
    lapply(list_schemes()$Scheme, function(scheme) {
        levels <- .scheme_levels(scheme)
        lapply(seq_along(levels), function(i) {
            se <- get_tf_counts(
                list(sp1 = toy), scheme = scheme, level = names(levels)[i]
            )
            rd <- SummarizedExperiment::rowData(se)
            expect_true(all(levels[seq_len(i - 1)] %in% names(rd)))
            expect_false(any(levels[seq(i, length(levels))] %in% names(rd)))
        })
    })

    # Values of higher levels are those of genes in each feature
    se <- get_tf_counts(
        list(sp1 = toy), scheme = "TAPscan", level = "subfamily"
    )
    rd <- as.data.frame(SummarizedExperiment::rowData(se))
    expect_true(all(c("TAPscan_family", "TAP_class") %in% names(rd)))
    lapply(rownames(rd), function(f) {
        genes <- toy$TAPscan_subfamily %in% f
        expect_true(rd[f, "TAPscan_family"] %in% toy$TAPscan_family[genes])
    })
})


test_that("filtered tables and empty tables can be counted", {
    tfs <- lapply(families, function(x) x[x$TAP_class %in% "TF", ])
    se <- get_tf_counts(tfs, scheme = "TAPscan")
    expect_true(all(SummarizedExperiment::rowData(se)$TAP_class == "TF"))

    empty <- list(sp1 = toy[0, ], sp2 = toy)
    se <- get_tf_counts(empty, scheme = "TAPscan", level = "subfamily")
    counts <- colSums(SummarizedExperiment::assay(se))
    expect_true(counts[["sp1"]] == 0)
    expect_true(counts[["sp2"]] > 0)

    se <- get_tf_counts(list(sp1 = toy[0, ]))
    expect_true(is(se, "SummarizedExperiment"))
    expect_equal(nrow(se), 0)
})


test_that("get_tf_counts() checks input", {
    # Invalid families
    expect_error(get_tf_counts(data.frame()), "named list")
    expect_error(get_tf_counts(unname(families)), "named list")
    expect_error(get_tf_counts(list(sp1 = "a")), "named list")
    expect_error(
        get_tf_counts(list(sp1 = toy[, c("Gene", "TAPscan_family")])),
        "Missing in: sp1"
    )

    # Invalid levels
    expect_error(
        get_tf_counts(families, scheme = "TAPscan", level = "class"),
        "Valid levels are"
    )
    expect_error(get_tf_counts(families, level = c("family", "subfamily")))
    expect_error(get_tf_counts(families, scheme = "fake"))

    # Species metadata that does not match species
    expect_error(get_tf_counts(families, species_metadata[1:2, ]))
})


test_that(".check_level() validates levels for each scheme", {
    lapply(list_schemes()$Scheme, function(scheme) {
        levels <- names(.scheme_levels(scheme))
        lapply(levels, function(l) expect_true(.check_level(scheme, l)))
        expect_error(.check_level(scheme, "fake_level"), "Valid levels")
    })
})


test_that(".count_rowdata() combines values of higher levels", {
    long <- data.frame(
        Feature = c("f1", "f1", "f2", "f3"), Species = "sp1",
        Parent = c("p1", "p2", NA, "p1")
    )
    rd <- .count_rowdata(long, c("f1", "f2", "f3"))
    expect_true(is(rd, "data.frame"))
    expect_equal(rownames(rd), c("f1", "f2", "f3"))

    # Multiple values are separated by ';', and missing values are NA
    expect_true(setequal(strsplit(rd["f1", "Parent"], ";")[[1]], c("p1", "p2")))
    expect_true(is.na(rd["f2", "Parent"]))
})
