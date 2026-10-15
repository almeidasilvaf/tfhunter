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

# Helper function to get counts as a named vector for a single species
toy_counts <- function(...) {
    se <- get_tf_counts(list(sp1 = toy), ...)
    return(SummarizedExperiment::assay(se)[, "sp1"])
}

#----Start tests----------------------------------------------------------------
test_that("get_tf_counts() returns a SummarizedExperiment object", {
    se <- get_tf_counts(families, species_metadata)
    expect_true(is(se, "SummarizedExperiment"))
    expect_equal(colnames(se), names(families))
    expect_equal(se$Origin, c("US", "Belgium", "China", "Brazil"))
    expect_true(is.integer(SummarizedExperiment::assay(se)))
    expect_null(names(dimnames(SummarizedExperiment::assay(se))))

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


test_that("get_tf_counts() counts TFs at each level", {
    # PlantTFDB: genes with multiple families are counted for each family
    expect_equal(
        toy_counts(scheme = "PlantTFDB", level = "family"),
        c(bHLH = 1, C2C2 = 1, GARP = 1, `MYB superfamily` = 1)
    )
    expect_equal(
        toy_counts(scheme = "PlantTFDB", level = "subfamily"),
        c(bHLH = 1, Dof = 1, `G2-like` = 1, `MYB-related` = 1)
    )

    # TAPscan
    expect_equal(
        toy_counts(scheme = "TAPscan", level = "family"),
        c(bHLH = 1, FHA = 1, MYB = 1, `RWP-RK` = 1)
    )
    expect_equal(
        toy_counts(scheme = "TAPscan", level = "subfamily"),
        c(bHLH = 1, FHA = 1, `MYB-related` = 1, RKD = 1)
    )

    # Plant-TFClass
    expect_equal(
        toy_counts(scheme = "PlantTFClass", level = "family"),
        c(bHLH = 1, `MYB-related` = 1, `RWP-RK` = 1)
    )
    expect_equal(
        toy_counts(scheme = "PlantTFClass", level = "class"),
        c(`Basic helix-loop-helix factors (bHLH)` = 1, `RWP-RK` = 1,
          `Tryptophan cluster factors` = 1)
    )
    expect_equal(
        toy_counts(scheme = "PlantTFClass", level = "superclass"),
        c(`Basic domains` = 1, `Helix-turn-helix domains` = 2)
    )

    # CIS-BP: families with multiple DBDs are counted as a single family
    expect_equal(
        toy_counts(scheme = "CISBP"),
        c(bHLH = 1, `Dof,Myb/SANT` = 1, `Myb/SANT` = 1, `RWP-RK` = 1)
    )

    # TFClass: genes without a TFClass family are not counted per family
    expect_equal(
        toy_counts(scheme = "TFClass", level = "class"),
        c(`Basic helix-loop-helix factors (bHLH)` = 1,
          `Tryptophan cluster factors` = 2)
    )
    se <- get_tf_counts(list(sp1 = toy), scheme = "TFClass", level = "family")
    expect_equal(rownames(se), "Myb-SANT")
    expect_equal(as.numeric(SummarizedExperiment::assay(se)), 2)
})


test_that("get_tf_counts() adds higher levels to rowData", {
    rd <- function(...) {
        as.data.frame(SummarizedExperiment::rowData(
            get_tf_counts(list(sp1 = toy), ...)
        ))
    }

    expect_equal(ncol(rd(scheme = "PlantTFDB", level = "family")), 0)
    ptfdb <- rd(scheme = "PlantTFDB", level = "subfamily")
    expect_equal(ptfdb["Dof", "PlantTFDB_family"], "C2C2")
    expect_equal(ptfdb["MYB-related", "PlantTFDB_family"], "MYB superfamily")

    expect_equal(names(rd(scheme = "TAPscan")), "TAP_class")
    tap <- rd(scheme = "TAPscan", level = "subfamily")
    expect_equal(names(tap), c("TAPscan_family", "TAP_class"))
    expect_equal(tap["RKD", "TAPscan_family"], "RWP-RK")
    expect_equal(tap["FHA", "TAP_class"], "TR")

    tfclass <- rd(scheme = "PlantTFClass", level = "family")
    expect_equal(
        names(tfclass), c("PlantTFClass_superclass", "PlantTFClass_class")
    )
    expect_equal(tfclass["RWP-RK", "PlantTFClass_class"], "RWP-RK")
    expect_equal(
        names(rd(scheme = "PlantTFClass", level = "class")),
        "PlantTFClass_superclass"
    )
})


test_that("filtered tables and empty tables can be counted", {
    tfs <- lapply(families, function(x) x[x$TAP_class %in% "TF", ])
    se <- get_tf_counts(tfs, scheme = "TAPscan")
    expect_true(all(SummarizedExperiment::rowData(se)$TAP_class == "TF"))

    empty <- list(sp1 = toy[0, ], sp2 = toy)
    se <- get_tf_counts(empty, scheme = "TAPscan", level = "subfamily")
    expect_equal(unname(colSums(SummarizedExperiment::assay(se))), c(0, 4))

    se <- get_tf_counts(list(sp1 = toy[0, ]))
    expect_equal(dim(se), c(0, 1))
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
    expect_true(.check_level("PlantTFDB", "subfamily"))
    expect_true(.check_level("TAPscan", "subfamily"))
    expect_true(.check_level("PlantTFClass", "superclass"))
    expect_true(.check_level("CISBP", "family"))
    expect_true(.check_level("TFClass", "superclass"))
    expect_error(.check_level("TFClass", "subfamily"))
    expect_error(.check_level("PlantTFDB", "class"))
    expect_error(.check_level("CISBP", "subfamily"))
    expect_error(.check_level("TAPscan", "superclass"))
    expect_error(.check_level("PlantTFClass", "subfamily"))
})
