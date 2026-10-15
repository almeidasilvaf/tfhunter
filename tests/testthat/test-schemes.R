#----Start tests----------------------------------------------------------------
test_that("list_schemes() returns a data frame of schemes", {
    schemes <- list_schemes()
    expect_true(is(schemes, "data.frame"))
    expect_equal(
        schemes$Scheme,
        c("PlantTFDB", "TAPscan", "PlantTFClass", "CISBP", "TFClass")
    )

    # Levels match the variables in the output of classify_tfs()
    levels <- lapply(schemes$Scheme, function(s) names(.scheme_levels(s)))
    expect_equal(
        vapply(levels, paste, character(1), collapse = ", "), schemes$Levels
    )
})


test_that(".check_scheme() validates and expands schemes", {
    all_schemes <- list_schemes()$Scheme
    expect_equal(.check_scheme("all"), all_schemes)
    expect_equal(.check_scheme(c("CISBP", "TAPscan")), c("TAPscan", "CISBP"))
    expect_equal(
        .check_scheme("all", c("PlantTFDB", "TAPscan")),
        c("PlantTFDB", "TAPscan", "PlantTFClass")
    )
    expect_equal(.check_scheme("all", "CISBP"), c("CISBP", "TFClass"))
    expect_equal(.check_scheme("TFClass", "CISBP"), "TFClass")

    expect_error(.check_scheme("fake"), "must be 'all'")
    expect_error(.check_scheme(character(0)), "must be 'all'")
    expect_error(.check_scheme(1), "must be 'all'")
    expect_error(.check_scheme("CISBP", "TAPscan"), "Missing")
    expect_error(.check_scheme("TFClass", "TAPscan"), "Missing")
    expect_error(.check_scheme("all", character(0)), "No classification")
})


test_that(".annotation_sets() returns annotations required by schemes", {
    expect_equal(.annotation_sets("PlantTFClass"), c("PlantTFDB", "TAPscan"))
    expect_equal(.annotation_sets(c("CISBP", "TAPscan")), c("TAPscan", "CISBP"))
    expect_equal(.annotation_sets(c("CISBP", "TFClass")), "CISBP")
})


test_that(".scheme_hmm_path() decompresses HMMs once per session", {
    path <- .scheme_hmm_path("CISBP")
    expect_true(file.exists(path))
    expect_equal(.scheme_hmm_path("CISBP"), path)
    expect_true(startsWith(readLines(path, n = 1), "HMMER3"))
})
