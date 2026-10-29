#----Load data------------------------------------------------------------------
# Helper function to create input for .classify_planttfclass()
make_families <- function(ptfdb, tap_fam, tap_sub = tap_fam) {
    data.frame(
        Gene = sprintf("g%d", seq_along(tap_fam)),
        PlantTFDB_subfamily = ptfdb,
        TAPscan_family = tap_fam,
        TAPscan_subfamily = tap_sub
    )
}

# Helper function to check that all classifications are valid
expect_valid_planttfclass <- function(res) {
    cols <- .scheme_columns("PlantTFClass")
    expect_true(all(cols %in% names(res)))
    valid <- list(
        planttfclass_map$superclass, planttfclass_map$class,
        planttfclass_map$family
    )
    lapply(seq_along(cols), function(i) {
        values <- unlist(strsplit(res[[cols[i]]], ";"))
        expect_true(all(values[!is.na(values)] %in% valid[[i]]))
    })
}

#----Start tests----------------------------------------------------------------
test_that("correspondence table has one superclass and class per family", {
    map <- unique(planttfclass_map[, c("family", "class", "superclass")])
    expect_false(any(duplicated(map$family)))
    expect_false(any(is.na(planttfclass_map$family)))
    expect_false(any(is.na(planttfclass_map$class)))
    expect_false(any(is.na(planttfclass_map$superclass)))

    data(planttfclass_scheme)
    expect_true(is(planttfclass_scheme, "data.frame"))
    expect_equal(nrow(planttfclass_scheme), nrow(planttfclass_map))
})


test_that(".classify_planttfclass() adds Plant-TFClass classifications", {
    fams <- make_families(
        ptfdb = c("WRKY", NA, NA, NA, NA, NA),
        tap_fam = c("WRKY", "bHLH", "bHLH", "HDZ", "LFY", "FHA"),
        tap_sub = c("WRKY", "bHLH", "bHLH_TCP", "C3HDZ", "LFY", "FHA")
    )
    res <- .classify_planttfclass(fams)

    expect_true(is(res, "data.frame"))
    expect_equal(nrow(res), nrow(fams))
    expect_equal(names(res), c(names(fams), .scheme_columns("PlantTFClass")))
    expect_valid_planttfclass(res)

    # TFs of TAPscan families in the table are classified at all levels
    expect_false(any(is.na(res$PlantTFClass_family[1:5])))
    expect_equal(
        is.na(res$PlantTFClass_superclass), is.na(res$PlantTFClass_family)
    )

    # TAPscan families that are not TFs (e.g., FHA) are not classified
    expect_true(is.na(res$PlantTFClass_family[6]))
})


test_that("subfamily-specific matches take precedence over family matches", {
    fams <- make_families(
        ptfdb = NA, tap_fam = c("C2H2", "C2H2"), tap_sub = c("C2H2", "C2H2_IDD")
    )
    res <- .classify_planttfclass(fams)
    expect_valid_planttfclass(res)
    expect_false(any(is.na(res$PlantTFClass_family)))
    expect_false(res$PlantTFClass_family[1] == res$PlantTFClass_family[2])
})


test_that("ambiguous classifications are resolved with PlantTFDB families", {
    fams <- make_families(
        ptfdb = c("ERF", "AP2", NA, "bZIP", "G2-like;ERF"),
        tap_fam = "AP2"
    )
    res <- .classify_planttfclass(fams)
    expect_valid_planttfclass(res)

    # Resolved with PlantTFDB subfamilies: a single family
    expect_false(any(grepl(";", res$PlantTFClass_family[c(1, 2, 5)])))
    expect_false(res$PlantTFClass_family[1] == res$PlantTFClass_family[2])
    expect_equal(res$PlantTFClass_family[1], res$PlantTFClass_family[5])

    # Not resolved (no or uninformative PlantTFDB subfamily): all families
    expect_true(all(grepl(";", res$PlantTFClass_family[3:4])))
    expect_equal(res$PlantTFClass_family[3], res$PlantTFClass_family[4])
    families <- strsplit(res$PlantTFClass_family[3], ";")[[1]]
    expect_true(all(res$PlantTFClass_family[1:2] %in% families))
})


test_that(".resolve_planttfclass() handles empty input", {
    res <- .resolve_planttfclass(planttfclass_map[0, ], NA)
    expect_equal(res, rep(NA_character_, 3))

    fams <- make_families(character(0), character(0), character(0))
    res <- .classify_planttfclass(fams)
    expect_true(is(res, "data.frame"))
    expect_equal(nrow(res), 0)
})
