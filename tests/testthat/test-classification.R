#----Load data------------------------------------------------------------------
data(gsu_annotation)

# Helper function to classify a single gene with a set of domains
classify_one <- function(domains) {
    fam <- .classify_planttfdb(data.frame(Gene = "g1", Domain = domains))
    if(nrow(fam) == 0) { return(NULL) }
    return(fam$Subfamily)
}

#----Start tests----------------------------------------------------------------
test_that("classify_tfs() returns a data frame with families", {
    families <- classify_tfs(gsu_annotation)

    expect_true(is(families, "data.frame"))
    cols <- unlist(lapply(list_schemes()$Scheme, .scheme_columns))
    expect_equal(names(families), c("Gene", cols))
    expect_true(nrow(families) > 0)
    expect_false(any(duplicated(families$Gene)))

    # More specific levels are never missing for classified genes
    expect_equal(
        is.na(families$PlantTFDB_family), is.na(families$PlantTFDB_subfamily)
    )
    expect_equal(
        is.na(families$TAPscan_family), is.na(families$TAPscan_subfamily)
    )
    expect_equal(
        is.na(families$PlantTFClass_class), is.na(families$PlantTFClass_family)
    )
})


test_that("classify_tfs() returns classifications for a single scheme", {
    all_fams <- classify_tfs(gsu_annotation)

    lapply(list_schemes()$Scheme, function(scheme) {
        fams <- classify_tfs(gsu_annotation, scheme = scheme)
        cols <- .scheme_columns(scheme)
        expect_true(is(fams, "data.frame"))
        expect_equal(names(fams), c("Gene", cols))
        expect_true(nrow(fams) > 0)

        # Genes classified by the scheme, as when using all schemes
        expect_false(any(is.na(fams[[cols[1]]])))
        expect_true(setequal(
            fams$Gene, all_fams$Gene[!is.na(all_fams[[cols[1]]])]
        ))
    })

    tap <- classify_tfs(gsu_annotation, scheme = "TAPscan")
    expect_true(all(tap$TAP_class %in% c("TF", "TR", "PT", NA)))

    expect_error(classify_tfs(gsu_annotation, scheme = "fake"))
})


test_that("classify_tfs() uses multiple schemes and available annotations", {
    all_fams <- classify_tfs(gsu_annotation)

    two <- classify_tfs(gsu_annotation, scheme = c("CISBP", "PlantTFDB"))
    expect_equal(
        names(two),
        c("Gene", .scheme_columns("PlantTFDB"), .scheme_columns("CISBP"))
    )
    classified <- !is.na(all_fams$PlantTFDB_family) |
        !is.na(all_fams$CISBP_family)
    expect_true(setequal(two$Gene, all_fams$Gene[classified]))

    # 'all' uses only schemes with available domain annotation
    lapply(list("CISBP", c("CISBP", "TFClass"), c("PlantTFDB", "TAPscan")),
           function(sets) {
        fams <- classify_tfs(gsu_annotation[sets])
        schemes <- .check_scheme("all", sets)
        expect_equal(fams, classify_tfs(gsu_annotation, scheme = schemes))
    })

    # Schemes without domain annotation
    msg <- "Missing domain annotation"
    lapply(c("CISBP", "PlantTFClass", "TFClass"), function(scheme) {
        expect_error(
            classify_tfs(gsu_annotation["PlantTFDB"], scheme = scheme), msg
        )
    })
})


test_that("classify_tfs() checks input", {
    msg <- "as returned by annotate_domains"
    expect_error(classify_tfs(gsu_annotation$PlantTFDB), msg)
    expect_error(classify_tfs(unname(gsu_annotation)), msg)
    expect_error(classify_tfs(list(fake = gsu_annotation$CISBP)), msg)
    expect_error(classify_tfs(list(PlantTFDB = 1, TAPscan = 2)), msg)

    bad <- gsu_annotation
    bad$TAPscan <- bad$TAPscan[, c("Gene", "Domain")]
    expect_error(classify_tfs(bad), msg)

    expect_error(classify_tfs(gsu_annotation, min_shared = -1), "min_shared")
    expect_error(classify_tfs(gsu_annotation, min_shared = NA), "min_shared")
    expect_error(
        classify_tfs(gsu_annotation, min_shared = c(0.1, 0.2)), "min_shared"
    )
    expect_error(classify_tfs(gsu_annotation, min_shared = "a"), "min_shared")
})


test_that(".list_domains() returns a named character vector", {
    dom <- .list_domains()
    expect_true(is.character(dom))
    expect_true(length(dom) > 0)
    expect_false(is.null(names(dom)))
    expect_false(any(duplicated(dom)))
})


test_that("all domains used for classification have profile HMMs", {
    lines <- readLines(.scheme_hmm_path("PlantTFDB"))
    hmm_names <- gsub("^NAME +", "", lines[startsWith(lines, "NAME")])
    expect_true(all(.list_domains() %in% hmm_names))
})


test_that(".classify_planttfdb() assigns domains to valid subfamilies", {
    data(planttfdb_scheme)

    # Single domains: valid subfamilies (or not classified)
    subfams <- unlist(lapply(.list_domains(), classify_one))
    expect_true(length(subfams) > 0)
    expect_true(all(
        unlist(strsplit(subfams, ";")) %in% planttfdb_scheme$Subfamily
    ))

    # Order and copies of domains do not change the classification
    doms <- c("PF02362", "PF00847")
    expect_equal(classify_one(doms), classify_one(rev(doms)))
    expect_equal(
        classify_one(c("PF00643", "PF06203")),
        classify_one(c("PF00643", "PF06203", "PF06203"))
    )

    # Number of copies of a domain can define different subfamilies
    expect_false(
        identical(classify_one("PF00249"), classify_one(rep("PF00249", 3)))
    )

    # Additional domains can define different subfamilies or exclude genes
    expect_false(identical(
        classify_one("PF00046"), classify_one(c("PF03789", "PF00046"))
    ))
    expect_null(classify_one(c("PF00096", "PF00929")))

    # Domains that are not used for classification
    expect_null(classify_one("fakedomain"))
})


test_that(".classify_planttfdb() handles multiple genes and families", {
    annot <- data.frame(
        Gene = c("g1", "g1", "g2", "g3", "g3"),
        Domain = c("PF00847", "PF00847", "PF00010", "G2-like", "PF00249")
    )
    fams <- .classify_planttfdb(annot)
    expect_true(is(fams, "data.frame"))
    expect_equal(names(fams), c("Gene", "Family", "Subfamily"))
    expect_true(setequal(fams$Gene, c("g1", "g2", "g3")))

    # Genes with multiple families have families and subfamilies separated
    # by ';' (in the same order)
    g3 <- fams[fams$Gene == "g3", ]
    expect_true(grepl(";", g3$Subfamily))
    expect_equal(
        lengths(strsplit(g3$Family, ";")), lengths(strsplit(g3$Subfamily, ";"))
    )

    empty <- .classify_planttfdb(data.frame(Gene = "g1", Domain = "fake"))
    expect_true(is(empty, "data.frame"))
    expect_equal(nrow(empty), 0)
})


test_that("all PlantTFDB subfamilies have a family in the scheme", {
    data(planttfdb_scheme)
    m <- matrix(
        0L, 1, length(.list_domains()),
        dimnames = list("g1", unname(.list_domains()))
    )
    labels <- unlist(lapply(.planttfdb_rules(m), names))
    expect_true(all(labels %in% planttfdb_scheme$Subfamily))
    expect_equal(
        planttfdb_families,
        unique(planttfdb_scheme[, c("Family", "Subfamily")]),
        ignore_attr = TRUE
    )
})
