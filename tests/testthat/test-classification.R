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
    expect_equal(
        names(families),
        c("Gene", "PlantTFDB_family", "PlantTFDB_subfamily", "TAPscan_family",
          "TAPscan_subfamily", "TAP_class", "PlantTFClass_superclass",
          "PlantTFClass_class", "PlantTFClass_family", "CISBP_family",
          "TFClass_superclass", "TFClass_class", "TFClass_family")
    )
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

    ptfdb <- classify_tfs(gsu_annotation, scheme = "PlantTFDB")
    expect_equal(
        names(ptfdb), c("Gene", "PlantTFDB_family", "PlantTFDB_subfamily")
    )
    expect_equal(nrow(ptfdb), sum(!is.na(all_fams$PlantTFDB_family)))

    tap <- classify_tfs(gsu_annotation, scheme = "TAPscan")
    expect_equal(
        names(tap),
        c("Gene", "TAPscan_family", "TAPscan_subfamily", "TAP_class")
    )
    expect_false(any(is.na(tap$TAPscan_family)))
    expect_true(all(tap$TAP_class %in% c("TF", "TR", "PT", NA)))
    expect_true(all(c("TF", "TR") %in% tap$TAP_class))

    tfclass <- classify_tfs(gsu_annotation, scheme = "PlantTFClass")
    expect_equal(
        names(tfclass),
        c("Gene", "PlantTFClass_superclass", "PlantTFClass_class",
          "PlantTFClass_family")
    )
    expect_false(any(is.na(tfclass$PlantTFClass_superclass)))

    cisbp <- classify_tfs(gsu_annotation, scheme = "CISBP")
    expect_equal(names(cisbp), c("Gene", "CISBP_family"))
    expect_false(any(is.na(cisbp$CISBP_family)))

    tfc <- classify_tfs(gsu_annotation, scheme = "TFClass")
    expect_equal(
        names(tfc),
        c("Gene", "TFClass_superclass", "TFClass_class", "TFClass_family")
    )
    expect_false(any(is.na(tfc$TFClass_class)))
    expect_equal(nrow(tfc), sum(!is.na(all_fams$TFClass_superclass)))

    expect_error(classify_tfs(gsu_annotation, scheme = "fake"))
})


test_that("classify_tfs() uses multiple schemes and available annotations", {
    all_fams <- classify_tfs(gsu_annotation)

    two <- classify_tfs(gsu_annotation, scheme = c("CISBP", "PlantTFDB"))
    expect_equal(
        names(two),
        c("Gene", "PlantTFDB_family", "PlantTFDB_subfamily", "CISBP_family")
    )
    classified <- !is.na(all_fams$PlantTFDB_family) |
        !is.na(all_fams$CISBP_family)
    expect_equal(two$Gene, all_fams$Gene[classified])

    # 'all' uses only schemes with available domain annotation
    cisbp_only <- classify_tfs(gsu_annotation["CISBP"])
    expect_equal(
        cisbp_only,
        classify_tfs(gsu_annotation, scheme = c("CISBP", "TFClass"))
    )
    plant <- classify_tfs(gsu_annotation[c("PlantTFDB", "TAPscan")])
    expect_equal(
        names(plant),
        setdiff(names(all_fams), c("CISBP_family", .scheme_columns("TFClass")))
    )

    # Schemes without domain annotation
    msg <- "Missing domain annotation"
    expect_error(
        classify_tfs(gsu_annotation["PlantTFDB"], scheme = "CISBP"), msg
    )
    expect_error(
        classify_tfs(gsu_annotation["PlantTFDB"], scheme = "PlantTFClass"),
        msg
    )
    expect_error(
        classify_tfs(gsu_annotation["PlantTFDB"], scheme = "TFClass"), msg
    )
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
})


test_that(".list_domains() returns a named character vector", {
    dom <- .list_domains()
    expect_true(is.character(dom))
    expect_equal(length(dom), 63)
    expect_false(any(duplicated(dom)))
})


test_that("all domains used for classification have profile HMMs", {
    lines <- readLines(.scheme_hmm_path("PlantTFDB"))
    hmm_names <- gsub("^NAME +", "", lines[startsWith(lines, "NAME")])
    expect_true(all(.list_domains() %in% hmm_names))
})


test_that(".classify_planttfdb() classifies AP2/ERF and B3 families", {
    expect_equal(classify_one("PF00847"), "ERF")
    expect_equal(classify_one(rep("PF00847", 2)), "AP2")
    expect_equal(classify_one(c("PF02362", "PF00847")), "RAV")
    expect_equal(classify_one("PF02362"), "B3")
    expect_equal(classify_one(c("PF06507", "PF02362")), "ARF")
    expect_null(classify_one("fakedomain"))
})


test_that(".classify_planttfdb() classifies C2C2 and GARP families", {
    expect_equal(classify_one("PF00320"), "GATA")
    expect_equal(classify_one(c("PF00643", "PF06203")), "CO-like")
    expect_equal(classify_one(c("PF00643", "PF00643", "PF06203")), "CO-like")
    expect_equal(classify_one("PF02701"), "Dof")
    expect_equal(classify_one("PF06943"), "LSD")
    expect_null(classify_one(c("PF06943", "PF00656")))
    expect_equal(classify_one("PF04690"), "YABBY")
    expect_equal(classify_one("G2-like"), "G2-like")
    expect_equal(classify_one(c("PF00072", "G2-like")), "ARR-B")
})


test_that(".classify_planttfdb() classifies HB families", {
    expect_equal(classify_one("PF00046"), "HB-other")
    expect_equal(classify_one(c("PF00046", "HD-ZIP_I/II")), "HD-ZIP")
    expect_equal(classify_one(c("PF01852", "PF00046")), "HD-ZIP")
    expect_equal(classify_one(c("PF03789", "PF00046")), "TALE")
    expect_equal(classify_one(c("BELL", "PF00046")), "TALE")
    expect_equal(classify_one(c("Wus_type_Homeobox", "PF00046")), "WOX")
    expect_equal(classify_one(c("PF00628", "PF00046")), "HB-PHD")
})


test_that(".classify_planttfdb() classifies MADS, MYB and NF-Y families", {
    expect_equal(classify_one("PF00319"), "M-type")
    expect_equal(classify_one(c("PF01486", "PF00319")), "MIKC")
    expect_equal(classify_one("PF00249"), "MYB-related")
    expect_equal(classify_one(rep("PF00249", 3)), "MYB")
    expect_null(classify_one(c("PF00249", "PF04433")))
    expect_equal(classify_one("PF02045"), "NF-YA")
    expect_equal(classify_one("NF-YB"), "NF-YB")
    expect_equal(classify_one("NF-YC"), "NF-YC")
})


test_that(".classify_planttfdb() classifies smaller families", {
    small <- c(
        PF06217 = "BBR-BPC", PF05687 = "BES1", PF00010 = "bHLH",
        PF00170 = "bZIP", PF00096 = "C2H2", PF00642 = "C3H",
        PF03859 = "CAMTA", PF03638 = "CPP", PF02319 = "E2F/DP",
        PF04873 = "EIL", PF03101 = "FAR1", PF04504 = "GeBP",
        PF03514 = "GRAS", `HRT-like` = "HRT-like", PF00447 = "HSF",
        PF03195 = "LBD", PF01698 = "LFY", PF02365 = "NAC",
        PF01422 = "NF-X1", PF02042 = "Nin-like", PF08744 = "NZZ/SPL",
        PF04689 = "S1Fa-like", SAP = "SAP", PF03110 = "SBP",
        PF05142 = "SRS", STAT = "STAT", PF03634 = "TCP",
        trihelix = "Trihelix", VOZ = "VOZ", PF08536 = "Whirly",
        PF03106 = "WRKY", PF04770 = "ZF-HD"
    )
    fams <- vapply(names(small), classify_one, character(1))
    expect_equal(fams, small)

    expect_equal(classify_one(c("PF08879", "PF08880")), "GRF")
    expect_equal(classify_one(rep("PF00643", 2)), "DBB")
    expect_null(classify_one(c("PF00096", "PF00929")))
    expect_null(classify_one(c("PF00642", "PF00076")))
})


test_that(".classify_planttfdb() handles multiple genes and families", {
    annot <- data.frame(
        Gene = c("g1", "g1", "g2", "g3", "g3"),
        Domain = c("PF00847", "PF00847", "PF00010", "G2-like", "PF00249")
    )
    fams <- .classify_planttfdb(annot)
    expect_equal(fams$Gene, c("g1", "g2", "g3"))
    expect_equal(fams$Subfamily, c("AP2", "bHLH", "G2-like;MYB-related"))
    expect_equal(fams$Family, c("AP2/ERF", "bHLH", "GARP;MYB superfamily"))

    empty <- .classify_planttfdb(data.frame(Gene = "g1", Domain = "fake"))
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
