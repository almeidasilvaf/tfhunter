#----Load data------------------------------------------------------------------
data(cisbp_scheme)

# Helper function to classify a single gene with a set of domains
classify_one <- function(domains) {
    fam <- .classify_cisbp(data.frame(Gene = "g1", Domain = domains))
    if(nrow(fam) == 0) { return(NULL) }
    return(fam$Family)
}

#----Start tests----------------------------------------------------------------
test_that(".classify_cisbp() assigns families based on sets of DBDs", {
    # Single DBDs (and multiple copies of the same DBD)
    expect_equal(classify_one("PF00010"), "bHLH")
    expect_equal(classify_one(rep("PF00096", 5)), "C2H2 ZF")
    expect_equal(classify_one("PF00172"), "Zinc cluster")

    # Multiple DBDs: sorted alphabetically (case-insensitive), regardless
    # of their order in the protein
    expect_equal(classify_one(c("PF00157", "PF00046")), "Homeodomain,POU")
    expect_equal(classify_one(c("PF00046", "PF00157")), "Homeodomain,POU")
    expect_equal(classify_one(c("PF00096", "PF00010")), "bHLH,C2H2 ZF")
    expect_equal(classify_one(c("PF00046", "PF00096")), "C2H2 ZF,Homeodomain")

    # DBDs with the same family are counted once
    expect_equal(classify_one(c("PF01698", "PF17538")), "Leafy")

    # Domains that are not DBDs in CIS-BP are ignored
    expect_equal(classify_one(c("PF00010", "fake")), "bHLH")
    expect_null(classify_one("fake"))
})


test_that(".classify_cisbp() handles multiple genes", {
    annot <- data.frame(
        Gene = c("g2", "g1", "g2", "g3", "g1"),
        Domain = c("PF00249", "PF00010", "PF02701", "fake", "PF00010")
    )
    fams <- .classify_cisbp(annot)
    expect_equal(fams$Gene, c("g1", "g2"))
    expect_equal(fams$Family, c("bHLH", "Dof,Myb/SANT"))

    empty <- .classify_cisbp(annot[annot$Gene == "g3", ])
    expect_equal(names(empty), c("Gene", "Family"))
    expect_equal(nrow(empty), 0)
})


test_that("CIS-BP scheme and profile HMMs are consistent", {
    expect_equal(cisbp_map$pfam, cisbp_scheme$Pfam)
    expect_equal(cisbp_map$family, cisbp_scheme$Family)
    expect_false(any(duplicated(cisbp_scheme$Pfam)))

    lines <- readLines(.scheme_hmm_path("CISBP"))
    hmm_names <- gsub("^NAME +", "", lines[startsWith(lines, "NAME")])
    expect_true(setequal(hmm_names, cisbp_scheme$Pfam))
})
