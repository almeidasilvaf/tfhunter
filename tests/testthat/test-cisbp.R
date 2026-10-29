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
    # Single DBDs: families in the scheme, regardless of number of copies
    single <- vapply(cisbp_scheme$Pfam, classify_one, character(1))
    expect_true(all(single %in% cisbp_scheme$Family))
    expect_equal(classify_one(rep("PF00096", 5)), classify_one("PF00096"))

    # Multiple DBDs: family names sorted (case-insensitive) and separated
    # by ',', regardless of their order in the protein
    multi <- classify_one(c("PF00157", "PF00046"))
    expect_equal(multi, classify_one(c("PF00046", "PF00157")))
    parts <- strsplit(multi, ",")[[1]]
    expect_equal(parts, parts[order(tolower(parts))])
    expect_true(setequal(
        parts, c(classify_one("PF00157"), classify_one("PF00046"))
    ))

    # DBDs with the same family are counted once
    dup <- cisbp_scheme$Family[duplicated(cisbp_scheme$Family)]
    lapply(unique(dup), function(f) {
        pfam <- cisbp_scheme$Pfam[cisbp_scheme$Family == f]
        expect_equal(classify_one(pfam), f)
    })

    # Domains that are not DBDs in CIS-BP are ignored
    expect_equal(classify_one(c("PF00010", "fake")), classify_one("PF00010"))
    expect_null(classify_one("fake"))
})


test_that(".classify_cisbp() handles multiple genes", {
    annot <- data.frame(
        Gene = c("g2", "g1", "g2", "g3", "g1"),
        Domain = c("PF00249", "PF00010", "PF02701", "fake", "PF00010")
    )
    fams <- .classify_cisbp(annot)
    expect_true(is(fams, "data.frame"))
    expect_equal(names(fams), c("Gene", "Family"))
    expect_true(setequal(fams$Gene, c("g1", "g2")))
    expect_false(any(duplicated(fams$Gene)))
    expect_true(grepl(",", fams$Family[fams$Gene == "g2"]))

    empty <- .classify_cisbp(annot[annot$Gene == "g3", ])
    expect_true(is(empty, "data.frame"))
    expect_equal(names(empty), c("Gene", "Family"))
    expect_equal(nrow(empty), 0)
})


test_that("CIS-BP scheme and profile HMMs are consistent", {
    expect_true(is(cisbp_scheme, "data.frame"))
    expect_true(nrow(cisbp_scheme) > 0)
    expect_equal(cisbp_map$pfam, cisbp_scheme$Pfam)
    expect_equal(cisbp_map$family, cisbp_scheme$Family)
    expect_false(any(duplicated(cisbp_scheme$Pfam)))

    lines <- readLines(.scheme_hmm_path("CISBP"))
    hmm_names <- gsub("^NAME +", "", lines[startsWith(lines, "NAME")])
    expect_true(setequal(hmm_names, cisbp_scheme$Pfam))
})
