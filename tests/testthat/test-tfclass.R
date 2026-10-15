#----Load data------------------------------------------------------------------
data(tfclass_scheme)

# Helper function to get the TFClass classification of a CIS-BP family
tfclass_of <- function(cisbp_family) {
    res <- .classify_tfclass(
        data.frame(Gene = "g1", CISBP_family = cisbp_family)
    )
    return(unlist(res[1, c("TFClass_superclass", "TFClass_class",
                           "TFClass_family")], use.names = FALSE))
}

#----Start tests----------------------------------------------------------------
test_that(".classify_tfclass() uses the correspondence table", {
    expect_equal(
        tfclass_of("Forkhead"),
        c("Helix-turn-helix domains", "Fork head / winged helix factors",
          "FOX")
    )

    # CIS-BP families that correspond to multiple TFClass families
    expect_equal(
        tfclass_of("Homeodomain"),
        c("Helix-turn-helix domains", "Homeo domain factors", NA)
    )

    # Manually curated families with multiple DBDs
    expect_equal(tfclass_of("Homeodomain,POU")[3], "POU")
    expect_equal(tfclass_of("C2H2 ZF,Homeodomain")[3], "HD-ZF")
    expect_equal(tfclass_of("Homeodomain,Paired box")[2], "Paired box factors")
})


test_that(".classify_tfclass() resolves other families with multiple DBDs", {
    # DBDs in the same TFClass class and family
    expect_equal(tfclass_of("BED ZF,C2H2 ZF")[2], "C2H2 zinc finger factors")
    expect_true(is.na(tfclass_of("BED ZF,C2H2 ZF")[3]))
    expect_equal(tfclass_of("CxxC,MBD")[3], "CpG-binding proteins")

    # DBDs in different TFClass classes (ordered by TFClass ID)
    expect_equal(
        tfclass_of("bHLH,T-box"),
        c("Basic domains;Immunoglobulin fold",
          "Basic helix-loop-helix factors (bHLH);T-Box factors", NA)
    )

    # Families not classified by TFClass
    expect_equal(tfclass_of("Zinc cluster"), rep(NA_character_, 3))
    expect_equal(tfclass_of("AP2,B3"), rep(NA_character_, 3))
    expect_equal(tfclass_of(NA), rep(NA_character_, 3))
})


test_that(".classify_tfclass() handles multiple genes", {
    fams <- data.frame(
        Gene = c("g1", "g2", "g3", "g4"),
        CISBP_family = c("Forkhead", "WRKY", "Forkhead", NA)
    )
    res <- .classify_tfclass(fams)
    expect_equal(names(res), c(names(fams), .scheme_columns("TFClass")))
    expect_equal(res$TFClass_family, c("FOX", NA, "FOX", NA))
})


test_that("TFClass correspondence table is consistent", {
    expect_equal(tfclass_map$cisbp_family, tfclass_scheme$CISBP_family)
    expect_false(any(duplicated(tfclass_scheme$CISBP_family)))
    expect_false(any(is.na(tfclass_scheme$TFClass_class)))

    # All single-DBD CIS-BP families in the table are in cisbp_scheme
    single <- tfclass_scheme$CISBP_family[
        !grepl(",", tfclass_scheme$CISBP_family)
    ]
    expect_true(all(single %in% cisbp_scheme$Family))
})
