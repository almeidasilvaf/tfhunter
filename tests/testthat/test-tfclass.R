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
    hd <- tfclass_of("Homeodomain")
    expect_equal(hd[1:2], c("Helix-turn-helix domains", "Homeo domain factors"))
    hd_families <- strsplit(hd[3], ";")[[1]]
    expect_true(all(c("HOX", "NK", "TALE-type HD") %in% hd_families))
    expect_equal(
        tfclass_of("SMAD"),
        c("beta-Hairpin exposed by an alpha/beta-scaffold",
          "SMAD/NF-1 DNA-binding domain factors", "SMAD;NF-1")
    )

    # AT hook: TFClass family is the class (not only HMGA in CIS-BP)
    expect_equal(tfclass_of("AT hook")[2:3], rep("A.T hook factors", 2))

    # Manually curated families with multiple DBDs
    expect_equal(tfclass_of("Homeodomain,POU")[3], "POU")
    expect_equal(tfclass_of("C2H2 ZF,Homeodomain")[3], "HD-ZF")
    expect_equal(tfclass_of("Homeodomain,Paired box")[2], "Paired box factors")
})


test_that(".classify_tfclass() resolves other families with multiple DBDs", {
    # DBDs in the same TFClass class: all possible families
    bed <- tfclass_of("BED ZF,C2H2 ZF")
    expect_equal(bed[2], "C2H2 zinc finger factors")
    expect_equal(
        bed[3],
        paste(tfclass_of("C2H2 ZF")[3], "BED zinc finger factors", sep = ";")
    )
    expect_equal(tfclass_of("CxxC,MBD")[3], "CpG-binding proteins")

    # DBDs in different TFClass classes (ordered by TFClass ID)
    expect_equal(
        tfclass_of("bHLH,T-box"),
        c("Basic domains;Immunoglobulin fold",
          "Basic helix-loop-helix factors (bHLH);T-Box factors",
          paste(tfclass_of("bHLH")[3], tfclass_of("T-box")[3], sep = ";"))
    )

    expect_equal(
        tfclass_of("AT hook,Homeodomain")[3],
        paste(tfclass_of("Homeodomain")[3], "A.T hook factors", sep = ";")
    )

    # All TFs with a TFClass class have a TFClass family
    expect_false(any(is.na(tfclass_scheme$TFClass_family)))

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
    expect_equal(
        lengths(strsplit(tfclass_scheme$TFClass_family, ";")),
        lengths(strsplit(tfclass_scheme$TFClass_ID, ";"))
    )

    # All single-DBD CIS-BP families in the table are in cisbp_scheme
    single <- tfclass_scheme$CISBP_family[
        !grepl(",", tfclass_scheme$CISBP_family)
    ]
    expect_true(all(single %in% cisbp_scheme$Family))
})


test_that(".tfclass_families() returns names of TFClass families", {
    fams <- .tfclass_families()
    expect_equal(names(fams), c("superclass", "class", "family", "id"))
    expect_false(any(duplicated(fams$id)))
    expect_true(all(.tfclass_hmm_ids() %in% fams$id))
    expect_true(all(paste0("2.3.", 1:4) %in% fams$id))
    expect_equal(fams$family[fams$id == "3.1.1"], "HOX")
    expect_false(any(grepl(";", unlist(fams))))
})


test_that(".kmer_table() encodes unique k-mers of each sequence", {
    seq <- Biostrings::AAStringSet(c(a = "ACDEFA", b = "ACDXACD", c = "AC"))
    km <- .kmer_table(seq, k = 3)

    # a: ACD, CDE, DEF, EFA; b: ACD (once; k-mers with X are ignored);
    # c: shorter than k. K-mers spanning two sequences are ignored
    expect_equal(km$seq, c(1, 1, 1, 1, 2))
    expect_equal(km$code[1], 0 * 400 + 1 * 20 + 2)
    expect_equal(km$code[5], km$code[1])
    expect_equal(nrow(.kmer_table(seq["c"], k = 3)), 0)
})


test_that(".shared_kmers() counts k-mers shared by queries and references", {
    ref <- Biostrings::AAStringSet(c(
        r1 = "MKVLAAGIRS", r2 = "MKVWWPPQQH", r3 = "MKVLAWWPPQ"
    ))
    query <- Biostrings::AAStringSet(c(q1 = "MKVLAAG", q2 = "TTTTT"))
    res <- .shared_kmers(query, ref, k = 3)

    # q1: MKV, KVL, VLA, LAA, AAG; MKV is in 3 references (> 2 references),
    # so it is ignored
    expect_equal(res$n, c(4, 1))
    expect_equal(res$counts, matrix(c(4, 0, 0, 0, 2, 0), nrow = 2))
})


test_that(".best_tfclass_family() chooses a single TFClass family", {
    hd <- .resolve_tfclass("Homeodomain")[4]
    c2h2 <- .resolve_tfclass("C2H2 ZF")[4]
    bed <- .resolve_tfclass("BED ZF,C2H2 ZF")[4]
    evidence <- data.frame(
        Gene = c("g1", "g1", "g2", "g2", "g3", "g3", "g4", "g5", "g5"),
        Domain = c(
            "3.1.1", "3.1.2", "3.1.1", "3.1.2", "2.3.3", "2.3.4", "2.3.1",
            "2.3.3", "2.3.4"
        ),
        score = c(50, 30, 50, 45, 0.3, 0.1, 0.05, 50, 10),
        method = c(rep("HMM", 4), rep("k-mer", 3), rep("HMM", 2))
    )
    genes <- c("g1", "g2", "g3", "g4", "g5", "g6")
    cisbp <- c(
        "Homeodomain", "Homeodomain", "C2H2 ZF", "C2H2 ZF", "BED ZF,C2H2 ZF",
        "Forkhead"
    )
    ids <- c(hd, hd, c2h2, c2h2, bed, "3.3.1")
    best <- .best_tfclass_family(genes, cisbp, ids, evidence)

    # g1: gap >= 10 bits; g2: gap < 10 bits; g3: shared >= 0.1; g4: shared
    # < 0.1; g5: no profile HMMs for C2H2 ZF families; g6: single family
    expect_equal(best, c("3.1.1", NA, "2.3.3", NA, NA, NA))
    expect_equal(
        .best_tfclass_family(genes, cisbp, ids, evidence, min_shared = 0.01),
        c("3.1.1", NA, "2.3.3", "2.3.1", NA, NA)
    )
    expect_equal(
        .best_tfclass_family(genes, cisbp, ids, evidence[0, ]),
        rep(NA_character_, 6)
    )
})


test_that("TFs are assigned to single TFClass families", {
    skip_if_not(hmmer_is_installed())

    seq <- Biostrings::readAAStringSet(test_path("fixtures", "cow_tfs.fa.gz"))
    annot <- annotate_domains(seq, scheme = "TFClass")
    evidence <- annot$TFClass
    expect_true(all(evidence$method %in% c("HMM", "k-mer")))
    expect_true(all(evidence$Domain[evidence$Gene == "HOXA1"] %in% c(
        .tfclass_hmm_ids()
    )))
    expect_true(all(evidence$Gene %in% c(
        "HOXA1", "SOX2", "SP1", "ZNF263", "ENSBTAG00000003339"
    )))

    fams <- classify_tfs(annot)
    fams <- setNames(fams$TFClass_family, fams$Gene)
    c2h2 <- strsplit(tfclass_of("C2H2 ZF")[3], ";")[[1]]
    expect_equal(fams[["FOXA1"]], "FOX")
    expect_equal(fams[["HOXA1"]], "HOX")
    expect_equal(fams[["SOX2"]], "SOX-related")
    expect_equal(fams[["SP1"]], c2h2[1])
    expect_equal(fams[["ZNF263"]], c2h2[3])
    expect_equal(fams[["ENSBTAG00000003339"]], tfclass_of("C2H2 ZF")[3])

    # Higher minimum fraction of shared k-mers
    strict <- classify_tfs(annot, min_shared = 0.7)
    strict <- setNames(strict$TFClass_family, strict$Gene)
    expect_equal(strict[["SP1"]], fams[["SP1"]])
    expect_equal(strict[["ZNF263"]], tfclass_of("C2H2 ZF")[3])

    expect_error(classify_tfs(annot, min_shared = 2), "min_shared")
    expect_error(classify_tfs(annot, min_shared = NA), "min_shared")
})
