#----Load data------------------------------------------------------------------
data(tfclass_scheme)
ev_cols <- c(
    "Gene", "Domain", "qlen", "c_evalue", "hmm_from", "hmm_to", "score",
    "method"
)

# Helper function to get the TFClass classification of a CIS-BP family
tfclass_of <- function(cisbp_family) {
    res <- .classify_tfclass(
        data.frame(Gene = "g1", CISBP_family = cisbp_family)
    )
    return(unlist(res[1, .scheme_columns("TFClass")], use.names = FALSE))
}

# Helper function to split ';'-separated values
split_values <- function(x) strsplit(x, ";", fixed = TRUE)[[1]]

# Helper function to check that TFClass families are valid
expect_valid_families <- function(families) {
    values <- unlist(strsplit(families[!is.na(families)], ";", fixed = TRUE))
    expect_true(all(values %in% .tfclass_families()$family))
}

#----Start tests----------------------------------------------------------------
test_that(".classify_tfclass() uses the correspondence table", {
    # CIS-BP families corresponding to a single TFClass family
    fox <- tfclass_of("Forkhead")
    expect_false(any(is.na(fox)))
    expect_false(any(grepl(";", fox)))

    # CIS-BP families corresponding to multiple TFClass families
    hd <- tfclass_of("Homeodomain")
    expect_false(any(grepl(";", hd[1:2])))
    expect_true(length(split_values(hd[3])) > 1)
    expect_valid_families(hd[3])

    # All CIS-BP families in the table have TFClass families
    expect_false(any(is.na(tfclass_scheme$TFClass_family)))
    fams <- lapply(tfclass_scheme$CISBP_family, tfclass_of)
    expect_false(any(is.na(unlist(fams))))

    # CIS-BP families not classified by TFClass
    expect_true(all(is.na(tfclass_of("Zinc cluster"))))
    expect_true(all(is.na(tfclass_of("AP2,B3"))))
    expect_true(all(is.na(tfclass_of(NA))))
})


test_that(".classify_tfclass() resolves other families with multiple DBDs", {
    # All families of each DBD (DBDs without TFClass classes are ignored)
    combos <- list(
        c("BED ZF", "C2H2 ZF"), c("bHLH", "T-box"), c("AT hook", "Homeodomain"),
        c("bHLH", "Zinc cluster")
    )
    lapply(combos, function(dbds) {
        res <- tfclass_of(paste(dbds, collapse = ","))
        parts <- lapply(dbds, tfclass_of)
        parts <- parts[!vapply(parts, function(x) all(is.na(x)), logical(1))]
        expected <- unique(unlist(lapply(parts, function(x) {
            split_values(x[3])
        })))
        expect_true(setequal(split_values(res[3]), expected))
    })

    # DBDs in different TFClass classes: all classes, separated by ';'
    expect_true(grepl(";", tfclass_of("bHLH,T-box")[2]))
})


test_that(".classify_tfclass() handles multiple genes", {
    fams <- data.frame(
        Gene = c("g1", "g2", "g3", "g4"),
        CISBP_family = c("Forkhead", "WRKY", "Forkhead", NA)
    )
    res <- .classify_tfclass(fams)
    expect_true(is(res, "data.frame"))
    expect_equal(names(res), c(names(fams), .scheme_columns("TFClass")))
    expect_equal(nrow(res), nrow(fams))
    expect_equal(res$TFClass_family[1], res$TFClass_family[3])
    expect_false(is.na(res$TFClass_family[1]))
    expect_true(all(is.na(res$TFClass_family[c(2, 4)])))
})


test_that("TFClass correspondence table is consistent", {
    expect_true(is(tfclass_scheme, "data.frame"))
    expect_true(nrow(tfclass_scheme) > 0)
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
    expect_true(is(fams, "data.frame"))
    expect_equal(names(fams), c("superclass", "class", "family", "id"))
    expect_true(nrow(fams) > 0)
    expect_false(any(duplicated(fams$id)))
    expect_false(any(grepl(";", unlist(fams))))
    expect_true(all(.tfclass_hmm_ids() %in% fams$id))
    expect_true(all(split_values(.resolve_tfclass("C2H2 ZF")[4]) %in% fams$id))
})


test_that(".kmer_table() encodes unique k-mers of each sequence", {
    seq <- Biostrings::AAStringSet(c(a = "ACDEFA", b = "ACDXACD", c = "AC"))
    km <- .kmer_table(seq, k = 3)

    expect_true(is(km, "data.frame"))
    expect_equal(names(km), c("code", "seq"))
    expect_true(nrow(km) > 0)
    expect_true(all(km$code >= 0))
    expect_true(all(km$seq %in% seq_along(seq)))

    # Unique k-mers per sequence; k-mers with non-standard amino acids and
    # sequences shorter than k are ignored
    expect_false(any(duplicated(km)))
    n <- tabulate(km$seq, nbins = length(seq))
    expect_true(all(n <= pmax(0, Biostrings::width(seq) - 3 + 1)))
    expect_true(n[2] < Biostrings::width(seq)[2] - 3 + 1)
    expect_equal(n[3], 0)
    expect_equal(nrow(.kmer_table(seq["c"], k = 3)), 0)

    # Same k-mers have the same code
    expect_true(any(km$code[km$seq == 1] %in% km$code[km$seq == 2]))
})


test_that(".shared_kmers() counts k-mers shared by queries and references", {
    ref <- Biostrings::AAStringSet(c(
        r1 = "MKVLAAGIRS", r2 = "PPWWQQHHEE", r3 = "MKVLAWWPPQ"
    ))
    query <- Biostrings::AAStringSet(c(q1 = "MKVLAAGIRS", q2 = "TTTTT"))
    res <- .shared_kmers(query, ref, k = 3)

    expect_true(is(res, "list"))
    expect_equal(dim(res$counts), c(length(query), length(ref)))
    expect_equal(length(res$n), length(query))
    expect_true(all(res$counts >= 0))
    expect_true(all(res$counts <= res$n))

    # A query identical to a reference shares all its k-mers with it, and
    # this is the most similar reference
    expect_equal(res$counts[1, 1], res$n[1])
    expect_equal(which.max(res$counts[1, ]), 1)
    expect_true(all(res$counts[2, ] == 0))

    # K-mers found in many references are ignored
    many <- .shared_kmers(query, ref[c(1, 1, 1)], k = 3)
    expect_true(all(many$counts == 0))
})


test_that(".best_tfclass_family() chooses a single TFClass family", {
    hd <- .resolve_tfclass("Homeodomain")[4]
    c2h2 <- .resolve_tfclass("C2H2 ZF")[4]
    bed <- .resolve_tfclass("BED ZF,C2H2 ZF")[4]
    hd1 <- split_values(hd)[1:2]
    c2h21 <- split_values(c2h2)[1:2]
    evidence <- data.frame(
        Gene = c("g1", "g1", "g2", "g2", "g3", "g3", "g4", "g5", "g5"),
        Domain = c(hd1, hd1, c2h21, c2h21[1], c2h21),
        score = c(500, 1, 50, 49, 0.9, 0.01, 0.001, 500, 1),
        method = c(rep("HMM", 4), rep("k-mer", 3), rep("HMM", 2))
    )
    genes <- c("g1", "g2", "g3", "g4", "g5", "g6")
    cisbp <- c(
        "Homeodomain", "Homeodomain", "C2H2 ZF", "C2H2 ZF", "BED ZF,C2H2 ZF",
        "Forkhead"
    )
    ids <- c(hd, hd, c2h2, c2h2, bed, .resolve_tfclass("Forkhead")[4])
    best <- .best_tfclass_family(genes, cisbp, ids, evidence)
    expect_true(is.character(best))
    expect_equal(length(best), length(genes))

    # Assigned: clear best HMM (g1) or high fraction of shared k-mers (g3)
    expect_true(best[1] %in% split_values(hd))
    expect_true(best[3] %in% split_values(c2h2))

    # Not assigned: similar HMM scores (g2), low fraction of shared k-mers
    # (g4), no profile HMMs for all families (g5), single family (g6)
    expect_true(all(is.na(best[c(2, 4, 5, 6)])))

    # Lower thresholds never assign fewer genes
    lenient <- .best_tfclass_family(
        genes, cisbp, ids, evidence, min_shared = 0, min_gap = 0
    )
    expect_true(all(!is.na(lenient[!is.na(best)])))
    expect_true(sum(!is.na(lenient)) >= sum(!is.na(best)))

    # Without evidence (or only evidence of other methods)
    none <- .best_tfclass_family(genes, cisbp, ids, evidence[0, ])
    expect_true(all(is.na(none)))
    other <- evidence[evidence$Gene %in% c("g1", "g3"), ]
    other$method <- ifelse(other$method == "HMM", "k-mer", "HMM")
    expect_true(all(is.na(.best_tfclass_family(genes, cisbp, ids, other))))
})


test_that("TFs are assigned to single TFClass families", {
    skip_if_not(hmmer_is_installed())

    seq <- Biostrings::readAAStringSet(test_path("fixtures", "cow_tfs.fa.gz"))
    annot <- annotate_domains(seq, scheme = "TFClass")
    evidence <- annot$TFClass
    expect_true(is(evidence, "data.frame"))
    expect_true(all(ev_cols %in% names(evidence)))
    expect_true(nrow(evidence) > 0)
    expect_true(all(evidence$method %in% c("HMM", "k-mer")))
    expect_true(all(evidence$Gene %in% names(seq)))
    expect_true(all(evidence$Domain %in% .tfclass_families()$id))

    fams <- classify_tfs(annot)
    expect_true(is(fams, "data.frame"))
    expect_true(nrow(fams) > 0)
    expect_valid_families(fams$TFClass_family)

    # Evidence reduces the number of genes with multiple families
    no_evidence <- .classify_tfclass(fams[, c("Gene", "CISBP_family")])
    multi <- function(x) sum(grepl(";", x$TFClass_family))
    expect_true(multi(fams) < multi(no_evidence))

    # Higher minimum fraction of shared k-mers never assigns more genes
    strict <- classify_tfs(annot, min_shared = 1)
    expect_true(multi(strict) >= multi(fams))
})


test_that(".empty_tfclass_evidence() returns an empty evidence table", {
    ev <- .empty_tfclass_evidence()
    expect_true(is(ev, "data.frame"))
    expect_equal(nrow(ev), 0)
    expect_equal(names(ev), ev_cols)
})


test_that(".tfclass_c2h2_refs() reads reference C2H2 ZF TFs once", {
    ref <- .tfclass_c2h2_refs()
    expect_true(is(ref, "AAStringSet"))
    expect_true(length(ref) > 0)

    # Names start with IDs of TFClass families of C2H2 ZF TFs
    ids <- sub("\\|.*", "", names(ref))
    expect_true(all(ids %in% split_values(.resolve_tfclass("C2H2 ZF")[4])))
    expect_identical(.tfclass_c2h2_refs(), ref)
})


test_that(".tfclass_kmer_evidence() scores families by shared k-mers", {
    local_mocked_bindings(.tfclass_c2h2_refs = function() {
        Biostrings::AAStringSet(c(
            `2.3.1|r1` = "MKVLAAGIRSTW", `2.3.1|r2` = "MKVLAAGPPPPP",
            `2.3.3|r3` = "WWWHHHEEEKKK"
        ))
    })
    seq <- Biostrings::AAStringSet(c(
        `q1 a description` = "MKVLAAGIRS", q2 = "WWWHHHEEECCC",
        q3 = "NNNNNNNNNN"
    ))
    ev <- .tfclass_kmer_evidence(seq)

    expect_true(is(ev, "data.frame"))
    expect_equal(names(ev), ev_cols)
    expect_true(nrow(ev) > 0)
    expect_true(all(ev$method == "k-mer"))
    expect_true(all(ev$score > 0 & ev$score <= 1))

    # Gene IDs without descriptions; genes without shared k-mers are absent
    expect_true(all(ev$Gene %in% c("q1", "q2")))
    expect_false("q3" %in% ev$Gene)

    # Highest score for the family of the most similar reference
    best <- function(g) {
        x <- ev[ev$Gene == g, ]
        x$Domain[which.max(x$score)]
    }
    expect_equal(best("q1"), "2.3.1")
    expect_equal(best("q2"), "2.3.3")

    expect_equal(.tfclass_kmer_evidence(seq[0]), .empty_tfclass_evidence())
})


test_that(".tfclass_hmm_evidence() keeps the best hit per gene and family", {
    local_mocked_bindings(.run_hmmsearch = function(...) {
        fake_hits(
            target_name = c("g1", "g1", "g1", "g1", "g2"),
            query_name = c("3.1.1", "3.1.1", "3.1.2", "3.1.3", "3.1.1"),
            full_evalue = c(1e-20, 1e-20, 1e-10, 1, 1e-5),
            dom_score = c(30, 60, 40, 5, 20),
            hmm_from = c(1, 2, 1, 1, 1)
        )
    })
    seq <- Biostrings::AAStringSet(c(g1 = "MKV", g2 = "MKV"))
    ev <- .tfclass_hmm_evidence(seq)

    expect_true(is(ev, "data.frame"))
    expect_equal(names(ev), ev_cols)
    expect_true(nrow(ev) > 0)
    expect_true(all(ev$method == "HMM"))
    expect_true(is.numeric(ev$score))

    # One row per gene and family, without non-significant hits
    expect_false(any(duplicated(paste(ev$Gene, ev$Domain))))
    expect_false(any(ev$Gene == "g1" & ev$Domain == "3.1.3"))

    # Coordinates of the best domain hit
    g1 <- ev[ev$Gene == "g1" & ev$Domain == "3.1.1", ]
    expect_equal(g1$hmm_from, 2)

    expect_equal(.tfclass_hmm_evidence(seq[0]), .empty_tfclass_evidence())
})


test_that(".annotate_tfclass() gets evidence for ambiguous TFs only", {
    local_mocked_bindings(
        .tfclass_hmm_evidence = function(seq, threads) {
            data.frame(Gene = names(seq), method = rep("HMM", length(seq)))
        },
        .tfclass_kmer_evidence = function(seq) {
            data.frame(Gene = names(seq), method = rep("k-mer", length(seq)))
        }
    )
    # g1: Homeodomain; g2: C2H2 ZF; g3: Forkhead; g4: BED ZF,C2H2 ZF
    cisbp <- data.frame(
        Gene = c("g1", "g2", "g3", "g4", "g4"),
        Domain = c("PF00046", "PF00096", "PF00250", "PF00096", "PF02892")
    )
    seq <- Biostrings::AAStringSet(c(
        `g1 homeodomain` = "MKV", g2 = "MKV", g3 = "MKV", g4 = "MKV"
    ))
    ev <- .annotate_tfclass(seq, cisbp)
    expect_true(is(ev, "data.frame"))
    genes <- sub(" .*", "", ev$Gene)

    # Profile HMMs for ambiguous families with profile HMMs, k-mers for
    # C2H2 ZF, and nothing for single families or families without HMMs
    expect_equal(ev$method[genes == "g1"], "HMM")
    expect_equal(ev$method[genes == "g2"], "k-mer")
    expect_false(any(c("g3", "g4") %in% genes))
})


test_that("classify_tfs() uses TFClass evidence", {
    data(gsu_annotation)
    annot <- gsu_annotation[c("CISBP", "TFClass")]
    gene <- annot$CISBP$Gene[annot$CISBP$Domain == "PF00046"][1]
    hd <- split_values(.resolve_tfclass("Homeodomain")[4])
    annot$TFClass <- data.frame(
        Gene = gene, Domain = hd[1:2], qlen = 60, c_evalue = 1e-10,
        hmm_from = 1, hmm_to = 60, score = c(500, 1), method = "HMM"
    )
    fams <- classify_tfs(annot, scheme = "TFClass")
    res <- fams[fams$Gene == gene, ]
    expected <- .tfclass_families()[.tfclass_families()$id == hd[1], ]
    expect_equal(res$TFClass_family, expected$family)
    expect_equal(res$TFClass_class, expected$class)
    expect_equal(res$TFClass_superclass, expected$superclass)
})
