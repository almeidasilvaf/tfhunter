#----Load data------------------------------------------------------------------
## Output of hmmsearch against TAPscan domains for proteins of Ectocarpus sp. 7
domtbl <- .read_hmmsearch(test_path("fixtures", "ectocarpus.domtblout.gz"))
ecto_annot <- data.frame(
    Gene = domtbl$target_name, Domain = domtbl$query_name, qlen = domtbl$qlen,
    c_evalue = domtbl$c_evalue, hmm_from = domtbl$hmm_from,
    hmm_to = domtbl$hmm_to
)

## Families obtained with the original TAPscan code (v4.76)
expected <- read.table(
    test_path("fixtures", "ectocarpus_tapscan_families.txt.gz"),
    sep = ";", col.names = c("Gene", "Family", "Subfamily")
)

# Families without subfamilies ('-' in TAPscan's output)
no_sub <- expected$Subfamily == "-"
expected$Subfamily[no_sub] <- expected$Family[no_sub]

#----Start tests----------------------------------------------------------------
test_that(".classify_tapscan() agrees with the original TAPscan results", {
    fams <- .classify_tapscan(ecto_annot)
    expect_true(is(fams, "data.frame"))
    expect_equal(names(fams), c("Gene", "Family", "Subfamily"))
    expect_true(nrow(fams) > 0)
    expect_false(any(duplicated(fams$Gene)))

    merged <- merge(fams, expected, by = "Gene", all = TRUE)
    agree <- merged$Family.x %in% merged$Family.y &
        merged$Subfamily.x %in% merged$Subfamily.y
    expect_true(mean(agree) > 0.95)
})


test_that(".tapscan_filter_coverage() filters hits by coverage", {
    annot <- data.frame(
        Gene = c("g1", "g2", "g3"), Domain = c("AP2", "AP2", "unknown"),
        qlen = 100, c_evalue = 1e-10, hmm_from = 1, hmm_to = c(2, 100, 1)
    )
    hits <- .tapscan_filter_coverage(annot)
    expect_true(is(hits, "data.frame"))
    expect_true(all(c("Gene", "Domain", "c_evalue") %in% names(hits)))

    # Low coverage: removed; domains without coverage cutoffs: kept
    expect_false("g1" %in% hits$Gene)
    expect_true(all(c("g2", "g3") %in% hits$Gene))
})


test_that(".tapscan_collapse_hits() collapses hits of the same domain", {
    hits <- data.frame(
        Gene = c("g2", "g2", "g2", "g1", "g1"),
        Domain = c("Myb_DNA-binding", "Myb_DNA-binding", "Myb_DNA-binding",
                   "AP2", "B3"),
        c_evalue = c(1e-05, 1e-10, 1e-08, 1e-03, 1e-20)
    )
    collapsed <- .tapscan_collapse_hits(hits)
    expect_true(is(collapsed, "data.frame"))
    expect_true("count" %in% names(collapsed))

    # One row per gene and domain, with the number of hits
    key <- paste(collapsed$Gene, collapsed$Domain)
    expect_false(any(duplicated(key)))
    expect_true(setequal(key, paste(hits$Gene, hits$Domain)))
    expect_equal(sum(collapsed$count), nrow(hits))
    expect_true(all(collapsed$c_evalue %in% hits$c_evalue))
})


test_that(".tapscan_exclude_similar() keeps only the best similar domain", {
    hits <- data.frame(
        Gene = c("g1", "g1", "g1", "g1", "g2", "g2"),
        Domain = c("G2-like_Domain", "Myb_DNA-binding", "WD40",
                   "FIE_clipped_for_HMM", "zf-Dof", "GATA"),
        c_evalue = 1e-10, count = 1
    )
    filtered <- .tapscan_exclude_similar(hits)
    expect_true(is(filtered, "data.frame"))
    expect_equal(names(filtered), names(hits))
    expect_true(nrow(filtered) > 0 && nrow(filtered) < nrow(hits))
    expect_true(all(filtered$Domain %in% hits$Domain))
    expect_true(all(c("g1", "g2") %in% filtered$Gene))
})


test_that(".tapscan_tokens() adds pseudo-domains", {
    hits <- data.frame(
        Gene = c("g1", "g1", "g2", "g3"),
        Domain = c("Myb_DNA-binding", "LIM", "Myb_DNA-binding", "LIM"),
        c_evalue = 1e-10, count = c(3, 2, 5, 1)
    )
    tokens <- .tapscan_tokens(hits)
    expect_true(is(tokens, "data.frame"))
    expect_true(all(c("Gene", "Token") %in% names(tokens)))

    # All domains are kept as tokens, and pseudo-domains are added
    expect_true(all(paste(hits$Gene, hits$Domain) %in%
                        paste(tokens$Gene, tokens$Token)))
    expect_true(nrow(tokens) > nrow(hits))
    expect_true(setequal(tokens$Gene, hits$Gene))
})


test_that(".tapscan_tiebreak() chooses a single family", {
    expect_equal(.tapscan_tiebreak("A", ";d1;", "d1"), "A")

    res <- .tapscan_tiebreak(c("A", "B"), c(";d1;", ";d2;d1;"), c("d2", "d1"))
    expect_true(length(res) == 1 && res %in% c("A", "B"))

    # Tokens shared by both candidates: keep the first candidate
    expect_equal(
        .tapscan_tiebreak(c("A", "B"), c(";d1;", ";d1;"), "d1"), "A"
    )
})


test_that(".tapscan_family_map() merges families and adds subfamilies", {
    input <- c("bZIP2", "MYB-2R", "NAC", "HRT")
    fams <- .tapscan_family_map(input)
    expect_true(is(fams, "data.frame"))
    expect_equal(names(fams), c("Family", "Subfamily"))
    expect_equal(nrow(fams), length(input))
    expect_false(any(is.na(fams$Family)))
    expect_false(any(is.na(fams$Subfamily)))
})


test_that(".classify_tapscan() handles input without classified genes", {
    annot <- data.frame(
        Gene = "g1", Domain = "AP2", qlen = 100, c_evalue = 1e-10,
        hmm_from = 1, hmm_to = 2
    )
    fams <- .classify_tapscan(annot)
    expect_true(is(fams, "data.frame"))
    expect_equal(nrow(fams), 0)
    expect_equal(names(fams), c("Gene", "Family", "Subfamily"))

    # Domain that passes the coverage filter, but is not enough for a family
    # (Argonaute requires PAZ and Piwi)
    annot$Domain <- "PAZ"
    annot$hmm_to <- 100
    fams <- .classify_tapscan(annot)
    expect_equal(nrow(fams), 0)
    expect_equal(names(fams), c("Gene", "Family", "Subfamily"))
})


test_that(".add_tap_class() adds TAP classes", {
    fams <- data.frame(
        TAPscan_family = c("MYB", "HAT", "NAC", "fake"),
        TAPscan_subfamily = c("MYB-2R", "GNAT", "NAC", "fake")
    )
    res <- .add_tap_class(fams)
    expect_true("TAP_class" %in% names(res))
    expect_true(all(res$TAP_class %in% c("TF", "TR", "PT", NA)))
    expect_true(is.na(res$TAP_class[4]))
})


test_that("tap_class data set has classes for TAPscan families", {
    data(tap_class)
    expect_true(is(tap_class, "data.frame"))
    expect_true(nrow(tap_class) > 0)
    expect_true(all(c("TAP_family", "TAP_class") %in% names(tap_class)))
    expect_true(all(tap_class$TAP_class %in% c("TF", "TR", "PT")))
    expect_false(any(duplicated(tap_class$TAP_family)))
    expect_equal(tap_class[, 1:2], tap_class_map)

    # Most families produced by TAPscan rules have a TAP class
    fams <- unique(.tapscan_family_map(unique(tapscan_rules$family)))
    names(fams) <- c("TAPscan_family", "TAPscan_subfamily")
    fams <- .add_tap_class(fams)
    expect_true(mean(!is.na(fams$TAP_class)) > 0.9)
})
