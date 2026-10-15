#----Load data------------------------------------------------------------------
data(gsu)
data(gsu_annotation)
seq <- gsu[1:5]

#----Start tests----------------------------------------------------------------
test_that("annotate_domains() returns a list of domain annotations", {
    skip_if_not(hmmer_is_installed())

    annot <- annotate_domains(seq, threads = 2)
    cols <- c("Gene", "Domain", "qlen", "c_evalue", "hmm_from", "hmm_to")
    expect_equal(names(annot), c("PlantTFDB", "TAPscan", "CISBP"))
    expect_equal(names(annot$PlantTFDB), cols)
    expect_equal(names(annot$TAPscan), cols)
    expect_equal(names(annot$CISBP), cols)
    expect_true(nrow(annot$PlantTFDB) >= 2)
    expect_true(nrow(annot$TAPscan) > 0)
    expect_true(nrow(annot$CISBP) > 0)
    expect_true(all(annot$CISBP$c_evalue < 0.01))
    expect_true(all(annot$CISBP$Domain %in% cisbp_map$pfam))

    ref <- gsu_annotation$TAPscan[gsu_annotation$TAPscan$Gene %in% names(seq), ]
    key <- function(x) sort(paste(x$Gene, x$Domain, x$hmm_from, x$hmm_to))
    expect_equal(key(annot$TAPscan), key(ref))

    expect_error(annotate_domains(as.character(seq)), "AAStringSet")
})


test_that("annotate_domains() annotates domains for the chosen schemes", {
    skip_if_not(hmmer_is_installed())

    cisbp <- annotate_domains(seq, scheme = "CISBP")
    expect_equal(names(cisbp), "CISBP")

    tfclass <- annotate_domains(seq, scheme = "PlantTFClass")
    expect_equal(names(tfclass), c("PlantTFDB", "TAPscan"))

    both <- annotate_domains(seq, scheme = c("CISBP", "TAPscan"))
    expect_equal(names(both), c("TAPscan", "CISBP"))

    expect_error(annotate_domains(seq, scheme = "fake"), "scheme")
})


test_that("annotate_domains() fails if HMMER is not installed", {
    path <- Sys.getenv("PATH")
    Sys.setenv(PATH = "")
    installed <- hmmer_is_installed()
    err <- tryCatch(annotate_domains(seq), error = conditionMessage)
    Sys.setenv(PATH = path)

    expect_false(installed)
    expect_match(err, "Could not find HMMER")
})


test_that(".resolve_overlaps() keeps the most significant overlapping hits", {
    hits <- data.frame(
        target_name = c("g1", "g1", "g1", "g1", "g2", "g2", "g3"),
        query_name = c("A", "B", "C", "A", "A", "A", "B"),
        i_evalue = c(1e-5, 1e-10, 1e-3, 1e-4, 1e-3, 1e-6, 1e-2),
        dom_score = c(20, 40, 10, 15, 10, 25, 5),
        ali_from = c(10, 50, 200, 120, 1, 80, 1),
        ali_to = c(60, 100, 250, 205, 100, 150, 50)
    )
    res <- .resolve_overlaps(hits)

    # g1: A (10-60) overlaps the more significant B (50-100) and is removed;
    # C (200-250) overlaps A (120-205), which is more significant and kept;
    # g2: overlapping hits of the same domain are kept
    expect_equal(res$target_name, c("g1", "g1", "g2", "g2", "g3"))
    expect_equal(res$query_name, c("B", "A", "A", "A", "B"))
    expect_equal(res$ali_from, c(50, 120, 1, 80, 1))

    # Order of rows does not change results
    shuffled <- .resolve_overlaps(hits[c(7, 3, 1, 6, 2, 5, 4), ])
    key <- function(x) sort(paste(x$target_name, x$query_name, x$ali_from))
    expect_equal(key(shuffled), key(res))

    expect_equal(.resolve_overlaps(hits[1, ]), hits[1, ])
    expect_equal(nrow(.resolve_overlaps(hits[0, ])), 0)
})
