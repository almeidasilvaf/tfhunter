#----Load data------------------------------------------------------------------
data(gsu)
seq <- gsu[1:5]
cols <- c("Gene", "Domain", "qlen", "c_evalue", "hmm_from", "hmm_to")

# Helper function to check that an object is a valid domain annotation
expect_annotation <- function(annot, sets) {
    expect_true(is(annot, "list"))
    expect_equal(names(annot), sets)
    lapply(annot, function(x) {
        expect_true(is(x, "data.frame"))
        expect_true(all(cols %in% names(x)))
    })
}

# Helper function to check whether hits of different domains overlap
has_overlaps <- function(hits) {
    pairs <- merge(hits, hits, by = "target_name")
    overlap <- pairs$query_name.x != pairs$query_name.y &
        pairs$ali_from.x <= pairs$ali_to.y &
        pairs$ali_to.x >= pairs$ali_from.y
    return(any(overlap))
}

#----Start tests----------------------------------------------------------------
test_that("annotate_domains() returns a list of domain annotations", {
    skip_if_not(hmmer_is_installed())

    annot <- annotate_domains(seq, threads = 2)
    expect_annotation(annot, .annotation_sets(list_schemes()$Scheme))
    expect_true(nrow(annot$TAPscan) > 0)
    expect_true(nrow(annot$CISBP) > 0)
    expect_true(all(annot$CISBP$Gene %in% names(seq)))
    expect_true(all(annot$CISBP$c_evalue < 0.01))
    expect_true(all(annot$CISBP$Domain %in% cisbp_map$pfam))
    expect_true(all(c("score", "method") %in% names(annot$TFClass)))

    expect_error(annotate_domains(as.character(seq)), "AAStringSet")
})


test_that("annotate_domains() annotates domains for the chosen schemes", {
    skip_if_not(hmmer_is_installed())

    lapply(list("CISBP", "PlantTFClass", c("CISBP", "TAPscan"), "TFClass"),
           function(scheme) {
        annot <- annotate_domains(seq, scheme = scheme)
        expect_annotation(annot, .annotation_sets(scheme))
    })

    expect_error(annotate_domains(seq, scheme = "fake"), "scheme")
})


test_that("annotate_domains() filters and formats hits of each set", {
    local_mocked_bindings(
        hmmer_is_installed = function() TRUE,
        .run_hmmsearch = function(seq_path, hmm_path, ...) {
            set <- sub("[0-9a-f]+\\.hmm$", "", basename(hmm_path))
            switch(
                set,
                # g1: score far above PlantTFDB's cutoff; g2: far below it;
                # g3: domain without cutoff, filtered by E-value
                PlantTFDB = fake_hits(
                    c("g1", "g2", "g3"), c("PF00046", "PF00046", "PF00076"),
                    full_evalue = c(1, 1, 1e-50), dom_score = c(500, 0, 0)
                ),
                TAPscan = fake_hits("g1", "Homeobox"),
                # g1: significant; g2: conditional E-value above threshold
                CISBP = fake_hits(
                    c("g1", "g2"), "PF00046", c_evalue = c(1e-50, 1)
                ),
                TFClass = fake_hits("g1", c("3.1.1", "3.1.2"))
            )
        }
    )

    seq <- Biostrings::AAStringSet(c(g1 = "MKV", g2 = "MKV", g3 = "MKV"))
    annot <- annotate_domains(seq)
    expect_annotation(annot, .annotation_sets(list_schemes()$Scheme))
    expect_true(all(c("g1", "g3") %in% annot$PlantTFDB$Gene))
    expect_false("g2" %in% annot$PlantTFDB$Gene)
    expect_true("g1" %in% annot$CISBP$Gene)
    expect_false("g2" %in% annot$CISBP$Gene)
    expect_true(nrow(annot$TFClass) > 0)

    expect_error(annotate_domains(as.character(seq)), "AAStringSet")
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


test_that(".resolve_overlaps() removes overlapping hits of different domains", {
    hits <- data.frame(
        target_name = c("g1", "g1", "g1", "g1", "g2", "g2", "g3"),
        query_name = c("A", "B", "C", "A", "A", "A", "B"),
        i_evalue = c(1e-5, 1e-10, 1e-3, 1e-4, 1e-3, 1e-6, 1e-2),
        dom_score = c(20, 40, 10, 15, 10, 25, 5),
        ali_from = c(10, 50, 200, 120, 1, 80, 1),
        ali_to = c(60, 100, 250, 205, 100, 150, 50)
    )
    expect_true(has_overlaps(hits))
    res <- .resolve_overlaps(hits)

    expect_true(is(res, "data.frame"))
    expect_equal(names(res), names(hits))
    expect_true(nrow(res) > 0 && nrow(res) < nrow(hits))
    expect_false(has_overlaps(res))

    # The most significant hit of each gene is always kept
    best <- hits[order(hits$i_evalue), ]
    best <- best[!duplicated(best$target_name), ]
    expect_true(all(paste(best$target_name, best$ali_from) %in%
                        paste(res$target_name, res$ali_from)))

    # Overlapping hits of the same domain are kept
    expect_true(sum(res$target_name == "g2") == sum(hits$target_name == "g2"))

    # Order of rows does not change results
    shuffled <- .resolve_overlaps(hits[rev(seq_len(nrow(hits))), ])
    key <- function(x) sort(paste(x$target_name, x$query_name, x$ali_from))
    expect_equal(key(shuffled), key(res))

    expect_equal(.resolve_overlaps(hits[1, ]), hits[1, ])
    expect_equal(nrow(.resolve_overlaps(hits[0, ])), 0)
})
