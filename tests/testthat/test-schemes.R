#----Start tests----------------------------------------------------------------
test_that("list_schemes() returns a data frame of schemes", {
    schemes <- list_schemes()
    expect_true(is(schemes, "data.frame"))
    expect_true(nrow(schemes) > 0)
    expect_true(all(
        c("Scheme", "Taxa", "Levels", "Annotation", "Reference") %in%
            names(schemes)
    ))
    expect_false(any(duplicated(schemes$Scheme)))

    # Levels match the variables in the output of classify_tfs()
    levels <- lapply(schemes$Scheme, function(s) names(.scheme_levels(s)))
    expect_equal(
        vapply(levels, paste, character(1), collapse = ", "), schemes$Levels
    )
})


test_that(".check_scheme() validates and expands schemes", {
    all_schemes <- list_schemes()$Scheme
    expect_equal(.check_scheme("all"), all_schemes)

    # Schemes are returned in the order of list_schemes()
    two <- .check_scheme(rev(all_schemes[1:2]))
    expect_equal(two, all_schemes[1:2])

    # 'all' is expanded to schemes with all their annotations available
    lapply(all_schemes, function(s) {
        expanded <- .check_scheme("all", .annotation_sets(s))
        expect_true(s %in% expanded)
        needs <- lapply(expanded, .annotation_sets)
        expect_true(all(unlist(needs) %in% .annotation_sets(s)))
    })

    expect_error(.check_scheme("fake"), "must be 'all'")
    expect_error(.check_scheme(character(0)), "must be 'all'")
    expect_error(.check_scheme(1), "must be 'all'")
    expect_error(.check_scheme("CISBP", "TAPscan"), "Missing")
    expect_error(.check_scheme("TFClass", "CISBP"), "Missing")
    expect_error(.check_scheme("all", character(0)), "No classification")
})


test_that(".annotation_sets() returns annotations required by schemes", {
    schemes <- list_schemes()
    lapply(seq_len(nrow(schemes)), function(i) {
        sets <- .annotation_sets(schemes$Scheme[i])
        expect_true(is.character(sets) && length(sets) > 0)
        expect_true(all(sets %in% strsplit(schemes$Annotation[i], ", ")[[1]]))
    })
    all_sets <- .annotation_sets(schemes$Scheme)
    expect_false(any(duplicated(all_sets)))
})


test_that(".scheme_hmm_path() decompresses HMMs once per session", {
    lapply(.annotation_sets(list_schemes()$Scheme), function(set) {
        path <- .scheme_hmm_path(set)
        expect_true(file.exists(path))
        expect_equal(.scheme_hmm_path(set), path)
        expect_true(startsWith(readLines(path, n = 1), "HMMER3"))
    })

    # TFClass: one profile HMM per TFClass family (except C2H2 ZF families)
    hmm <- readLines(.scheme_hmm_path("TFClass"))
    hmm_names <- sub("^NAME +", "", hmm[startsWith(hmm, "NAME")])
    expect_true(setequal(hmm_names, .tfclass_hmm_ids()))
})
