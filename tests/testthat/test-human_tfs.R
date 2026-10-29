#----Load data------------------------------------------------------------------
data(human_tfs)

#----Start tests----------------------------------------------------------------
test_that("human_tfs has the curated human TFs with classifications", {
    expect_true(is(human_tfs, "data.frame"))
    expect_true(nrow(human_tfs) > 0)
    expect_false(any(duplicated(human_tfs$Gene)))
    expect_false(any(is.na(human_tfs$CISBP_family)))

    # Same variables as the output of classify_tfs() for CIS-BP and TFClass
    cols <- c("Gene", .scheme_columns("CISBP"), .scheme_columns("TFClass"))
    expect_true(all(cols %in% names(human_tfs)))

    # TFClass classifications are available for all TFs with a source
    has_source <- !is.na(human_tfs$TFClass_source)
    expect_equal(has_source, !is.na(human_tfs$TFClass_class))
    expect_true(all(
        human_tfs$TFClass_source[has_source] %in%
            c("TFClass", "Inferred from CIS-BP")
    ))
})


test_that("human_tfs can be used in get_tf_counts()", {
    se <- get_tf_counts(list(Human = human_tfs))
    expect_true(is(se, "SummarizedExperiment"))
    expect_equal(sum(SummarizedExperiment::assay(se)), nrow(human_tfs))

    se <- get_tf_counts(
        list(Human = human_tfs), scheme = "TFClass", level = "class"
    )
    expect_true(sum(SummarizedExperiment::assay(se)) > 0)
})
