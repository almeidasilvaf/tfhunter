#---Load and create data for testing--------------------------------------------
## Example output of hmmsearch (.txt format)
example_output <- c(
    "#                                                                            --- full sequence --- -------------- this domain -------------   hmm coord   ali coord   env coord",
    "# target name        accession   tlen query name           accession   qlen   E-value  score  bias   #  of  c-Evalue  i-Evalue  score  bias  from    to  from    to  from    to  acc description of target",
    "#------------------- ---------- ----- -------------------- ---------- ----- --------- ------ ----- --- --- --------- --------- ------ ----- ----- ----- ----- ----- ----- ----- ---- ---------------------",
    "gsu06140.1           -            340 G2-like              -             56   5.5e-09   24.7   0.0   1   1   2.1e-09   1.1e-08   23.8   0.0     1    53    95   141    95   144 0.89 -",
    "gsu06160.1           -            120 NF-YC                -            102   1.1e-07   20.5   0.0   1   1   2.7e-08   1.3e-07   20.2   0.0    19    90    18    89     7    98 0.90 a description with spaces"
)
file <- file.path(tempdir(), "example_file_unit_test.txt")
writeLines(example_output, file)

empty_file <- file.path(tempdir(), "example_empty_file_unit_test.txt")
writeLines(example_output[1:3], empty_file)

#----Start tests----------------------------------------------------------------
test_that(".read_hmmsearch() reads and parses hmmsearch output", {
    output <- .read_hmmsearch(file)

    expect_true(is(output, "data.frame"))
    expect_equal(ncol(output), 23)
    expect_true(nrow(output) > 0)
    expect_true(is.character(output$target_name))
    expect_true(is.numeric(output$dom_score))
    expect_true(is.numeric(output$c_evalue))

    # Descriptions with spaces are kept in a single variable
    expect_true(any(grepl(" ", output$description)))

    empty <- .read_hmmsearch(empty_file)
    expect_true(is(empty, "data.frame"))
    expect_equal(nrow(empty), 0)
})


test_that(".filter_hmmer() filters hmmsearch output", {
    output <- .read_hmmsearch(file)
    foutput <- .filter_hmmer(output)

    expect_true(is(foutput, "data.frame"))
    expect_equal(names(foutput), names(output))
    expect_true(nrow(foutput) <= nrow(output))

    # Domains without pre-defined cutoffs are filtered by E-value
    output$query_name <- "PF00076"
    lenient <- .filter_hmmer(output, evalue = 1)
    strict <- .filter_hmmer(output, evalue = 1e-100)
    expect_true(nrow(lenient) > 0)
    expect_equal(nrow(strict), 0)
})


test_that("hmmer_is_installed() works", {
    h <- hmmer_is_installed()
    expect_true(is.logical(h))
    expect_equal(length(h), 1)
})


test_that(".run_hmmsearch() runs hmmsearch and parses its output", {
    skip_on_os("windows")

    path <- Sys.getenv("PATH")
    Sys.setenv(PATH = paste(fake_hmmsearch(status = 0), path, sep = ":"))
    hits <- tryCatch(
        .run_hmmsearch("seqs.fasta", "domains.hmm", threads = 2),
        finally = Sys.setenv(PATH = path)
    )
    expect_true(is(hits, "data.frame"))
    expect_equal(names(hits), names(.read_hmmsearch(file)))
    expect_true(nrow(hits) > 0)

    Sys.setenv(PATH = paste(fake_hmmsearch(status = 1), path, sep = ":"))
    err <- tryCatch(
        .run_hmmsearch("seqs.fasta", "domains.hmm"),
        error = conditionMessage, finally = Sys.setenv(PATH = path)
    )
    expect_match(err, "hmmsearch failed")
})
