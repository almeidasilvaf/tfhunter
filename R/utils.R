#' Package-level environment to cache objects within an R session
#' @noRd
.pkg_env <- new.env(parent = emptyenv())


#' Read and parse hmmsearch output
#'
#' @param path Path to hmmsearch output in --domtblout format.
#'
#' @return A data frame with the 23 columns of the --domtblout format, in
#' the same row order as in the file. If the file has no hits, a data frame
#' with 0 rows is returned.
#'
#' @noRd
.read_hmmsearch <- function(path = NULL) {

    cols <- c(
        "target_name", "target_accession", "tlen",
        "query_name", "query_accession", "qlen",
        "full_evalue", "full_score", "full_bias",
        "dom_n", "dom_of", "c_evalue", "i_evalue", "dom_score", "dom_bias",
        "hmm_from", "hmm_to", "ali_from", "ali_to", "env_from", "env_to",
        "acc", "description"
    )
    num_cols <- setdiff(cols, c(
        "target_name", "target_accession", "query_name", "query_accession",
        "description"
    ))

    lines <- readLines(path)
    lines <- lines[!startsWith(lines, "#") & nzchar(trimws(lines))]

    fields <- strsplit(trimws(lines), "[[:space:]]+")
    mat <- vapply(fields, function(x) {
        c(x[seq_len(22)], paste(x[-seq_len(22)], collapse = " "))
    }, character(23))
    mat <- matrix(mat, ncol = 23, byrow = TRUE, dimnames = list(NULL, cols))

    res <- as.data.frame(mat, stringsAsFactors = FALSE)
    res[num_cols] <- lapply(res[num_cols], as.numeric)
    return(res)
}


#' Run hmmsearch and parse its output
#'
#' @param seq_path Path to a FASTA file with protein sequences.
#' @param hmm_path Path to a file with profile HMMs.
#' @param extra_args Character vector of additional arguments to hmmsearch.
#' @param threads Numeric indicating the number of threads to use.
#'
#' @return A data frame with parsed hmmsearch output as returned
#' by \code{.read_hmmsearch()}.
#' @noRd
.run_hmmsearch <- function(
        seq_path, hmm_path, extra_args = NULL, threads = 1
) {

    out_file <- tempfile(pattern = "hmmsearch", fileext = ".domtblout")
    on.exit(unlink(out_file))

    args <- c(
        "--noali", "--cpu", threads, extra_args,
        "--domtblout", shQuote(out_file), shQuote(hmm_path), shQuote(seq_path)
    )
    status <- system2("hmmsearch", args = args, stdout = FALSE)
    if(!identical(as.integer(status), 0L)) {
        stop("hmmsearch failed for HMM file '", basename(hmm_path), "'.")
    }

    return(.read_hmmsearch(out_file))
}


#' Filter HMMER results based on domain cutoffs for each domain
#'
#' @param hmmer_results A data frame with HMMER results as returned
#' by \code{.read_hmmsearch()}.
#' @param evalue Numeric indicating the E-value threshold to use for domains
#' without pre-defined domain cutoffs. Default: 1e-05.
#'
#' @return The same data frame passed as input, but filtered based
#' on domain cutoffs.
#' @noRd
.filter_hmmer <- function(hmmer_results = NULL, evalue = 1e-05) {

    idx <- match(hmmer_results$query_name, score_cutoff$domain)
    cutoff <- score_cutoff$domaincutoff[idx]

    # Domains with a cutoff are filtered by score; others, by E-value
    pass_score <- !is.na(cutoff) & hmmer_results$dom_score >= cutoff
    pass_evalue <- is.na(cutoff) & hmmer_results$full_evalue < evalue
    keep <- !is.na(idx) & (pass_score | pass_evalue)

    res_filtered <- hmmer_results[keep, ]
    rownames(res_filtered) <- NULL
    return(res_filtered)
}


#' Check if HMMER is installed
#'
#' @return Logical indicating whether HMMER is installed (i.e., whether
#' the program \code{hmmsearch} is in the PATH) or not.
#' @export
#' @rdname hmmer_is_installed
#' @examples
#' hmmer_is_installed()
hmmer_is_installed <- function() {
    installed <- nzchar(Sys.which("hmmsearch"))
    return(unname(installed))
}
