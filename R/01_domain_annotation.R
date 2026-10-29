#' Annotate protein sequences with domains used for TF classification
#'
#' Protein sequences are scanned with HMMER for the domains used by
#' the classification schemes of choice (see \code{list_schemes()}).
#' Each domain annotation set has its own profile HMMs and filtering
#' criteria:
#' \itemize{
#'   \item \strong{PlantTFDB}: PFAM and self-built profile HMMs, with the
#'   domain-specific score cutoffs defined by PlantTFDB (or an E-value
#'   threshold, for domains without cutoffs).
#'   \item \strong{TAPscan}: TAPscan v4 profile HMMs, with their gathering
#'   (GA) thresholds, as in the original TAPscan implementation.
#'   \item \strong{CISBP}: Pfam profile HMMs of the DNA-binding domains used
#'   by CIS-BP, with the thresholds used by CIS-BP (per-sequence E-value
#'   < 0.01 and per-domain conditional E-value < 0.01). As in CIS-BP,
#'   overlapping hits of different domains are resolved by keeping only
#'   the most significant hit (lowest independent E-value).
#'   \item \strong{TFClass}: evidence used to choose among multiple
#'   TFClass families for TFs in CIS-BP families that correspond to more than
#'   one TFClass family (e.g., CIS-BP's Homeodomain). TFs in the CIS-BP family
#'   C2H2 ZF are compared to human and mouse C2H2 ZF TFs with known TFClass
#'   families based on shared amino acid 5-mers, and other TFs are scanned
#'   with profile HMMs of the DNA-binding domains of each TFClass family
#'   (built from human and mouse TFs). See \code{classify_tfs()} for
#'   details.
#' }
#'
#' @param seq An `AAStringSet` object. The sequences in this object must
#' represent only the translated sequences of primary (or longest) transcripts.
#' To create such an object from a FASTA file,
#' use \code{Biostrings::readAAStringSet()}.
#' @param scheme Character vector indicating the classification scheme(s)
#' for which domains will be annotated. One or more of the schemes in
#' \code{list_schemes()}, or 'all' to annotate domains for all schemes.
#' Default: 'all'.
#' @param evalue Numeric indicating the E-value threshold for \code{hmmsearch}
#' to be used for PlantTFDB domains without pre-defined domain cutoffs.
#' Default: 1e-05.
#' @param threads Numeric indicating the number of threads to be used
#' by \code{hmmsearch}. Default: 1.
#'
#' @return A list of data frames with domain annotation for each domain
#' annotation set required by the chosen scheme(s) (\strong{PlantTFDB},
#' \strong{TAPscan}, \strong{CISBP}, and/or \strong{TFClass}; see variable
#' \strong{Annotation} in \code{list_schemes()}). All data frames have one
#' row per domain hit and the following variables:
#' \describe{
#'   \item{Gene}{Character, gene ID.}
#'   \item{Domain}{Character, domain ID or name.}
#'   \item{qlen}{Numeric, length of the profile HMM.}
#'   \item{c_evalue}{Numeric, conditional E-value of the domain hit.}
#'   \item{hmm_from, hmm_to}{Numeric, start and end coordinates of the alignment
#'   in the profile HMM.}
#' }
#' Rows are kept in the same order as in the HMMER output, as this order
#' is used by the TAPscan classification algorithm.
#'
#' The data frame \strong{TFClass} has one row per gene and candidate TFClass
#' family, with \strong{Domain} indicating the TFClass family ID, and two
#' additional variables:
#' \describe{
#'   \item{score}{Numeric, bit score of the profile HMM of the TFClass
#'   family (for the whole sequence), or fraction of the gene's 5-mers shared
#'   with the most similar reference C2H2 ZF TF in the TFClass family.}
#'   \item{method}{Character, method used to obtain the score ('HMM' or
#'   'k-mer'). For 'k-mer', variables \strong{qlen}, \strong{c_evalue},
#'   \strong{hmm_from}, and \strong{hmm_to} are NA.}
#' }
#'
#' @export
#' @rdname annotate_domains
#' @examples
#' data(gsu)
#' seq <- gsu[1:5]
#' if(hmmer_is_installed()) {
#'     annotate_domains(seq)
#'
#'     # Only domains used by CIS-BP's scheme
#'     annotate_domains(seq, scheme = "CISBP")
#' }
annotate_domains <- function(
        seq = NULL, scheme = "all", evalue = 1e-05, threads = 1
) {

    sets <- .annotation_sets(.check_scheme(scheme))
    seq_path <- .write_seqs(seq)
    on.exit(unlink(seq_path))

    hmm_sets <- setdiff(sets, "TFClass")
    annotation <- lapply(hmm_sets, function(set) {
        annot <- switch(
            set,
            PlantTFDB = .annotate_planttfdb(seq_path, evalue, threads),
            TAPscan = .annotate_tapscan(seq_path, threads),
            CISBP = .annotate_cisbp(seq_path, threads)
        )
        return(annot)
    })
    names(annotation) <- hmm_sets

    # TFClass evidence requires CIS-BP families
    if("TFClass" %in% sets) {
        annotation$TFClass <- .annotate_tfclass(seq, annotation$CISBP, threads)
    }
    return(annotation)
}


#' Write sequences to a temporary FASTA file after checking input
#'
#' @param seq An `AAStringSet` object.
#'
#' @return Path to the temporary FASTA file.
#' @importFrom Biostrings writeXStringSet
#' @importFrom methods is
#' @noRd
.write_seqs <- function(seq = NULL) {

    if(!hmmer_is_installed()) {
        stop("Could not find HMMER. Check if it is installed and in your PATH.")
    }
    if(!is(seq, "AAStringSet")) {
        stop("Argument to 'seq' must be an AAStringSet object.")
    }

    seq_path <- tempfile(pattern = "seq", fileext = ".fasta")
    Biostrings::writeXStringSet(seq, filepath = seq_path)
    return(seq_path)
}


#' Annotate sequences in a FASTA file with PlantTFDB domains
#'
#' @param seq_path Path to FASTA file.
#' @inheritParams annotate_domains
#'
#' @return A data frame as in the element \strong{PlantTFDB} of the output
#' of \code{annotate_domains()}.
#' @noRd
.annotate_planttfdb <- function(seq_path, evalue = 1e-05, threads = 1) {

    hits <- .run_hmmsearch(
        seq_path, .scheme_hmm_path("PlantTFDB"), threads = threads
    )
    hits <- .filter_hmmer(hits, evalue)

    return(.format_hits(hits))
}


#' Annotate sequences in a FASTA file with TAPscan domains
#'
#' @param seq_path Path to FASTA file.
#' @inheritParams annotate_domains
#'
#' @return A data frame as in the element \strong{TAPscan} of the output
#' of \code{annotate_domains()}.
#' @noRd
.annotate_tapscan <- function(seq_path, threads = 1) {

    hits <- .run_hmmsearch(
        seq_path, .scheme_hmm_path("TAPscan"), extra_args = "--cut_ga",
        threads = threads
    )

    return(.format_hits(hits))
}


#' Annotate sequences in a FASTA file with CIS-BP domains
#'
#' DNA-binding domains are kept if they have a per-sequence E-value
#' < 0.01 and a per-domain conditional E-value < 0.01, as in CIS-BP.
#' Then, overlapping hits of different domains are resolved by keeping only
#' the most significant hit, as in CIS-BP (build 2.00 onwards).
#'
#' @param seq_path Path to FASTA file.
#' @inheritParams annotate_domains
#'
#' @return A data frame as in the element \strong{CISBP} of the output
#' of \code{annotate_domains()}.
#' @noRd
.annotate_cisbp <- function(seq_path, threads = 1) {

    hits <- .run_hmmsearch(
        seq_path, .scheme_hmm_path("CISBP"), threads = threads
    )
    hits <- hits[hits$full_evalue < 0.01 & hits$c_evalue < 0.01, ]
    hits <- .resolve_overlaps(hits)

    return(.format_hits(hits))
}


#' Resolve overlapping domain hits by keeping the most significant hits
#'
#' For each sequence, domain hits are processed from the most to the least
#' significant (i.e., by increasing independent E-value), and a hit is
#' removed if its alignment overlaps the alignment of a more significant
#' hit of a different domain that was kept. Overlapping hits of the same
#' domain (e.g., in proteins with domain repeats) are kept.
#'
#' @param hits A data frame with hmmsearch hits as returned
#' by \code{.read_hmmsearch()}.
#'
#' @return The input data frame without overlapping hits, with rows in the
#' original order.
#' @noRd
.resolve_overlaps <- function(hits = NULL) {

    if(nrow(hits) < 2) { return(hits) }

    # Row indices per sequence, from the most to the least significant hit
    ord <- order(hits$i_evalue, -hits$dom_score, method = "radix")
    idx <- split(ord, factor(hits$target_name[ord]))

    keep <- unlist(lapply(idx, function(i) {
        Reduce(function(kept, j) {
            overlaps <- hits$query_name[kept] != hits$query_name[j] &
                hits$ali_from[j] <= hits$ali_to[kept] &
                hits$ali_to[j] >= hits$ali_from[kept]
            if(any(overlaps)) { return(kept) }
            return(c(kept, j))
        }, i[-1], i[1])
    }), use.names = FALSE)

    hits <- hits[sort(keep), ]
    rownames(hits) <- NULL
    return(hits)
}


#' Format hmmsearch hits as domain annotation
#'
#' @param hits A data frame with hmmsearch hits as returned
#' by \code{.read_hmmsearch()}.
#'
#' @return A data frame with variables \strong{Gene}, \strong{Domain},
#' \strong{qlen}, \strong{c_evalue}, \strong{hmm_from}, and
#' \strong{hmm_to}, as in the output of \code{annotate_domains()}.
#' @noRd
.format_hits <- function(hits = NULL) {

    annotation <- data.frame(
        Gene = hits$target_name,
        Domain = hits$query_name,
        qlen = hits$qlen,
        c_evalue = hits$c_evalue,
        hmm_from = hits$hmm_from,
        hmm_to = hits$hmm_to
    )
    return(annotation)
}
