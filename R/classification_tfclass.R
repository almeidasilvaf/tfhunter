#' Add TFClass classifications based on CIS-BP classifications
#'
#' TFClass superclasses, classes, and families are inferred from CIS-BP
#' families using a correspondence table between the two classification
#' systems (see \code{tfclass_scheme}). CIS-BP families with multiple types
#' of DNA-binding domains (DBDs) that are not in the correspondence table
#' are classified based on the TFClass classifications of each DBD's family
#' (see \code{.resolve_tfclass()}). If TFClass evidence is available, genes
#' with multiple possible TFClass families are assigned to a single family
#' when possible (see \code{.best_tfclass_family()}).
#'
#' @param families A data frame with variables \strong{Gene} and
#' \strong{CISBP_family}.
#' @param evidence (Optional) A data frame with TFClass evidence, as in the
#' element \strong{TFClass} of the output of \code{annotate_domains()}.
#' If NULL, genes with multiple possible TFClass families are not assigned
#' to a single family. Default: NULL.
#' @param min_shared Numeric indicating the minimum fraction of shared
#' 5-mers to assign C2H2 ZF TFs to a single TFClass family. Default: 0.1.
#'
#' @return The input data frame with three additional variables:
#' \strong{TFClass_superclass}, \strong{TFClass_class}, and
#' \strong{TFClass_family}.
#' @noRd
.classify_tfclass <- function(
        families = NULL, evidence = NULL, min_shared = 0.1
) {

    cisbp <- families$CISBP_family
    ucisbp <- unique(cisbp)
    res <- vapply(ucisbp, .resolve_tfclass, character(4), USE.NAMES = FALSE)
    res <- matrix(res, ncol = 4, byrow = TRUE)
    res <- res[match(cisbp, ucisbp), , drop = FALSE]

    if(!is.null(evidence)) {
        best <- .best_tfclass_family(
            families$Gene, cisbp, res[, 4], evidence, min_shared
        )
        tfclass_fams <- .tfclass_families()
        idx <- match(best, tfclass_fams$id)
        found <- which(!is.na(idx))
        res[found, ] <- as.matrix(tfclass_fams[idx[found], ])
    }

    families$TFClass_superclass <- res[, 1]
    families$TFClass_class <- res[, 2]
    families$TFClass_family <- res[, 3]
    return(families)
}


#' Get TFClass classification of a CIS-BP family
#'
#' CIS-BP families in the correspondence table between CIS-BP and TFClass
#' are classified as in the table, where CIS-BP families corresponding to
#' multiple TFClass families have all possible families separated by ';'.
#' Other CIS-BP families with multiple types of DBDs (e.g., 'bHLH,Sox')
#' are classified based on the TFClass classifications of the families of
#' each DBD (ignoring DBDs without TFClass classes), and all possible
#' superclasses, classes, and families are reported, separated by ';'.
#' The TFClass family is NA if the family of any of the DBDs is NA.
#'
#' @param cisbp_family Character scalar with a CIS-BP family, or NA.
#'
#' @return A character vector of length 4 with TFClass superclass, class,
#' family, and family ID (NA if not available).
#' @noRd
.resolve_tfclass <- function(cisbp_family = NA) {

    map <- tfclass_map
    res <- rep(NA_character_, 4)
    if(is.na(cisbp_family)) { return(res) }

    cand <- map[map$cisbp_family == cisbp_family, ]
    if(nrow(cand) == 0) {
        dbd_families <- strsplit(cisbp_family, ",", fixed = TRUE)[[1]]
        cand <- map[map$cisbp_family %in% dbd_families, ]
    }
    if(nrow(cand) == 0) { return(res) }

    cand <- cand[order(numeric_version(sub(";.*", "", cand$tfclass_id))), ]
    res[1] <- paste(unique(cand$superclass), collapse = ";")
    res[2] <- paste(unique(cand$class), collapse = ";")
    if(!anyNA(cand$family)) {
        fam <- unlist(strsplit(cand$family, ";", fixed = TRUE))
        ids <- unlist(strsplit(cand$tfclass_id, ";", fixed = TRUE))
        res[3] <- paste(unique(fam), collapse = ";")
        res[4] <- paste(unique(ids), collapse = ";")
    }
    return(res)
}


#' Get TFClass superclass, class, and family names of TFClass family IDs
#'
#' @return A data frame with variables \strong{superclass}, \strong{class},
#' \strong{family}, and \strong{id} (TFClass family ID), with one row per
#' TFClass family in the correspondence table between CIS-BP and TFClass.
#' @noRd
.tfclass_families <- function() {

    map <- tfclass_map[
        !is.na(tfclass_map$family) & !grepl(";", tfclass_map$class),
    ]
    ids <- strsplit(map$tfclass_id, ";", fixed = TRUE)
    n <- lengths(ids)
    fams <- data.frame(
        superclass = rep(map$superclass, n), class = rep(map$class, n),
        family = unlist(strsplit(map$family, ";", fixed = TRUE)),
        id = unlist(ids)
    )
    fams <- fams[!duplicated(fams$id), ]
    rownames(fams) <- NULL
    return(fams)
}


#' Get IDs of TFClass families with profile HMMs
#'
#' Profile HMMs are available for all TFClass families of CIS-BP families
#' with a single type of DBD that correspond to multiple TFClass families,
#' except for C2H2 ZF (see \emph{inst/extdata/TFClass.hmm.gz}).
#'
#' @return A character vector with TFClass family IDs.
#' @noRd
.tfclass_hmm_ids <- function() {

    keep <- grepl(";", tfclass_map$tfclass_id) &
        tfclass_map$cisbp_family %in% cisbp_map$family &
        tfclass_map$cisbp_family != "C2H2 ZF"
    ids <- unlist(strsplit(tfclass_map$tfclass_id[keep], ";", fixed = TRUE))
    return(ids)
}


#' Choose a single TFClass family for genes with multiple possible families
#'
#' For genes in CIS-BP family C2H2 ZF, the TFClass family is the family of
#' the reference C2H2 ZF TF that shares the largest fraction of the
#' gene's 5-mers, if this fraction is >= \strong{min_shared}. For other
#' genes, if profile HMMs are available for all possible TFClass families,
#' the TFClass family is the family whose profile HMM has the highest
#' score, if this score is >= \strong{min_gap} bits higher than the score
#' of any other possible family.
#'
#' @param genes Character vector of gene IDs.
#' @param cisbp Character vector of CIS-BP families, in the same order
#' as \strong{genes}.
#' @param ids Character vector of TFClass family IDs (separated by ';' for
#' genes with multiple possible families), in the same order
#' as \strong{genes}.
#' @param evidence A data frame with TFClass evidence, as in the element
#' \strong{TFClass} of the output of \code{annotate_domains()}.
#' @param min_shared Numeric indicating the minimum fraction of shared
#' 5-mers for C2H2 ZF TFs. Default: 0.1.
#' @param min_gap Numeric indicating the minimum difference in profile HMM
#' scores (in bits) between the best and the second best TFClass families.
#' Default: 10.
#'
#' @return A character vector with the TFClass family ID chosen for each
#' gene, or NA if no family could be chosen (or if genes have a
#' single possible family).
#' @noRd
.best_tfclass_family <- function(
        genes, cisbp, ids, evidence, min_shared = 0.1, min_gap = 10
) {

    best <- rep(NA_character_, length(genes))
    amb <- which(grepl(";", ids) & genes %in% evidence$Gene)
    if(length(amb) == 0) { return(best) }

    hmm_ids <- .tfclass_hmm_ids()
    rows <- split(seq_len(nrow(evidence)), evidence$Gene)
    best[amb] <- vapply(amb, function(i) {
        cand <- strsplit(ids[i], ";", fixed = TRUE)[[1]]
        ev <- evidence[rows[[genes[i]]], ]
        ev <- ev[ev$Domain %in% cand, ]

        kmer <- cisbp[i] == "C2H2 ZF"
        if(!kmer && !all(cand %in% hmm_ids)) { return(NA_character_) }
        ev <- ev[ev$method == ifelse(kmer, "k-mer", "HMM"), ]
        if(nrow(ev) == 0) { return(NA_character_) }

        ev <- ev[order(-ev$score), ]
        pass <- ifelse(
            kmer, ev$score[1] >= min_shared,
            ev$score[1] - c(ev$score[-1], 0)[1] >= min_gap
        )
        return(ifelse(pass, ev$Domain[1], NA_character_))
    }, character(1))

    return(best)
}


#' Get evidence to choose among multiple TFClass families
#'
#' @param seq An `AAStringSet` object.
#' @param cisbp A data frame with CIS-BP domain annotation, as in the element
#' \strong{CISBP} of the output of \code{annotate_domains()}.
#' @inheritParams annotate_domains
#'
#' @return A data frame as in the element \strong{TFClass} of the output
#' of \code{annotate_domains()}.
#' @noRd
.annotate_tfclass <- function(seq, cisbp, threads = 1) {

    fam <- .classify_cisbp(cisbp)
    ids <- vapply(fam$Family, function(f) {
        .resolve_tfclass(f)[4]
    }, character(1), USE.NAMES = FALSE)
    amb <- grepl(";", ids)
    has_hmm <- vapply(strsplit(ids, ";", fixed = TRUE), function(x) {
        all(x %in% .tfclass_hmm_ids())
    }, logical(1))

    # Sequence IDs as in HMMER output (i.e., up to the first space)
    seq_ids <- sub(" .*", "", names(seq))
    c2h2 <- fam$Gene[amb & fam$Family == "C2H2 ZF"]
    hmm <- fam$Gene[amb & fam$Family != "C2H2 ZF" & has_hmm]

    evidence <- rbind(
        .tfclass_hmm_evidence(seq[match(hmm, seq_ids)], threads),
        .tfclass_kmer_evidence(seq[match(c2h2, seq_ids)])
    )
    rownames(evidence) <- NULL
    return(evidence)
}


#' Score sequences with profile HMMs of TFClass families
#'
#' @param seq An `AAStringSet` object.
#' @inheritParams annotate_domains
#'
#' @return A data frame as in the element \strong{TFClass} of the output
#' of \code{annotate_domains()}, with one row per gene and TFClass family
#' with a hit (full sequence E-value < 0.01). Variables \strong{qlen},
#' \strong{c_evalue}, \strong{hmm_from}, and \strong{hmm_to} refer to the
#' best domain hit.
#' @importFrom Biostrings writeXStringSet
#' @noRd
.tfclass_hmm_evidence <- function(seq, threads = 1) {

    if(length(seq) == 0) { return(.empty_tfclass_evidence()) }

    seq_path <- tempfile(pattern = "seq", fileext = ".fasta")
    on.exit(unlink(seq_path))
    Biostrings::writeXStringSet(seq, filepath = seq_path)
    hits <- .run_hmmsearch(
        seq_path, .scheme_hmm_path("TFClass"), threads = threads
    )
    hits <- hits[hits$full_evalue < 0.01, ]

    # Keep the best domain hit for each gene and family
    hits <- hits[order(
        hits$target_name, hits$query_name, -hits$dom_score, method = "radix"
    ), ]
    hits <- hits[!duplicated(hits[, c("target_name", "query_name")]), ]

    evidence <- .format_hits(hits)
    evidence$score <- hits$full_score
    evidence$method <- rep("HMM", nrow(hits))
    return(evidence)
}


#' Compare C2H2 ZF TFs to reference C2H2 ZF TFs based on shared k-mers
#'
#' @param seq An `AAStringSet` object.
#'
#' @return A data frame as in the element \strong{TFClass} of the output
#' of \code{annotate_domains()}, with one row per gene and TFClass family
#' with shared 5-mers.
#' @noRd
.tfclass_kmer_evidence <- function(seq) {

    if(length(seq) == 0) { return(.empty_tfclass_evidence()) }

    ref <- .tfclass_c2h2_refs()
    ref_family <- sub("\\|.*", "", names(ref))
    shared <- .shared_kmers(seq, ref)

    # Highest fraction of shared k-mers per family
    families <- sort(unique(ref_family))
    score <- vapply(families, function(f) {
        apply(shared$counts[, ref_family == f, drop = FALSE], 1, max)
    }, numeric(length(seq)))
    score <- matrix(score, nrow = length(seq)) / pmax(1, shared$n)

    evidence <- data.frame(
        Gene = rep(sub(" .*", "", names(seq)), length(families)),
        Domain = rep(families, each = length(seq)),
        qlen = NA_real_, c_evalue = NA_real_,
        hmm_from = NA_real_, hmm_to = NA_real_,
        score = as.vector(score), method = "k-mer"
    )
    evidence <- evidence[evidence$score > 0, ]
    evidence <- evidence[order(evidence$Gene, evidence$Domain), ]
    return(evidence)
}


#' Create an empty data frame of TFClass evidence
#'
#' @return A data frame with 0 rows and the variables of the element
#' \strong{TFClass} of the output of \code{annotate_domains()}.
#' @noRd
.empty_tfclass_evidence <- function() {

    evidence <- data.frame(
        Gene = character(0), Domain = character(0), qlen = numeric(0),
        c_evalue = numeric(0), hmm_from = numeric(0), hmm_to = numeric(0),
        score = numeric(0), method = character(0)
    )
    return(evidence)
}


#' Read reference C2H2 ZF TFs
#'
#' Reference sequences are human and mouse C2H2 ZF TFs with known TFClass
#' families, named '<TFClass family ID>|<Ensembl gene ID>'. Sequences are
#' read once per R session.
#'
#' @return An `AAStringSet` object.
#' @importFrom Biostrings readAAStringSet
#' @noRd
.tfclass_c2h2_refs <- function() {

    ref <- .pkg_env$TFClass_C2H2
    if(is.null(ref)) {
        ref <- Biostrings::readAAStringSet(system.file(
            "extdata", "TFClass_C2H2.fa.gz", package = "tfhunter"
        ))
        .pkg_env$TFClass_C2H2 <- ref
    }
    return(ref)
}


#' Count k-mers shared between query and reference sequences
#'
#' K-mers that occur in more than 2% of the reference sequences (or in
#' more than 2 sequences, whichever is higher) are ignored, as they are
#' not informative (e.g., k-mers in the linkers between C2H2 zinc fingers).
#'
#' @param query,ref `AAStringSet` objects with query and reference
#' sequences.
#' @param k Numeric indicating the k-mer size. Default: 5.
#' @param max_ref_freq Numeric indicating the maximum fraction of reference
#' sequences that can contain a k-mer. Default: 0.02.
#'
#' @return A list with elements:
#' \describe{
#'   \item{counts}{A matrix with queries in rows and references in columns,
#'   with the number of informative k-mers shared by each pair.}
#'   \item{n}{A numeric vector with the number of informative k-mers in
#'   each query.}
#' }
#' @noRd
.shared_kmers <- function(query, ref, k = 5, max_ref_freq = 0.02) {

    rk <- .kmer_table(ref, k)
    qk <- .kmer_table(query, k)

    # Remove uninformative k-mers
    ucode <- unique(rk$code)
    nseqs <- tabulate(match(rk$code, ucode))
    common <- ucode[nseqs > max(2, max_ref_freq * length(ref))]
    rk <- rk[!rk$code %in% common, ]
    rk <- rk[order(rk$code), ]
    qk <- qk[!qk$code %in% common, ]
    n <- tabulate(qk$seq, nbins = length(query))

    # Pairs of query and reference sequences for each shared k-mer
    from <- match(qk$code, rk$code)
    qk <- qk[!is.na(from), ]
    from <- from[!is.na(from)]
    len <- findInterval(qk$code, rk$code) - from + 1
    ridx <- rk$seq[sequence(len, from = from)]
    qidx <- rep(qk$seq, len)

    counts <- tabulate(
        (ridx - 1) * length(query) + qidx,
        nbins = length(query) * length(ref)
    )
    counts <- matrix(counts, nrow = length(query))
    return(list(counts = counts, n = n))
}


#' Get unique amino acid k-mers of each sequence
#'
#' K-mers are encoded as numbers (base 20), and k-mers with characters
#' other than the 20 standard amino acids are ignored.
#'
#' @param seq An `AAStringSet` object.
#' @param k Numeric indicating the k-mer size. Default: 5.
#'
#' @return A data frame with variables \strong{code} (numeric code of the
#' k-mer) and \strong{seq} (index of the sequence in \strong{seq}), with
#' one row per unique k-mer in each sequence.
#' @importFrom Biostrings width
#' @noRd
.kmer_table <- function(seq, k = 5) {

    lut <- rep(NA_real_, 256)
    lut[as.integer(charToRaw("ACDEFGHIKLMNPQRSTVWY")) + 1] <- seq_len(20) - 1
    x <- lut[as.integer(charToRaw(paste(as.character(seq), collapse = ""))) + 1]
    sid <- rep(seq_along(seq), Biostrings::width(seq))

    n <- length(x) - k + 1
    if(n < 1) { return(data.frame(code = numeric(0), seq = integer(0))) }
    start <- seq_len(n)
    code <- Reduce(function(code, j) {
        code * 20 + x[start + j]
    }, seq_len(k - 1), x[start])

    # Remove k-mers with non-standard amino acids or spanning two sequences
    keep <- !is.na(code) & sid[start] == sid[start + k - 1]
    kmers <- data.frame(code = code[keep], seq = sid[start][keep])
    kmers <- kmers[!duplicated((kmers$seq - 1) * 20^k + kmers$code), ]
    return(kmers)
}
