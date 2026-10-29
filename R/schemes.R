#' List available TF classification schemes
#'
#' @return A data frame with one row per classification scheme and the
#' following variables:
#' \describe{
#'   \item{Scheme}{Name of the classification scheme, to be passed to
#'   the argument \strong{scheme} of \code{annotate_domains()},
#'   \code{classify_tfs()}, and \code{get_tf_counts()}.}
#'   \item{Taxa}{Taxonomic group for which the scheme was designed.}
#'   \item{Levels}{Classification levels, from the most general to
#'   the most specific.}
#'   \item{Annotation}{Domain annotation(s) required by the scheme
#'   (i.e., elements of the output of \code{annotate_domains()}), separated
#'   by ', '.}
#'   \item{Reference}{Reference of the classification scheme.}
#' }
#'
#' @export
#' @rdname list_schemes
#' @examples
#' list_schemes()
list_schemes <- function() {

    schemes <- data.frame(
        Scheme = c(
            "PlantTFDB", "TAPscan", "PlantTFClass", "CISBP", "TFClass"
        ),
        Taxa = c(
            "Plants", "Plants and algae", "Plants", "Eukaryotes", "Mammals"
        ),
        Levels = c(
            "family, subfamily", "family, subfamily",
            "superclass, class, family", "family",
            "superclass, class, family"
        ),
        Annotation = c(
            "PlantTFDB", "TAPscan", "PlantTFDB, TAPscan", "CISBP",
            "CISBP, TFClass"
        ),
        Reference = c(
            "Jin et al. (2017), doi:10.1093/nar/gkw982",
            "Petroll et al. (2025), doi:10.1111/tpj.17184",
            "Blanc-Mathieu et al. (2024), doi:10.1016/j.tplants.2023.06.023",
            "Weirauch et al. (2014), doi:10.1016/j.cell.2014.08.009",
            "Wingender et al. (2018), doi:10.1093/nar/gkx987"
        )
    )
    return(schemes)
}


#' Check and expand classification schemes
#'
#' @param scheme Character vector with classification scheme(s), or 'all'.
#' @param available (Optional) Character vector with the names of the
#' available domain annotations. If not NULL, 'all' is expanded to the schemes
#' whose domain annotations are all available, and an error is thrown if
#' a requested scheme requires a domain annotation that is not available.
#'
#' @return A character vector with the classification schemes, in the
#' same order as in \code{list_schemes()}.
#' @noRd
.check_scheme <- function(scheme = "all", available = NULL) {

    schemes <- list_schemes()$Scheme
    valid <- is.character(scheme) && length(scheme) > 0 &&
        all(scheme %in% c("all", schemes))
    if(!valid) {
        stop(
            "Argument to 'scheme' must be 'all' or one or more of: ",
            paste0("'", schemes, "'", collapse = ", "), "."
        )
    }

    if(is.null(available)) { available <- .annotation_sets(schemes) }
    has_annot <- vapply(schemes, function(s) {
        all(.annotation_sets(s) %in% available)
    }, logical(1))

    if("all" %in% scheme) {
        scheme <- schemes[has_annot]
        if(length(scheme) == 0) {
            stop("No classification scheme can be used with this annotation.")
        }
    }
    missing <- setdiff(scheme, schemes[has_annot])
    if(length(missing) > 0) {
        stop(
            "Missing domain annotation for scheme(s) ",
            paste0("'", missing, "'", collapse = ", "), ". Run ",
            "annotate_domains() with these schemes first."
        )
    }

    scheme <- schemes[schemes %in% scheme]
    return(scheme)
}


#' Get domain annotations required by classification schemes
#'
#' @param scheme Character vector with classification schemes.
#'
#' @return A character vector with the names of the domain annotations
#' (and profile HMM sets) required by the classification schemes.
#' @noRd
.annotation_sets <- function(scheme) {

    schemes <- list_schemes()
    annot <- schemes$Annotation[schemes$Scheme %in% scheme]
    annot <- unique(unlist(strsplit(annot, ", ")))
    return(annot)
}


#' Get variables with classification levels of a classification scheme
#'
#' @param scheme Character indicating the classification scheme. One of
#' the schemes in \code{list_schemes()}.
#'
#' @return A named character vector with the names of the variables in the
#' output of \code{classify_tfs()} for each classification level, from
#' the most general to the most specific level. Vector names are the
#' classification levels.
#' @noRd
.scheme_levels <- function(scheme) {

    levels <- switch(
        scheme,
        PlantTFDB = c(
            family = "PlantTFDB_family", subfamily = "PlantTFDB_subfamily"
        ),
        TAPscan = c(
            family = "TAPscan_family", subfamily = "TAPscan_subfamily"
        ),
        PlantTFClass = c(
            superclass = "PlantTFClass_superclass",
            class = "PlantTFClass_class", family = "PlantTFClass_family"
        ),
        CISBP = c(family = "CISBP_family"),
        TFClass = c(
            superclass = "TFClass_superclass",
            class = "TFClass_class", family = "TFClass_family"
        )
    )
    return(levels)
}


#' Get all variables of a classification scheme
#'
#' @param scheme Character indicating the classification scheme. One of
#' the schemes in \code{list_schemes()}.
#'
#' @return A character vector with the names of the variables of the
#' classification scheme in the output of \code{classify_tfs()}, including
#' classification levels and additional variables (i.e., TAP class for
#' TAPscan).
#' @noRd
.scheme_columns <- function(scheme) {

    cols <- unname(.scheme_levels(scheme))
    if(scheme == "TAPscan") { cols <- c(cols, "TAP_class") }
    return(cols)
}


#' Get path to the uncompressed profile HMMs of a domain annotation set
#'
#' Profile HMMs are stored gzip-compressed in the package (one file per
#' domain annotation set, named \emph{<set>.hmm.gz}). This function
#' decompresses them to a temporary file once per R session and returns the
#' path to the uncompressed file.
#'
#' @param set Character indicating the domain annotation set. One of
#' 'PlantTFDB', 'TAPscan', 'CISBP', or 'TFClass'.
#'
#' @return Character with the path to the uncompressed HMM file.
#' @noRd
.scheme_hmm_path <- function(set) {

    key <- paste0(set, "_hmm")
    path <- .pkg_env[[key]]
    if(is.null(path) || !file.exists(path)) {
        gz <- system.file(
            "extdata", paste0(set, ".hmm.gz"), package = "tfhunter"
        )
        path <- tempfile(pattern = set, fileext = ".hmm")
        con <- gzfile(gz, open = "rt")
        on.exit(close(con))
        writeLines(readLines(con), path)
        .pkg_env[[key]] <- path
    }
    return(path)
}
