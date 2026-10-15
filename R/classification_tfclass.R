#' Add TFClass classifications based on CIS-BP classifications
#'
#' TFClass superclasses, classes, and families are inferred from CIS-BP
#' families using a correspondence table between the two classification
#' systems (see \code{tfclass_scheme}). CIS-BP families with multiple types
#' of DNA-binding domains (DBDs) that are not in the correspondence table
#' are classified based on the TFClass classifications of each DBD's family
#' (see \code{.resolve_tfclass()}).
#'
#' @param families A data frame with variables \strong{Gene} and
#' \strong{CISBP_family}.
#'
#' @return The input data frame with three additional variables:
#' \strong{TFClass_superclass}, \strong{TFClass_class}, and
#' \strong{TFClass_family}.
#' @noRd
.classify_tfclass <- function(families = NULL) {

    cisbp <- families$CISBP_family
    ucisbp <- unique(cisbp)
    res <- vapply(ucisbp, .resolve_tfclass, character(3), USE.NAMES = FALSE)
    res <- matrix(res, ncol = 3, byrow = TRUE)
    res <- res[match(cisbp, ucisbp), , drop = FALSE]

    families$TFClass_superclass <- res[, 1]
    families$TFClass_class <- res[, 2]
    families$TFClass_family <- res[, 3]
    return(families)
}


#' Get TFClass classification of a CIS-BP family
#'
#' CIS-BP families in the correspondence table between CIS-BP and TFClass
#' are classified as in the table. Other CIS-BP families with multiple
#' types of DBDs (e.g., 'bHLH,Sox') are classified based on the TFClass
#' classes of the families of each DBD (ignoring DBDs without TFClass
#' classes). If all of them are in the same TFClass class, TFs are
#' assigned to this class, and to a TFClass family if all of them are
#' in the same TFClass family. If they are in different TFClass classes,
#' all possible superclasses and classes are reported, separated by ';',
#' and the TFClass family is NA.
#'
#' @param cisbp_family Character scalar with a CIS-BP family, or NA.
#'
#' @return A character vector of length 3 with TFClass superclass, class,
#' and family (NA if not available).
#' @noRd
.resolve_tfclass <- function(cisbp_family = NA) {

    map <- tfclass_map
    res <- rep(NA_character_, 3)
    if(is.na(cisbp_family)) { return(res) }

    cand <- map[map$cisbp_family == cisbp_family, ]
    if(nrow(cand) == 0) {
        dbd_families <- strsplit(cisbp_family, ",", fixed = TRUE)[[1]]
        cand <- map[map$cisbp_family %in% dbd_families, ]
    }
    if(nrow(cand) == 0) { return(res) }

    cand <- cand[order(numeric_version(cand$tfclass_id)), ]
    res[1] <- paste(unique(cand$superclass), collapse = ";")
    res[2] <- paste(unique(cand$class), collapse = ";")
    if(length(unique(cand$class)) == 1 && length(unique(cand$family)) == 1) {
        res[3] <- cand$family[1]
    }
    return(res)
}
