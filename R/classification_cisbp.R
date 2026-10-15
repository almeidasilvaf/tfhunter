#' Identify TFs and classify them in families using CIS-BP's scheme
#'
#' As in CIS-BP, each gene is assigned to a family defined by the set of
#' DNA-binding domains (DBDs) it contains, regardless of the number and order
#' of DBDs. Each DBD corresponds to a family name (see \code{cisbp_scheme}),
#' and genes with more than one type of DBD are assigned to a family whose
#' name is the family names of all DBDs sorted alphabetically
#' (case-insensitive) and separated by ',' (e.g., 'Homeodomain,POU').
#'
#' @param domain_annotation A data frame with variables
#' \strong{Gene} and \strong{Domain}, as in the element \strong{CISBP}
#' of the output of \code{annotate_domains()}.
#'
#' @return A 2-column data frame with variables \strong{Gene} and
#' \strong{Family}.
#' @noRd
.classify_cisbp <- function(domain_annotation = NULL) {

    dbd_family <- cisbp_map$family[
        match(domain_annotation$Domain, cisbp_map$pfam)
    ]
    keep <- !is.na(dbd_family)
    if(!any(keep)) {
        return(data.frame(Gene = character(0), Family = character(0)))
    }

    # Unique DBD families per gene, sorted alphabetically (case-insensitive)
    dbds <- unique(data.frame(
        Gene = domain_annotation$Gene[keep], Family = dbd_family[keep]
    ))
    dbds <- dbds[order(dbds$Gene, tolower(dbds$Family), method = "radix"), ]

    genes <- factor(dbds$Gene, levels = unique(dbds$Gene))
    families <- vapply(
        split(dbds$Family, genes), paste, character(1), collapse = ","
    )

    fam_df <- data.frame(Gene = names(families), Family = unname(families))
    return(fam_df)
}
