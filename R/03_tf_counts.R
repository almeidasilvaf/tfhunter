
#' Get TF frequencies for each species as a `SummarizedExperiment` object
#'
#' This function counts TFs per family (or other classification level) in
#' each species from TF classifications obtained with \code{classify_tfs()},
#' and returns the counts as a `SummarizedExperiment` object.
#'
#' @param families A named list of data frames with TF classifications for
#' each species, as returned by \code{classify_tfs()}. List names
#' represent species names. 
#' @param species_metadata (Optional) A data frame containing species names
#' in row names (names must match element names in
#' the \strong{families} list), and species metadata
#' (e.g., taxonomic information, ecological information)
#' in columns. If NULL, the colData of the `SummarizedExperiment` object
#' will be empty.
#' @param scheme Character indicating which classification scheme to use.
#' One of 'CISBP', 'PlantTFDB', 'TAPscan', 'PlantTFClass', or 'TFClass'.
#' Default: 'CISBP'.
#' @param level Character indicating the classification level at which TFs
#' will be counted. Possible levels depend on the classification scheme:
#' 'family' or 'subfamily' for PlantTFDB and TAPscan, 'superclass',
#' 'class', or 'family' for Plant-TFClass and TFClass, and 'family' for
#' CIS-BP
#' (see variable \strong{Levels} in \code{list_schemes()}).
#' Default: 'family'.
#'
#' @return A SummarizedExperiment object containing TF frequencies per
#' family (or other level) in each species, as well as species metadata (if
#' \strong{species_metadata} is not NULL). The rowData of the object
#' contains the more general classification levels of each row (if any),
#' which can be used to aggregate counts, as well as the TAP class of each
#' row for TAPscan. For PlantTFDB, genes assigned to more than one family
#' are counted once for each family. For Plant-TFClass and TFClass, genes
#' with ambiguous classifications (i.e., multiple possible classifications
#' separated by ';') are counted in a separate row with all possible
#' classifications. Genes with NA at the chosen level (e.g., TFs
#' not classified into a TFClass family) are not counted.
#'
#' @export
#' @rdname get_tf_counts
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @examples
#' data(gsu_annotation)
#' gsu_families <- classify_tfs(gsu_annotation)
#'
#' # Simulate TF classifications for 4 species by sampling 100 genes
#' set.seed(123)
#' families <- list(
#'     Gsu1 = gsu_families[sample(nrow(gsu_families), 100), ],
#'     Gsu2 = gsu_families[sample(nrow(gsu_families), 100), ],
#'     Gsu3 = gsu_families[sample(nrow(gsu_families), 100), ],
#'     Gsu4 = gsu_families[sample(nrow(gsu_families), 100), ]
#' )
#'
#' # Create species metadata
#' species_metadata <- data.frame(
#'     row.names = names(families),
#'     Division = "Rhodophyta",
#'     Origin = c("US", "Belgium", "China", "Brazil")
#' )
#'
#' # Count TFs per CIS-BP family (default)
#' se <- get_tf_counts(families, species_metadata)
#' se
#'
#' # Count TAPscan TFs (excluding TRs and PTs) per subfamily
#' tapscan_tfs <- lapply(families, function(x) x[x$TAP_class %in% "TF", ])
#' se_tapscan <- get_tf_counts(
#'     tapscan_tfs, scheme = "TAPscan", level = "subfamily"
#' )
#' se_tapscan
#'
#' # Count TFs per Plant-TFClass superclass
#' se_superclass <- get_tf_counts(
#'     families, scheme = "PlantTFClass", level = "superclass"
#' )
#'
#' # Count TFs per PlantTFDB family
#' se_ptfdb <- get_tf_counts(families, scheme = "PlantTFDB")
#'
#' # Count TFs per TFClass class
#' se_tfclass <- get_tf_counts(families, scheme = "TFClass", level = "class")
get_tf_counts <- function(
        families, species_metadata = NULL,
        scheme = c("CISBP", "PlantTFDB", "TAPscan", "PlantTFClass", "TFClass"),
        level = "family"
) {

    scheme <- match.arg(scheme)
    .check_level(scheme, level)
    coldata <- .check_tf_counts_input(families, species_metadata, scheme)

    # Get a long data frame of TF classifications in each species
    long <- lapply(names(families), function(species) {
        fam <- .families_to_count(families[[species]], scheme, level)
        fam$Species <- rep(species, nrow(fam))
        return(fam)
    })
    long <- do.call(rbind, long)

    # Create a feature x species count matrix
    features <- unique(long$Feature)
    features <- features[order(tolower(features), features, method = "radix")]
    counts <- unclass(table(
        factor(long$Feature, levels = features),
        factor(long$Species, levels = names(families))
    ))
    counts <- matrix(
        counts, nrow = nrow(counts), ncol = ncol(counts),
        dimnames = unname(dimnames(counts))
    )

    # Create `SummarizedExperiment` object
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = counts), colData = coldata,
        rowData = .count_rowdata(long, rownames(counts))
    )
    return(se)
}


#' Check if a classification level is valid for a classification scheme
#'
#' @param scheme Character indicating the classification scheme.
#' @param level Character indicating the classification level.
#'
#' @return TRUE (invisibly) if the level is valid. Otherwise, an error is
#' thrown.
#' @noRd
.check_level <- function(scheme, level) {

    valid <- names(.scheme_levels(scheme))
    if(!is.character(level) || length(level) != 1 || !level %in% valid) {
        stop(
            "Invalid level for scheme '", scheme, "'. Valid levels are: ",
            paste0("'", valid, "'", collapse = ", "), "."
        )
    }
    return(invisible(TRUE))
}


#' Check input to \code{get_tf_counts()} and create colData
#'
#' @inheritParams get_tf_counts
#'
#' @return A data frame to be used as colData, with species in rows.
#' @importFrom methods is
#' @noRd
.check_tf_counts_input <- function(
        families, species_metadata = NULL, scheme = "CISBP"
) {

    # Check 1: is `families` a named list of data frames?
    valid <- is(families, "list") && length(families) > 0 &&
        !is.null(names(families)) &&
        all(vapply(families, is, logical(1), "data.frame"))
    if(!valid) {
        stop(
            "Input to 'families' must be a named list of data frames ",
            "as returned by classify_tfs()."
        )
    }

    # Check 2: do data frames have the variables of the chosen scheme?
    cols <- c("Gene", .scheme_columns(scheme))
    has_cols <- vapply(
        families, function(x) all(cols %in% names(x)), logical(1)
    )
    if(!all(has_cols)) {
        stop(
            "Data frames in 'families' must have variables ",
            paste0("'", cols, "'", collapse = ", "), " for scheme '",
            scheme, "'. Missing in: ",
            paste(names(families)[!has_cols], collapse = ", ")
        )
    }

    # Check 3: handle colData depending on `species_metadata`'s class
    coldata <- data.frame(row.names = names(families))
    if(is(species_metadata, "data.frame")) {
        if(!setequal(rownames(species_metadata), names(families))) {
            stop(
                "Names of 'families' must match rownames of ",
                "'species_metadata'."
            )
        }
        coldata <- species_metadata[names(families), , drop = FALSE]
    }
    return(coldata)
}


#' Extract classifications to count from the output of \code{classify_tfs()}
#'
#' @param families A data frame as returned by \code{classify_tfs()}.
#' @param scheme Character indicating the classification scheme.
#' @param level Character indicating the classification level.
#'
#' @return A data frame with variable \strong{Feature} (classification
#' to count, one row per gene and classification), and variables with
#' more general classification levels (and TAP class, for TAPscan).
#' For PlantTFDB, genes assigned to multiple families are included once for
#' each family.
#' @noRd
.families_to_count <- function(
        families, scheme = "CISBP", level = "family"
) {

    levels <- .scheme_levels(scheme)
    pos <- match(level, names(levels))
    parents <- unname(levels[seq_len(pos - 1)])
    if(scheme == "TAPscan") { parents <- c(parents, "TAP_class") }

    families <- families[!is.na(families[[levels[pos]]]), ]
    fam <- c(
        list(Feature = families[[levels[pos]]]),
        as.list(families[parents])
    )
    fam <- lapply(fam, as.character)

    if(scheme == "PlantTFDB") {
        fam <- lapply(fam, function(x) {
            as.character(unlist(strsplit(x, ";")))
        })
    }
    fam <- as.data.frame(fam, check.names = FALSE)
    return(fam)
}


#' Create rowData for the output of \code{get_tf_counts()}
#'
#' @param long A data frame as returned by \code{.families_to_count()}
#' for all species, with an additional variable \strong{Species}.
#' @param features Character vector of features (rows of the count matrix).
#'
#' @return A data frame with more general classification levels for each
#' feature, with features in row names. If a feature has multiple values for
#' a level, values are separated by ';'.
#' @noRd
.count_rowdata <- function(long, features) {

    parents <- setdiff(names(long), c("Feature", "Species"))
    feature <- factor(long$Feature, levels = features)
    parent_values <- lapply(parents, function(p) {
        vapply(split(long[[p]], feature), function(x) {
            x <- unique(x[!is.na(x)])
            if(length(x) == 0) { return(NA_character_) }
            return(paste(x, collapse = ";"))
        }, character(1), USE.NAMES = FALSE)
    })
    names(parent_values) <- parents

    rowdata <- data.frame(row.names = features)
    rowdata[parents] <- parent_values
    return(rowdata)
}
