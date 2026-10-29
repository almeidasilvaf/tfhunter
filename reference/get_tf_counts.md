# Get TF frequencies for each species as a `SummarizedExperiment` object

This function counts TFs per family (or other classification level) in
each species from TF classifications obtained with
[`classify_tfs()`](classify_tfs.md), and returns the counts as a
`SummarizedExperiment` object.

## Usage

``` r
get_tf_counts(
  families,
  species_metadata = NULL,
  scheme = c("CISBP", "PlantTFDB", "TAPscan", "PlantTFClass", "TFClass"),
  level = "family"
)
```

## Arguments

- families:

  A named list of data frames with TF classifications for each species,
  as returned by [`classify_tfs()`](classify_tfs.md). List names
  represent species names.

- species_metadata:

  (Optional) A data frame containing species names in row names (names
  must match element names in the **families** list), and species
  metadata (e.g., taxonomic information, ecological information) in
  columns. If NULL, the colData of the `SummarizedExperiment` object
  will be empty.

- scheme:

  Character indicating which classification scheme to use. One of
  'CISBP', 'PlantTFDB', 'TAPscan', 'PlantTFClass', or 'TFClass'.
  Default: 'CISBP'.

- level:

  Character indicating the classification level at which TFs will be
  counted. Possible levels depend on the classification scheme: 'family'
  or 'subfamily' for PlantTFDB and TAPscan, 'superclass', 'class', or
  'family' for Plant-TFClass and TFClass, and 'family' for CIS-BP (see
  variable **Levels** in [`list_schemes()`](list_schemes.md)). Default:
  'family'.

## Value

A SummarizedExperiment object containing TF frequencies per family (or
other level) in each species, as well as species metadata (if
**species_metadata** is not NULL). The rowData of the object contains
the more general classification levels of each row (if any), which can
be used to aggregate counts, as well as the TAP class of each row for
TAPscan. For PlantTFDB, genes assigned to more than one family are
counted once for each family. For Plant-TFClass and TFClass, genes with
ambiguous classifications (i.e., multiple possible classifications
separated by ';') are counted in a separate row with all possible
classifications. Genes with NA at the chosen level (e.g., TFs not
classified by TFClass) are not counted.

## Examples

``` r
data(gsu_annotation)
gsu_families <- classify_tfs(gsu_annotation)

# Simulate TF classifications for 4 species by sampling 100 genes
set.seed(123)
families <- list(
    Gsu1 = gsu_families[sample(nrow(gsu_families), 100), ],
    Gsu2 = gsu_families[sample(nrow(gsu_families), 100), ],
    Gsu3 = gsu_families[sample(nrow(gsu_families), 100), ],
    Gsu4 = gsu_families[sample(nrow(gsu_families), 100), ]
)

# Create species metadata
species_metadata <- data.frame(
    row.names = names(families),
    Division = "Rhodophyta",
    Origin = c("US", "Belgium", "China", "Brazil")
)

# Count TFs per CIS-BP family (default)
se <- get_tf_counts(families, species_metadata)
se
#> class: SummarizedExperiment 
#> dim: 20 4 
#> metadata(0):
#> assays(1): counts
#> rownames(20): ARID/BRIGHT bHLH ... TCR/CxC Zinc cluster
#> rowData names(0):
#> colnames(4): Gsu1 Gsu2 Gsu3 Gsu4
#> colData names(2): Division Origin

# Count TAPscan TFs (excluding TRs and PTs) per subfamily
tapscan_tfs <- lapply(families, function(x) x[x$TAP_class %in% "TF", ])
se_tapscan <- get_tf_counts(
    tapscan_tfs, scheme = "TAPscan", level = "subfamily"
)
se_tapscan
#> class: SummarizedExperiment 
#> dim: 23 4 
#> metadata(0):
#> assays(1): counts
#> rownames(23): ARID bHLH ... Sir2 Zn_clus
#> rowData names(2): TAPscan_family TAP_class
#> colnames(4): Gsu1 Gsu2 Gsu3 Gsu4
#> colData names(0):

# Count TFs per Plant-TFClass superclass
se_superclass <- get_tf_counts(
    families, scheme = "PlantTFClass", level = "superclass"
)

# Count TFs per PlantTFDB family
se_ptfdb <- get_tf_counts(families, scheme = "PlantTFDB")

# Count TFs per TFClass class
se_tfclass <- get_tf_counts(families, scheme = "TFClass", level = "class")
```
