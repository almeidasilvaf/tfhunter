Data creation
================

# Data in data/

## gsu.rda

This object contains the protein sequences of the algae *Galdieria
sulphuraria*, and it was downloaded from PLAZA Diatoms. Only genes that
have domains of interest (used for classification with any of the
schemes in `list_schemes()`) were kept for package size issues.

``` r
gsu_full <- Biostrings::readAAStringSet(
  "https://ftp.psb.ugent.be/pub/plaza/plaza_diatoms_01/Fasta/proteome.selected_transcript.gsu.fasta.gz"
)
names(gsu_full) <- gsub("\ .*", "", names(gsu_full))

# Annotate domains in the whole proteome (see gsu_annotation.rda below)
gsu_annotation <- annotate_domains(gsu_full)
genes <- unique(unlist(lapply(gsu_annotation, `[[`, "Gene")))
gsu <- gsu_full[names(gsu_full) %in% genes]

outfile <- paste0(tempdir(), "/gsu.fasta")
Biostrings::writeXStringSet(gsu, filepath = outfile)
gsu <- Biostrings::readAAStringSet(outfile)

usethis::use_data(
  gsu, compress = "xz", overwrite = TRUE
)
```

## gsu_annotation.rda

Domains were annotated in the whole proteome of *G. sulphuraria* (object
`gsu_full` created above).

``` r
gsu_annotation <- annotate_domains(gsu_full)

usethis::use_data(
  gsu_annotation, compress = "xz", overwrite = TRUE
)
```

## planttfdb_scheme.rda

The classification scheme is the same as implemented in PlantTFDB, but
three family/subfamily names were changed to match the names used in the
output of `classify_tfs()` (‘M_type’ to ‘M-type’, ‘MYB_related’ to
‘MYB-related’, and ‘LBD (AS2/LOB)’ to ‘LBD’).

The classification scheme is the same as implemented in PlantTFDB.

``` r
library(tidyverse)
table <- rvest::read_html("http://planttfdb.gao-lab.org/help_famschema.php") %>%
  rvest::html_table()
table <- as.data.frame(table[[3]])

table <- table[-1, -1]
names(table) <- c("Family", "Subfamily", "DBD", "Auxiliary", "Forbidden")
table[table == ""] <- NA

table$DBD <- gsub("\\r\\n\\t\\t\\r\\n\\t\\t", " or ", table$DBD)
table$DBD <- gsub("\ \\(self-build\\)", "", table$DBD)
table$Auxiliary <- gsub("\ \\(self-build\\)", "", table$Auxiliary)
table$Auxiliary <- gsub("BELLor", "BELL or", table$Auxiliary)
table$Forbidden <- gsub("\\r\\n\\t\\tor", " or ", table$Forbidden)

planttfdb_scheme <- table
planttfdb_scheme$DBD[3] <- gsub("or", "and", planttfdb_scheme$DBD[3])

# Use the same family and subfamily names as in the output of classify_tfs()
fix_names <- c(
    M_type = "M-type", MYB_related = "MYB-related", `LBD (AS2/LOB)` = "LBD"
)
fix <- function(x) ifelse(x %in% names(fix_names), fix_names[x], x)
planttfdb_scheme$Family <- fix(planttfdb_scheme$Family)
planttfdb_scheme$Subfamily <- fix(planttfdb_scheme$Subfamily)
rownames(planttfdb_scheme) <- NULL

usethis::use_data(planttfdb_scheme,
                  compress = "xz", overwrite = TRUE)
```

## tap_class.rda

Classes of transcription-associated proteins (TAPs) in TAPscan v4 (TF,
TR, or PT), as well as descriptions and references for each TAP family,
were obtained from the GitHub repository of TAPscan’s website
(<https://github.com/Rensing-Lab/TAPscan-v4-website>). Family MADS_MIKC,
which is not included in the original file, was manually added as a TF
family (as family MADS).

``` r
url <- paste0(
    "https://raw.githubusercontent.com/Rensing-Lab/TAPscan-v4-website/",
    "c13c073b1fda474101b95009fe176eb9cb5174d2/_data/import-tapinfo/",
    "tapinfo_v4.csv"
)
lines <- readLines(url, encoding = "UTF-8")[-1]

# Fields are tab-separated; descriptions and references contain quotes, so
# the file is parsed manually instead of with read.delim()
fields <- strsplit(lines, "\t")
field <- function(i) vapply(fields, `[`, character(1), i)
strip_html <- function(x) gsub("<[^>]+>", "", x)
unquote <- function(x) gsub('^"|"$', "", trimws(x))

refs <- vapply(field(3), function(x) {
    x <- strsplit(unquote(x), '","', fixed = TRUE)[[1]]
    paste(trimws(strip_html(x)), collapse = " | ")
}, character(1), USE.NAMES = FALSE)

tap_class <- data.frame(
    TAP_family = trimws(field(1)),
    TAP_class = trimws(field(4)),
    Description = trimws(strip_html(unquote(field(2)))),
    References = refs
)

# Add families missing from the original file
missing_fams <- data.frame(
    TAP_family = "MADS_MIKC", TAP_class = "TF",
    Description = NA_character_, References = NA_character_
)
tap_class <- rbind(tap_class, missing_fams)

usethis::use_data(tap_class, compress = "xz", overwrite = TRUE)
```

## tapscan_scheme.rda

TAPscan’s classification scheme, in a format similar to
`planttfdb_scheme`, with one row per classification rule. This data set
was created from TAPscan’s classification rules (v82) and the internal
data in *sysdata.rda* (see below), so the code below must be run with
`devtools::load_all()`.

``` r
rules <- tapscan_rules
rule_fams <- unique(rules$family)

# Make pseudo-domains used by TAPscan more readable
pseudo <- c(
    `MYB-2R` = "Myb_DNA-binding (2)", `MYB-3R` = "Myb_DNA-binding (3)",
    `MYB-4R` = "Myb_DNA-binding (4)", two_or_more_LIM = "LIM (>=2)"
)
rules$domain <- ifelse(
    rules$domain %in% names(pseudo), pseudo[rules$domain], rules$domain
)
collapse_domains <- function(type, sep) {
    vapply(rule_fams, function(f) {
        d <- rules$domain[rules$family == f & rules$type == type]
        if(length(d) == 0) { return(NA_character_) }
        paste(d, collapse = sep)
    }, character(1), USE.NAMES = FALSE)
}

tapscan_scheme <- data.frame(
    .tapscan_family_map(rule_fams),
    Rule = rule_fams,
    Should = collapse_domains("should", " and "),
    Should_not = collapse_domains("should not", " or ")
)
tapscan_scheme <- .add_tap_class(
    setNames(tapscan_scheme, c("TAPscan_family", "TAPscan_subfamily", 
                               "Rule", "Should", "Should_not"))
)
names(tapscan_scheme) <- c(
    "Family", "Subfamily", "Rule", "Should", "Should_not", "TAP_class"
)
tapscan_scheme <- tapscan_scheme[order(
    tolower(tapscan_scheme$Family), tapscan_scheme$Subfamily, 
    method = "radix"
), ]
rownames(tapscan_scheme) <- NULL

usethis::use_data(tapscan_scheme, compress = "xz", overwrite = TRUE)
```

## cisbp_scheme.rda

DNA-binding domains (DBDs) used by CIS-BP and their corresponding TF
families were obtained from the file *TF_Information.txt* of CIS-BP
(build 3.10). In CIS-BP, each TF is assigned to a family defined by the
set of DBDs it contains (variable *DBDs*), and families of TFs with
multiple types of DBDs are named after the families of each DBD, sorted
alphabetically (case-insensitive) and separated by ‘,’. This rule
reproduces the family of all 392,337 TFs in CIS-BP 3.10, so only TFs
with a single type of DBD were used to get the family of each DBD.

DBD names were matched to accessions in Pfam release 38.2. Seven Pfam
families were renamed since they were used by CIS-BP, and they were
matched manually. The Pfam family FLO_LFY was split into SAM_LFY
(PF01698) and C_LFY_FLO (PF17538, the C-terminal DBD), so both were used
for family ‘Leafy’.

The RNA-binding domains RRM_1 and zf-CCCH were excluded. CIS-BP
classifies TFs in families ‘RRM’ (1 TF) and ‘CCCH ZF’ (5 TFs) only if
there is experimental evidence of DNA binding, and these domains are
never used for TFs with other DBDs. Including them would identify
thousands of RNA-binding proteins as TFs.

``` r
# Download TF information from CIS-BP
url <- paste0(
    "https://cisbp.ccbr.utoronto.ca/data/3_10/DataFiles/Bulk_downloads/",
    "EntireDataset/TF_Information.zip"
)
zip_path <- file.path(tempdir(), "TF_Information.zip")
options(timeout = 3600)
download.file(url, destfile = zip_path, mode = "wb")
inner <- unzip(zip_path, "TF_Information.txt.zip", exdir = tempdir())
txt <- unzip(inner, "TF_Information.txt", exdir = tempdir())

tfinfo <- read.delim(txt, quote = "", na.strings = ".")
tfinfo <- tfinfo[!duplicated(tfinfo$TF_ID), ]

# Check that families are defined by sets of DBDs
dbd_sets <- vapply(strsplit(tfinfo$DBDs, ","), function(x) {
    paste(sort(x), collapse = ",")
}, character(1))
stopifnot(all(tapply(tfinfo$Family_Name, dbd_sets, function(x) {
    length(unique(x)) == 1
})))

single <- unique(tfinfo[
    !grepl(",", tfinfo$DBDs) & tfinfo$DBDs != "UNKNOWN",
    c("Family_Name", "DBDs")
])

# RNA-binding domains used by CIS-BP only for TFs with experimental evidence
# of DNA binding (never in TFs with multiple DBD types)
single <- single[!single$DBDs %in% c("RRM_1", "zf-CCCH"), ]

# Match DBD names to Pfam accessions
clans_gz <- file.path(tempdir(), "Pfam-A.clans.tsv.gz")
download.file(
    paste0(
        "https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam38.2/",
        "Pfam-A.clans.tsv.gz"
    ),
    destfile = clans_gz, mode = "wb"
)
clans <- read.delim(
    gzfile(clans_gz), header = FALSE, quote = "",
    col.names = c("Accession", "Clan", "Clan_name", "Pfam_name", "Description")
)
renamed <- c(
    DUF260 = "PF03195", DUF573 = "PF04504", E2F_TDP = "PF02319",
    EIN3 = "PF04873", FLO_LFY = "PF01698", Homeobox = "PF00046",
    `zf-Dof` = "PF02701"
)
acc <- clans$Accession[match(single$DBDs, clans$Pfam_name)]
acc[is.na(acc)] <- renamed[single$DBDs[is.na(acc)]]
stopifnot(!anyNA(acc))

cisbp_scheme <- data.frame(
    Family = single$Family_Name, CISBP_DBD = single$DBDs, Pfam = acc
)

# Pfam split FLO_LFY into SAM_LFY (PF01698) and the C-terminal DBD (PF17538)
cisbp_scheme <- rbind(cisbp_scheme, data.frame(
    Family = "Leafy", CISBP_DBD = "FLO_LFY", Pfam = "PF17538"
))
idx <- match(cisbp_scheme$Pfam, clans$Accession)
cisbp_scheme$Pfam_name <- clans$Pfam_name[idx]
cisbp_scheme$Pfam_description <- trimws(clans$Description[idx])
cisbp_scheme <- cisbp_scheme[order(
    tolower(cisbp_scheme$Family), cisbp_scheme$Pfam, method = "radix"
), ]
rownames(cisbp_scheme) <- NULL

usethis::use_data(cisbp_scheme, compress = "xz", overwrite = TRUE)
```

## tfclass_scheme.rda

Correspondence table between CIS-BP families and TFClass superclasses,
classes, and families. TFClass (downloaded on 2026-10-15) includes human
TFs and their mammalian orthologs, so the table was created by comparing
the CIS-BP families (build 3.10) of human TFs in CIS-BP with their
TFClass classifications. Human TFs in CIS-BP were matched to TFClass by
Ensembl gene ID or by gene symbol (1391 of 1639 TFs). For each CIS-BP
family with a single type of DBD, the TFClass class was the class of
most TFs (\>= 99% of the TFs for all families), and the TFClass family
was the family of \>= 90% of the TFs (otherwise, NA). The following
entries were manually curated:

1.  ‘AT hook’: TFClass family NA, as CIS-BP’s AT hook TFs include
    proteins other than HMGA (TFClass’ only family in class ‘A.T hook
    factors’).
2.  ‘CUT’ and ‘POU’: CUT and POU domains only occur in homeodomain TFs,
    so TFs in which only these domains were found are assigned to the
    TFClass families ‘HD-CUT’ and ‘POU’.
3.  CIS-BP families with multiple types of DBDs for which TFClass
    classifies human TFs based on one of the DBDs (e.g.,
    ‘Homeodomain,POU’).

The objects `tfinfo` and `cisbp_scheme` are created in the chunk
*cisbp_scheme* above.

``` r
# Download TFClass ontology, and parse TFClass nodes and human TFs (with
# genus and Ensembl gene ID)
ttl_path <- file.path(tempdir(), "tfclass.ttl.gz")
download.file(
    "http://tfclass.bioinf.med.uni-goettingen.de/suppl/tfclass.ttl.gz",
    destfile = ttl_path, mode = "wb"
)
ttl <- readLines(gzfile(ttl_path), encoding = "UTF-8")
starts <- grep("^<http://sybig\\.de/tfclass#", ttl)
ends <- c(starts[-1] - 1, length(ttl))
blocks <- vapply(seq_along(starts), function(i) {
    paste(ttl[starts[i]:ends[i]], collapse = "\n")
}, character(1))
ids <- sub("^<http://sybig\\.de/tfclass#([^>]+)>.*", "\\1", ttl[starts])
extract <- function(pattern, x) {
    m <- regmatches(x, regexec(pattern, x))
    vapply(m, function(y) if(length(y)) y[2] else NA_character_, character(1))
}
is_node <- grepl("^[0-9]+(\\.[0-9]+)*$", ids)
nodes <- data.frame(
    id = ids[is_node], label = extract('rdfs:label "([^"]*)"', blocks[is_node])
)
is_hs <- startsWith(ids, "Homo_sapiens_")
hs <- data.frame(
    factor = toupper(extract('rdfs:label "([^"]*)"', blocks[is_hs])),
    genus = extract(paste0(
        ":belongs_to ;\\s*owl:someValuesFrom ",
        "<http://sybig\\.de/tfclass#([0-9.]+)>"
    ), blocks[is_hs]),
    ensembl = extract("ENSEMBL_GeneID:(ENSG[0-9]+)", blocks[is_hs])
)
hs <- hs[!is.na(hs$genus), ]
glabel <- nodes$label[match(hs$genus, nodes$id)]
hs$symbol <- toupper(sub(".*\\(([^)]+)\\)$", "\\1", glabel))

# Match human TFs in CIS-BP to TFClass by Ensembl ID, then by gene symbol
hs_tfs <- tfinfo[tfinfo$TF_Species == "Homo_sapiens", ]
hs_tfs <- hs_tfs[!duplicated(hs_tfs$DBID), ]
idx <- match(hs_tfs$DBID, hs$ensembl)
miss <- is.na(idx)
idx[miss] <- match(toupper(hs_tfs$TF_Name[miss]), hs$symbol)
miss <- is.na(idx)
idx[miss] <- match(toupper(hs_tfs$TF_Name[miss]), hs$factor)
hs_tfs$genus <- hs$genus[idx]
hs_tfs <- hs_tfs[!is.na(hs_tfs$genus), ]
ancestor <- function(id, n) {
    vapply(strsplit(id, ".", fixed = TRUE), function(x) {
        paste(x[seq_len(min(n, length(x)))], collapse = ".")
    }, character(1))
}
hs_tfs$class <- ancestor(hs_tfs$genus, 2)
hs_tfs$family <- ancestor(hs_tfs$genus, 3)

# Single-DBD CIS-BP families (in cisbp_scheme): most frequent TFClass 
# class, and TFClass family if >= 90% of the TFs are in the same family
hs_single <- hs_tfs[hs_tfs$Family_Name %in% cisbp_scheme$Family, ]
map <- lapply(split(hs_single, hs_single$Family_Name), function(x) {
    cl <- sort(table(x$class), decreasing = TRUE)
    fam <- sort(table(x$family), decreasing = TRUE)
    data.frame(
        cisbp_family = x$Family_Name[1], class_id = names(cl)[1],
        family_id = ifelse(fam[[1]] / nrow(x) >= 0.9, names(fam)[1], NA)
    )
})
map <- do.call(rbind, map)

# Manual curation: (1) AT hook TFs in CIS-BP include proteins other than 
# HMGA (TFClass' only A.T hook family); (2) CUT and POU domains occur
# only in homeodomain TFs (TFClass families HD-CUT and POU), so TFs in
# which only CUT or POU domains were found are assigned to these families;
# (3) TFs with multiple types of DBDs, for which TFClass classifies TFs 
# based on one of the DBDs (based on human TFs in TFClass).
manual <- data.frame(
    cisbp_family = c(
        "AT hook", "CUT", "POU", "CUT,Homeodomain", "Homeodomain,POU",
        "C2H2 ZF,Homeodomain", "Homeodomain,Paired box", "bZIP,C2H2 ZF",
        "ARID/BRIGHT,RFX"
    ),
    class_id = c("8.2", "3.1", "3.1", "3.1", "3.1", "3.1", "3.2", "1.1", "3.7"),
    family_id = c(
        NA, "3.1.9", "3.1.10", "3.1.9", "3.1.10", "3.1.8", "3.2.1", 
        "1.1.1", "3.7.1"
    )
)
map <- rbind(map[!map$cisbp_family %in% manual$cisbp_family, ], manual)

# Add names of TFClass superclasses, classes, and families
label <- setNames(nodes$label, nodes$id)
map$superclass <- unname(label[ancestor(map$class_id, 1)])
map$class <- unname(label[map$class_id])
map$family <- unname(label[map$family_id])
map$tfclass_id <- ifelse(is.na(map$family_id), map$class_id, map$family_id)
map <- map[order(
    as.numeric_version(map$tfclass_id), map$cisbp_family, method = "radix"
), c("cisbp_family", "superclass", "class", "family", "tfclass_id")]
rownames(map) <- NULL
tfclass_map <- map

tfclass_scheme <- tfclass_map
names(tfclass_scheme) <- c(
    "CISBP_family", "TFClass_superclass", "TFClass_class", 
    "TFClass_family", "TFClass_ID"
)

usethis::use_data(tfclass_scheme, compress = "xz", overwrite = TRUE)
```

## human_tfs.rda

Curated human TFs from Lambert et al. (2018), downloaded from the Human
TFs database (version 1.01; <https://humantfs.ccbr.utoronto.ca/>). Only
proteins considered TFs by the curators were kept, and only categorical
annotations were included. CIS-BP families were obtained from CIS-BP
(build 3.10), whose human TFs are the TFs in this list. TFClass
classifications were obtained from TFClass for TFs in TFClass, and
inferred from CIS-BP families (as in `classify_tfs()`) for other TFs.

``` r
# Download the curated list of human TFs (Lambert et al., 2018)
lambert_path <- file.path(tempdir(), "DatabaseExtract_v_1.01.txt")
download.file(
    paste0(
        "https://humantfs.ccbr.utoronto.ca/download/v_1.01/",
        "DatabaseExtract_v_1.01.txt"
    ),
    destfile = lambert_path, mode = "wb"
)
lambert <- read.delim(lambert_path, check.names = FALSE, quote = "\"")
names(lambert) <- trimws(names(lambert))
lambert <- lambert[lambert$`Is TF?` == "Yes", ]

# CIS-BP families of human TFs (CIS-BP's human TFs are the curated TFs)
cisbp_hs <- tfinfo[tfinfo$TF_Species == "Homo_sapiens", ]
cisbp_hs <- cisbp_hs[!duplicated(cisbp_hs$DBID), ]
stopifnot(setequal(cisbp_hs$DBID, lambert$`Ensembl ID`))

human_tfs <- data.frame(
    Gene = lambert$`Ensembl ID`,
    Symbol = lambert$`HGNC symbol`,
    CISBP_family = cisbp_hs$Family_Name[
        match(lambert$`Ensembl ID`, cisbp_hs$DBID)
    ]
)

# TFClass classifications from TFClass 
idx <- match(human_tfs$Gene, hs_tfs$DBID)
label <- setNames(nodes$label, nodes$id)
from_tfclass <- data.frame(
    TFClass_superclass = unname(label[ancestor(hs_tfs$genus[idx], 1)]),
    TFClass_class = unname(label[hs_tfs$class[idx]]),
    TFClass_family = unname(label[hs_tfs$family[idx]])
)
inferred <- .classify_tfclass(human_tfs)[, names(from_tfclass)]
in_tfclass <- !is.na(idx)
human_tfs[names(from_tfclass)] <- inferred
human_tfs[in_tfclass, names(from_tfclass)] <- from_tfclass[in_tfclass, ]
human_tfs$TFClass_source <- ifelse(
    in_tfclass, "TFClass",
    ifelse(is.na(inferred$TFClass_class), NA, "Inferred from CIS-BP")
)

# Annotation from Lambert et al. (2018)
human_tfs$DBD <- lambert$DBD
human_tfs$TF_assessment <- lambert$`TF assessment`
human_tfs$Binding_mode <- lambert$`Binding mode`
human_tfs$Motif_status <- lambert$`Motif status`
human_tfs <- human_tfs[order(human_tfs$Symbol, method = "radix"), ]
rownames(human_tfs) <- NULL

usethis::use_data(human_tfs, compress = "xz", overwrite = TRUE)
```

## planttfclass_scheme.rda

Correspondence table between TAPscan and Plant-TFClass, as used by
`classify_tfs()` (see *Internal data* below for details on how it was
created). The code below must be run with `devtools::load_all()`.

``` r
planttfclass_scheme <- planttfclass_map[, c(
    "tapscan_family", "tapscan_subfamily", "superclass", "class", "family",
    "planttfdb_subfamily"
)]
names(planttfclass_scheme) <- c(
    "TAPscan_family", "TAPscan_subfamily", "PlantTFClass_superclass",
    "PlantTFClass_class", "PlantTFClass_family", "PlantTFDB_subfamily"
)

usethis::use_data(planttfclass_scheme, compress = "xz", overwrite = TRUE)
```

# Data in inst/extdata

## plantTFClass_to_tapscan.csv

Correspondence table between TAPscan families and Plant-TFClass
families, classes, and superclasses. Plant-TFClass was described in
Blanc-Mathieu et al. (2024), Trends in Plant Science, 29(1), 40-51.

## PlantTFDB.hmm.gz

The seed alignments for all TF families were downloaded from PFAM and
converted to profile HMMs using the command `hmmbuild` from HMMER (file
*PFAM.hmm*, code below). Domains that were self-built by PlantTFDB were
kindly provided by the database maintainers by email (file
*self_built.hmm*).

``` r
library(tidyverse)
domains <- c(
  AP2 = "PF00847", B3 = "PF02362", auxin_resp = "PF06507", 
  GAGA_bind = "PF06217", DUF822 = "PF05687", HLH = "PF00010",
  bZIP_1 = "PF00170", zf_B_box = "PF00643", CCT = "PF06203", 
  Zf_Dof = "PF02701", GATA_zf = "PF00320", Zf_LSD1 = "PF06943",
  Peptidase_C14 = "PF00656", YABBY = "PF04690", Zf_C2H2 =  "PF00096",
  RNase_T = "PF00929", Zf_CCCH =  "PF00642", RRM_1 = "PF00076",
  Helicase_C = "PF00271", CG1 = "PF03859", TCR = "PF03638",
  E2F_TDP = "PF02319", EIN3 = "PF04873", FAR1 = "PF03101",
  STAT = "PF02864", Response_reg = "PF00072", DUF573 = "PF04504",
  GRAS = "PF03514", WRC = "PF08879", QLQ = "PF08880",
  Homeobox = "PF00046", SMART = "PF01852", ELK = "PF03789",
  PHD = "PF00628", HSF = "PF00447", DUF260 = "PF03195",
  FLO_LFY = "PF01698", SRF_TF = "PF00319", K_box = "PF01486",
  MYB_DNAbind = "PF00249", SWIRM = "PF04433", NAM = "PF02365",
  Zf_NF_X1 = "PF01422", CBFB_NFYA = "PF02045", RWP_RK = "PF02042",
  NOZZLE = "PF08744", S1FA = "PF04689", SBP = "PF03110",
  DUF702 = "PF05142", TCP = "PF03634", WRKY = "PF03106",
  Whirly = "PF08536", ZF_HD_dimer = "PF04770"
)

# Define function to download seed alignments for PFAM domains
download_pfam_seeds <- function(domains = NULL) {
  
  d <- lapply(domains, function(x) {
    url <- paste0('http://pfam.xfam.org/family/', x,
                  '/alignment/seed/format?format=fasta')
    aln_path <- paste0(tempdir(), "/", x, ".fasta")
    if(!file.exists(aln_path)) {
      options(timeout = 5)
      result <- tryCatch({
        downloader::download(url = url, destfile = aln_path)
        res <- TRUE
      },
      error = function(e) return(FALSE),
      warning = function(w) return(FALSE)
      )
      
      if (!result) {
        file.remove(aln_path)
      }
    }
    return(NULL)
  })
  files <- list.files(tempdir(), pattern = ".fasta")
  ds <- sapply(domains, function(x) return(TRUE %in% 
                                             stringr::str_detect(files, x)))
  status <- data.frame(
    domain = domains,
    download_status = ds
  )
  return(status)
}

download_aln <- download_pfam_seeds(domains)


# Define function to build profile HMMs from alignments
build_profile_HMM <- function(domains = NULL) {
  
  f <- lapply(domains, function(x) {
    
    aln_path <- paste0(tempdir(), "/", x, ".fasta")
    hmm_path <- paste0(tempdir(), "/", x, ".hmm")
    build <- system2("hmmbuild", args = c(hmm_path, aln_path))
    sedargs <- paste0("s/aln$/", x, "/g")
    rename_aln <- system2("sed", args = c("-i", sedargs, hmm_path))
    return(NULL)
  })
  files <- list.files(tempdir(), pattern = ".hmm")
  build <- sapply(domains, function(x) return(TRUE %in% 
                                             stringr::str_detect(files, x)))
  status_df <- data.frame(
    domain = domains,
    build_status = build
  )
  return(status_df)
}

build_profile_HMM(domains)

# Combine .hmm files in a single one
combine_domains <- function() {
  files <- list.files(tempdir(), pattern = ".hmm", full.names = TRUE)
  files <- paste0(files, collapse = " ")
  out_hmm <- file.path(tempdir(), "PFAM.hmm")
  system2("cat", args = c(files, " > ", out_hmm))
  return(NULL)
}

combine_domains()
fs::file_move(path = file.path(tempdir(), "PFAM.hmm"),
              new_path = here::here("PFAM.hmm"))
```

Self-built HMMs (HMMER3/b format) were converted to the format of PFAM
HMMs (HMMER3/f), as all HMMs in a file must have the same format, and
both files were combined in a single file. Domain annotations obtained
with the combined file are identical to the ones obtained with each file
separately.

``` bash
hmmconvert self_built.hmm > self_built_3f.hmm
cat PFAM.hmm self_built_3f.hmm | gzip -9 -c > inst/extdata/PlantTFDB.hmm.gz
```

## TAPscan.hmm.gz

Profile HMMs used by TAPscan v4 (file *domains_v13.txt*) were obtained
from the TAPscan-classify repository
(<https://github.com/Rensing-Lab/TAPscan-classify>), and compressed with
`gzip -9`.

``` bash
gzip -9 -c TAPscan-classify/domains_v13.txt > inst/extdata/TAPscan.hmm.gz
```

## CISBP.hmm.gz

Profile HMMs of the Pfam families in `cisbp_scheme` (see above) were
extracted from the Pfam-A HMM library of Pfam release 38.2 with the
command `hmmfetch` from HMMER. HMM names were replaced with Pfam
accessions (without version), as in *PlantTFDB.hmm.gz*.

``` r
# Download Pfam-A HMMs (Pfam 38.2, ~420 MB)
pfam_gz <- file.path(tempdir(), "Pfam-A.hmm.gz")
download.file(
    paste0(
        "https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam38.2/",
        "Pfam-A.hmm.gz"
    ),
    destfile = pfam_gz, mode = "wb"
)

# Extract HMMs by name, in order of accession
pfam_hmm <- file.path(tempdir(), "Pfam-A.hmm")
system2("gzip", c("-dc", pfam_gz), stdout = pfam_hmm)
system2("hmmfetch", c("--index", pfam_hmm), stdout = FALSE)
scheme <- cisbp_scheme[order(cisbp_scheme$Pfam), ]
keys <- file.path(tempdir(), "cisbp_keys.txt")
writeLines(scheme$Pfam_name, keys)
hmm <- system2("hmmfetch", c("-f", pfam_hmm, keys), stdout = TRUE)

# Replace HMM names with Pfam accessions (without version)
name_idx <- which(startsWith(hmm, "NAME "))
acc <- sub("^ACC +(PF[0-9]+)\\..*", "\\1", hmm[name_idx + 1])
stopifnot(identical(acc, scheme$Pfam))
hmm[name_idx] <- paste0("NAME  ", acc)

con <- gzfile("inst/extdata/CISBP.hmm.gz", open = "w", compression = 9)
writeLines(hmm, con)
close(con)
```

# Internal data (sysdata.rda)

This data frame contains the sequence and domain cutoff for each family.
The values were retrieved from the Supplementary Materials of the
PlantTFDB v2 paper.

``` r
# Domain cutoffs for each domain (NA: no cutoff, use E-value instead)
cutoffs <- c(
    `PF00010` = 11, `PF00046` = 15, `PF00072` = 20, `PF00076` = NA,
    `PF00096` = 12.2, `PF00170` = 19.7, `PF00249` = 21.5, `PF00271` = NA,
    `PF00319` = 21, `PF00320` = 22, `PF00447` = 22.5, `PF00628` = NA,
    `PF00642` = 17, `PF00643` = 16, `PF00656` = NA, `PF00847` = NA,
    `PF00929` = NA, `PF01422` = 16, `PF01486` = NA, `PF01698` = 20.7,
    `PF01852` = NA, `PF02042` = 20.8, `PF02045` = 22, `PF02319` = 23,
    `PF02362` = 30.2, `PF02365` = 21.5, `PF02701` = 26, `PF02864` = 150,
    `PF03101` = 21.3, `PF03106` = 22, `PF03110` = 23, `PF03195` = 22,
    `PF03514` = 20.5, `PF03634` = 21, `PF03638` = 20.9, `PF03789` = NA,
    `PF03859` = 25, `PF04433` = NA, `PF04504` = 29.5, `PF04689` = 22,
    `PF04690` = 21, `PF04770` = 20.5, `PF04873` = 22, `PF05142` = 21.3,
    `PF05687` = 25, `PF06203` = NA, `PF06217` = 23, `PF06507` = NA,
    `PF06943` = 20.2, `PF08536` = 20.2, `PF08744` = 20.2, `PF08879` = 21.3,
    `PF08880` = NA, `G2-like` = 20, `HRT-like` = 24, `NF-YB` = 27,
    `NF-YC` = 24, `SAP` = 50, `STAT` = 150, `VOZ` = 100,
    `HD-ZIP_I/II` = NA, `BELL` = NA, `Wus_type_Homeobox` = NA, `trihelix` = 16
)

# Check that there is a cutoff for each profile HMM
lines <- readLines(gzfile(file.path("inst", "extdata", "PlantTFDB.hmm.gz")))
hmm_names <- gsub("^NAME +", "", lines[startsWith(lines, "NAME")])
stopifnot(setequal(names(cutoffs), hmm_names))

score_cutoff <- data.frame(
    domain = names(cutoffs), domaincutoff = unname(cutoffs)
)
```

Classification rules (v82) and domain-specific coverage cutoffs (v11)
used by TAPscan were obtained from the TAPscan-classify repository
(<https://github.com/Rensing-Lab/TAPscan-classify>). Rules are sorted as
in the original Perl script (i.e., in C locale), as the order of
families affects how genes with multiple candidate families are
classified.

``` r
rules <- readLines("TAPscan-classify/rules_v82.txt")
rules <- strsplit(sort(rules, method = "radix"), ";")
tapscan_rules <- data.frame(
    family = vapply(rules, `[`, character(1), 1),
    domain = vapply(rules, `[`, character(1), 2),
    type = vapply(rules, `[`, character(1), 3)
)

tapscan_coverage <- read.delim(
    "TAPscan-classify/coverage_values_v11.txt", header = FALSE,
    col.names = c("domain", "cutoff")
)
```

The correspondence table between TAPscan families and Plant-TFClass
families, classes, and superclasses is available in
*inst/extdata/plantTFClass_to_tapscan.csv*. This table was converted to
a tidy data frame with one row per TAPscan family/subfamily and
Plant-TFClass classification. The following rules were used:

1.  Entries in the format ‘Family:Subfamily’ (e.g., ‘MYB:MYB-2R’) match
    TAPscan families and subfamilies, and they take precedence over
    other entries.
2.  Other entries match TAPscan families (or subfamilies, when
    available) with the same name. If there is no family or subfamily
    with the same name, entries match all subfamilies of a TAPscan
    family (e.g., ‘HDZ’ matches C1HDZ, C2HDZ, C3HDZ, and C4HDZ).
3.  For each TAPscan family/subfamily, only the most specific matches
    are kept.
4.  As stated in the table, RAV TFs are classified as AP2 or ABI3/VP1 by
    TAPscan, so TAPscan families AP2 and ABI3/VP1 were also matched to
    Plant-TFClass family RAV.
5.  Plant-TFClass classes without families have the class name as
    family.
6.  When a TAPscan family/subfamily matches multiple Plant-TFClass
    families, PlantTFDB subfamilies (variable *planttfdb_subfamily*) are
    used to choose among them.

The code below must be run with `devtools::load_all()`, as it uses
internal functions and data from this package.

``` r
raw <- read.csv(
    "inst/extdata/plantTFClass_to_tapscan.csv",
    check.names = FALSE, na.strings = c("", "NA")
)

# RAV TFs are classified by TAPscan as AP2 or ABI3/VP1 (see note in the table)
rav <- raw$Family %in% "RAV"
raw[rav, 6:9] <- list("AP2", "ABI3/VP1", NA, NA)

# All possible TAPscan family/subfamily pairs
fams <- unique(tapscan_rules$family)
pairs <- unique(.tapscan_family_map(fams))
pairs$name <- ifelse(is.na(pairs$Subfamily), pairs$Family, pairs$Subfamily)

# Expand each TAPscan token in the table to TAPscan family/subfamily pairs.
# Specificity: 2 = 'Family:Subfamily'; 1 = name of family (or subfamily, if
# any) as in TAPscan's family statistics; 0 = parent family of subfamilies.
expand_token <- function(token) {
    if(grepl(":", token)) {
        x <- strsplit(token, ":")[[1]]
        m <- pairs[pairs$Family == x[1] & pairs$Subfamily %in% x[2], 1:2]
        spec <- rep(2, nrow(m))
    } else {
        m1 <- pairs$name == token
        m0 <- pairs$Family == token & !m1
        m <- pairs[m1 | m0, 1:2]
        spec <- ifelse(m1, 1, 0)[m1 | m0]
    }
    data.frame(m, specificity = spec)
}

map <- do.call(rbind, lapply(seq_len(nrow(raw)), function(i) {
    tokens <- unlist(raw[i, 6:9])
    tokens <- trimws(tokens[!is.na(tokens)])
    m <- do.call(rbind, lapply(tokens, expand_token))
    if(is.null(m) || nrow(m) == 0) { return(NULL) }
    data.frame(
        tapscan_family = m$Family, tapscan_subfamily = m$Subfamily,
        specificity = m$specificity,
        family = raw$Family[i], class = raw$Class[i],
        superclass = raw$SuperClass[i]
    )
}))

# Report tokens that do not match any TAPscan family
all_tokens <- trimws(unlist(raw[, 6:9]))
all_tokens <- unique(all_tokens[!is.na(all_tokens)])
unmatched <- all_tokens[vapply(all_tokens, function(t) nrow(expand_token(t)) == 0, logical(1))]
message("Unmatched tokens: ", paste(unmatched, collapse = " | "))

# Keep only the most specific match for each TAPscan family/subfamily pair
key <- paste(map$tapscan_family, map$tapscan_subfamily)
max_spec <- ave(map$specificity, key, FUN = max)
map <- map[map$specificity == max_spec, ]
map <- unique(map[, c("tapscan_family", "tapscan_subfamily", "family", "class", "superclass")])

# Classes without families have the class name as family
map$family <- ifelse(is.na(map$family), map$class, map$family)

# PlantTFDB subfamilies used to resolve TAPscan families that match multiple
# Plant-TFClass families
resolve <- c(
    AP2 = "AP2", `ERF/DREB` = "ERF", RAV = "RAV", LAV = "B3",
    `Type I` = "M-type", `Type II` = "MIKC"
)
map$planttfdb_subfamily <- unname(resolve[map$family])
map <- map[order(map$superclass, map$class, map$family, map$tapscan_family, method = "radix"), ]
rownames(map) <- NULL
planttfclass_map <- map

# Lookup table of TAP classes for each TAPscan family/subfamily
# (from tap_class.rda, see above)
tap_class_map <- tap_class[, c("TAP_family", "TAP_class")]

# Lookup table of PlantTFDB families for each subfamily
# (from planttfdb_scheme.rda, see above)
planttfdb_families <- unique(planttfdb_scheme[, c("Family", "Subfamily")])
rownames(planttfdb_families) <- NULL

# Lookup table of CIS-BP families for each DBD (from cisbp_scheme.rda, 
# see above). The correspondence table between CIS-BP and TFClass
# (tfclass_map) was created in the chunk tfclass_scheme above.
cisbp_map <- data.frame(pfam = cisbp_scheme$Pfam, family = cisbp_scheme$Family)

usethis::use_data(
    score_cutoff, tapscan_rules, tapscan_coverage, planttfclass_map,
    tap_class_map, planttfdb_families, cisbp_map, tfclass_map,
    compress = "xz", internal = TRUE, overwrite = TRUE
)
```
