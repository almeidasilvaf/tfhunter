#' Protein sequences of the algae species Galdieria sulphuraria
#'
#' Data obtained from PLAZA Diatoms. Only genes containing domains used for
#' TF family classification (with any of the schemes in
#' \code{list_schemes()}) were kept for package size issues.
#'
#' @name gsu
#' @format An AAStringSet object as returned 
#' by \code{Biostrings::readAAStringSet()}.
#' @references
#' Osuna-Cruz, C. M., Bilcke, G., Vancaester, E., De Decker, S., Bones, A. M., 
#' Winge, P., ... & Vandepoele, K. (2020). The Seminavis robusta genome 
#' provides insights into the evolutionary adaptations of benthic diatoms. 
#' Nature communications, 11(1), 1-13.
#' @examples
#' data(gsu)
#' @usage data(gsu)
"gsu"


#' Data frame of PlantTFDB's TF family classification scheme
#'
#' The classification scheme is the same as the one used by PlantTFDB,
#' except for three family/subfamily names that were changed to match the
#' output of \code{classify_tfs()} ('M_type' to 'M-type', 'MYB_related' to
#' 'MYB-related', and 'LBD (AS2/LOB)' to 'LBD'). As in PlantTFDB, families
#' without subfamilies have the family name as subfamily.
#' 
#' @name planttfdb_scheme
#' @format A data frame with the following variables:
#' \describe{
#'  \item{Family}{TF family name.}
#'  \item{Subfamily}{TF subfamily name.}
#'  \item{DBD}{DNA-binding domain}
#'  \item{Auxiliary}{Auxiliary domain}
#'  \item{Forbidden}{Forbidden domain}
#' }
#' @references 
#' Jin, J., Tian, F., Yang, D. C., Meng, Y. Q., Kong, L., Luo, J., & 
#' Gao, G. (2016). PlantTFDB 4.0: toward a central hub for transcription 
#' factors and regulatory interactions in plants. 
#' Nucleic acids research, gkw982.
#' @examples 
#' data(planttfdb_scheme)
#' @usage data(planttfdb_scheme)
"planttfdb_scheme"


#' Domain annotation for the algae species Galdieria sulphuraria
#'  
#' The data set was created using the function \code{annotate_domains()}.
#'
#' @name gsu_annotation
#' @format A list of 3 data frames (\strong{PlantTFDB}, \strong{TAPscan},
#' and \strong{CISBP}) with domain annotation, as returned
#' by \code{annotate_domains()}.
#' @examples
#' data(gsu_annotation)
#' @usage data(gsu_annotation)
"gsu_annotation"


#' Classes of transcription-associated proteins (TAPs) in TAPscan
#'
#' Families of transcription-associated proteins (TAPs) in TAPscan v4, their
#' classes, descriptions, and references. Data were obtained from TAPscan's
#' website (file \emph{_data/import-tapinfo/tapinfo_v4.csv} in the
#' GitHub repository Rensing-Lab/TAPscan-v4-website). HTML tags in
#' descriptions and references were removed. Family MADS_MIKC, which is not
#' included in the original file, was manually added as a TF family (as
#' family MADS), without description and references.
#'
#' @name tap_class
#' @format A data frame with the following variables:
#' \describe{
#'  \item{TAP_family}{TAP family or subfamily, as in the variable
#'  \strong{TAPscan_subfamily} of the output of \code{classify_tfs()}.}
#'  \item{TAP_class}{TAP class. One of 'TF' (transcription factor),
#'  'TR' (transcriptional regulator), or 'PT' (putative TAP).}
#'  \item{Description}{Description of the TAP family.}
#'  \item{References}{References for the TAP family, separated by ' | '.}
#'}
#' @references
#' Petroll, R., Varshney, D., Hiltemann, S., Finke, H., Schreiber, M.,
#' de Vries, J., & Rensing, S. A. (2025). Enhanced sensitivity of TAPscan v4
#' enables comprehensive analysis of streptophyte transcription factor
#' evolution. The Plant Journal, 121(1), e17184.
#' @examples
#' data(tap_class)
#' @usage data(tap_class)
"tap_class"


#' Data frame of TAPscan's classification scheme
#'
#' TAPscan v4 classification rules, with one row per rule. Genes are
#' assigned to a family (and subfamily, if any) if they have all domains
#' in variable \strong{Should} and none of the domains in variable
#' \strong{Should_not}. For families with more than one rule (e.g., bZIP),
#' genes are assigned to the family if they satisfy any of the rules.
#'
#' @name tapscan_scheme
#' @format A data frame with the following variables:
#' \describe{
#'  \item{Family}{TAPscan family, as in variable \strong{TAPscan_family} of
#'  the output of \code{classify_tfs()}.}
#'  \item{Subfamily}{TAPscan subfamily, as in variable
#'  \strong{TAPscan_subfamily} of the output of \code{classify_tfs()}.
#'  Families without subfamilies have the family name as subfamily.}
#'  \item{Rule}{Name of the rule in TAPscan's rules file.}
#'  \item{Should}{Domains that genes should have, separated by ' and '.
#'  Numbers in parentheses (e.g., 'Myb_DNA-binding (2)') indicate the
#'  required number of copies of a domain.}
#'  \item{Should_not}{Domains that genes should not have, separated
#'  by ' or '.}
#'  \item{TAP_class}{TAP class. One of 'TF' (transcription factor),
#'  'TR' (transcriptional regulator), or 'PT' (putative TAP). NA for
#'  families without a TAP class in TAPscan (see \code{tap_class}).}
#'}
#' @references
#' Petroll, R., Varshney, D., Hiltemann, S., Finke, H., Schreiber, M.,
#' de Vries, J., & Rensing, S. A. (2025). Enhanced sensitivity of TAPscan v4
#' enables comprehensive analysis of streptophyte transcription factor
#' evolution. The Plant Journal, 121(1), e17184.
#' @examples
#' data(tapscan_scheme)
#' @usage data(tapscan_scheme)
"tapscan_scheme"


#' Correspondence table between TAPscan and Plant-TFClass
#'
#' Table used by \code{classify_tfs()} to infer Plant-TFClass superclasses,
#' classes, and families from TAPscan families and subfamilies.
#' TAPscan families/subfamilies with multiple rows correspond to more than
#' one Plant-TFClass family, and the PlantTFDB subfamily of each gene
#' (variable \strong{PlantTFDB_subfamily}) is used to choose among them.
#'
#' @name planttfclass_scheme
#' @format A data frame with the following variables:
#' \describe{
#'  \item{TAPscan_family}{TAPscan family.}
#'  \item{TAPscan_subfamily}{TAPscan subfamily (or family, for families
#'  without subfamilies).}
#'  \item{PlantTFClass_superclass}{Plant-TFClass superclass.}
#'  \item{PlantTFClass_class}{Plant-TFClass class.}
#'  \item{PlantTFClass_family}{Plant-TFClass family (or class, for classes
#'  without families).}
#'  \item{PlantTFDB_subfamily}{PlantTFDB subfamily used to choose among
#'  multiple Plant-TFClass families for the same TAPscan family/subfamily.
#'  NA for TAPscan families/subfamilies with a single Plant-TFClass family.}
#'}
#' @references
#' Blanc-Mathieu, R., Dumas, R., Turchi, L., Lucas, J., & Parcy, F. (2024).
#' Plant-TFClass: a structural classification for plant transcription
#' factors. Trends in Plant Science, 29(1), 40-51.
#' @examples
#' data(planttfclass_scheme)
#' @usage data(planttfclass_scheme)
"planttfclass_scheme"


#' Data frame of CIS-BP's TF family classification scheme
#'
#' DNA-binding domains (DBDs) used by CIS-BP to identify TFs, and the
#' TF family corresponding to each DBD. Genes with more than one type of DBD
#' are assigned to a family whose name is the family names of all DBDs
#' sorted alphabetically (case-insensitive) and separated by ',' (e.g.,
#' 'Homeodomain,POU'), as in CIS-BP. DBDs and families were obtained from
#' CIS-BP (build 3.10), and DBDs were matched to Pfam (release 38.2)
#' accessions. Some Pfam families were renamed since they were used by
#' CIS-BP (e.g., 'DUF260' is now 'LOB'), and the Pfam family
#' 'FLO_LFY' was split into 'SAM_LFY' (PF01698) and 'C_LFY_FLO' (PF17538),
#' both of which are used for family 'Leafy'. The RNA-binding domains
#' 'RRM_1' and 'zf-CCCH' are not included, as CIS-BP classifies TFs in
#' families 'RRM' and 'CCCH ZF' only if there is experimental evidence
#' of DNA binding (6 TFs in CIS-BP 3.10), and never uses these domains for
#' TFs with other DBDs.
#'
#' @name cisbp_scheme
#' @format A data frame with the following variables:
#' \describe{
#'  \item{Family}{CIS-BP family of TFs with the DBD.}
#'  \item{CISBP_DBD}{DBD name in CIS-BP.}
#'  \item{Pfam}{Pfam accession of the DBD.}
#'  \item{Pfam_name}{Pfam name of the DBD (Pfam 38.2).}
#'  \item{Pfam_description}{Pfam description of the DBD.}
#' }
#' @references
#' Weirauch, M. T., Yang, A., Albu, M., Cote, A. G., Montenegro-Montero, A.,
#' Drewe, P., ... & Hughes, T. R. (2014). Determination and inference of
#' eukaryotic transcription factor sequence specificity. Cell, 158(6),
#' 1431-1443.
#' @examples
#' data(cisbp_scheme)
#' @usage data(cisbp_scheme)
"cisbp_scheme"


#' Correspondence table between CIS-BP and TFClass
#'
#' Table used by \code{classify_tfs()} to infer TFClass superclasses,
#' classes, and families from CIS-BP families. The table was created by
#' comparing the CIS-BP families of human TFs in CIS-BP (build 3.10) with
#' their classifications in TFClass. CIS-BP families were assigned to the
#' TFClass class of most of their TFs, and to a TFClass family only 
#' if >= 90\% of their TFs are in the same TFClass family (otherwise, the
#' TFClass family is NA). Some entries were manually curated, including
#' CIS-BP families with multiple types of DNA-binding domains (e.g.,
#' 'Homeodomain,POU'). CIS-BP families that are not in this table (e.g.,
#' families that do not occur in mammals, or with DNA-binding domains that
#' are not classified by TFClass) are not classified by TFClass, except
#' for families with multiple types of DNA-binding domains (see Details in
#' \code{?classify_tfs}).
#'
#' @name tfclass_scheme
#' @format A data frame with the following variables:
#' \describe{
#'  \item{CISBP_family}{CIS-BP family.}
#'  \item{TFClass_superclass}{TFClass superclass.}
#'  \item{TFClass_class}{TFClass class.}
#'  \item{TFClass_family}{TFClass family (NA if the CIS-BP family
#'  corresponds to more than one TFClass family).}
#'  \item{TFClass_ID}{TFClass identifier of the most specific
#'  classification level (class or family).}
#'}
#' @references
#' Wingender, E., Schoeps, T., Haubrock, M., Krull, M., & Dönitz, J. (2018).
#' TFClass: expanding the classification of human transcription factors to
#' their mammalian orthologs. Nucleic Acids Research, 46(D1), D343-D347.
#' @examples
#' data(tfclass_scheme)
#' @usage data(tfclass_scheme)
"tfclass_scheme"


#' Curated human transcription factors and their classifications
#'
#' Human TFs curated by Lambert et al. (2018), classified with the CIS-BP and
#' TFClass schemes. As this list was manually curated by experts, it is
#' more accurate than TFs identified from the human proteome with
#' \code{annotate_domains()} and \code{classify_tfs()}, and it can be used
#' as the set of human TFs in comparative analyses (e.g., with
#' \code{get_tf_counts()}). CIS-BP families were obtained from CIS-BP
#' (build 3.10), whose human TFs are the TFs in this list. TFClass
#' classifications were obtained from TFClass for TFs in TFClass, and
#' inferred from CIS-BP families (as in \code{classify_tfs()}) for
#' other TFs.
#'
#' @name human_tfs
#' @format A data frame with one row per TF and the following variables:
#' \describe{
#'  \item{Gene}{Ensembl gene ID (or UniProt entry name, for DUX1 and DUX3,
#'  which are not in Ensembl).}
#'  \item{Symbol}{HGNC gene symbol.}
#'  \item{CISBP_family}{TF family according to CIS-BP. TFs without a
#'  known DNA-binding domain are in family 'Unknown'.}
#'  \item{TFClass_superclass}{TF superclass according to TFClass.}
#'  \item{TFClass_class}{TF class according to TFClass.}
#'  \item{TFClass_family}{TF family according to TFClass.}
#'  \item{TFClass_source}{Source of TFClass classifications. One of
#'  'TFClass' (TF classified by TFClass) or 'Inferred from CIS-BP'
#'  (TFClass classifications inferred from the CIS-BP family). NA for TFs
#'  without TFClass classifications.}
#'  \item{DBD}{DNA-binding domain(s), as in Lambert et al. (2018).}
#'  \item{TF_assessment}{Assessment of DNA-binding activity. One of
#'  'Known motif', 'Inferred motif' (motif inferred from a closely
#'  related TF), or 'Likely to be sequence specific TF'.}
#'  \item{Binding_mode}{Mode of interacting with DNA. One of 'Monomer or
#'  homomultimer', 'Obligate heteromer', or 'Low specificity DNA-binding
#'  protein'.}
#'  \item{Motif_status}{Availability of DNA-binding motifs. One of
#'  'High-throughput in vitro', '100 perc ID - in vitro' (motif available
#'  for a TF with an identical DNA-binding domain), 'In vivo/Misc source',
#'  or 'No motif'.}
#'}
#' @references
#' Lambert, S. A., Jolma, A., Campitelli, L. F., Das, P. K., Yin, Y.,
#' Albu, M., ... & Weirauch, M. T. (2018). The human transcription factors.
#' Cell, 172(4), 650-665.
#'
#' Weirauch, M. T., Yang, A., Albu, M., Cote, A. G., Montenegro-Montero, A.,
#' Drewe, P., ... & Hughes, T. R. (2014). Determination and inference of
#' eukaryotic transcription factor sequence specificity. Cell, 158(6),
#' 1431-1443.
#'
#' Wingender, E., Schoeps, T., Haubrock, M., Krull, M., & Dönitz, J. (2018).
#' TFClass: expanding the classification of human transcription factors to
#' their mammalian orthologs. Nucleic Acids Research, 46(D1), D343-D347.
#' @examples
#' data(human_tfs)
#' head(human_tfs)
#'
#' # Count human TFs per TFClass class
#' se <- get_tf_counts(
#'     list(Human = human_tfs), scheme = "TFClass", level = "class"
#' )
#' @usage data(human_tfs)
"human_tfs"
