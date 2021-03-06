
#' Protein sequences of the algae species Galdieria sulphuraria
#'
#' Data obtained from PLAZA Diatoms. Only genes containing domains used for
#' TF family classification were kept for package size issues.
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


#' Data frame of TF family classification scheme
#'
#' The classification scheme is the same as the one used by PlantTFDB.
#' 
#' @name classification_scheme
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
#' data(classification_scheme)
#' @usage data(classification_scheme)
"classification_scheme"


#' Data frame of domain cutoff for each domain
#'
#' The cutoffs is the same as the one used by PlantTFDB. When domain cutoffs
#' are not specified (NA), the evalue cutoff in the \strong{evalue} parameter
#' of \code{annotate_pfam} is used.
#' 
#' @name score_cutoff
#' @format A data frame with the following variables:
#' \describe{
#'  \item{domain}{Character with the domain name.}
#'  \item{domaincutoff}{Numeric with the domain score cutoff.}
#' }
#' @references 
#' Jin, J., Tian, F., Yang, D. C., Meng, Y. Q., Kong, L., Luo, J., & 
#' Gao, G. (2016). PlantTFDB 4.0: toward a central hub for transcription 
#' factors and regulatory interactions in plants. 
#' Nucleic acids research, gkw982.
#' @examples 
#' data(score_cutoff)
#' @usage data(score_cutoff)
"score_cutoff"


#' Domain annotation for the algae species Galdieria sulphuraria
#'  
#' The data set was created using the funcion \code{annotate_pfam()} in local
#' mode.
#'
#' @name gsu_annotation
#' @format A 2-column data frame with the following variables:
#' \describe{
#'  \item{Gene}{Gene ID}
#'  \item{Annotation}{Domain ID or domain name when ID is not available in PFAM}
#'}
#' @examples
#' data(gsu_annotation)
#' @usage data(gsu_annotation)
"gsu_annotation"


#' TFs families of the algae species Galdieria sulphuraria
#'  
#' The data set was created using the funcion \code{classify_tfs()}.
#'
#' @name gsu_families
#' @format A 2-column data frame with the following variables:
#' \describe{
#'  \item{Gene}{Gene ID}
#'  \item{Family}{TF family}
#'}
#' @examples
#' data(gsu_families)
#' @usage data(gsu_families)
"gsu_families"

