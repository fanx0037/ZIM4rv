#' An example dataset of genedata
#'
#' Small, artificially generated toy data set that
#' provides artificial information of genotypes for 200 individuals
#' on 3 rs locations to illustrate the analysis with the use of the package.
#'
#' @docType data
#'
#' @usage data(Ex1_genedata)
#'
#' @format An object of class \code{"data.frame"}
#' \describe{
#'  \item{FID}{Family IDs}
#'  \item{IID}{Individual IDs}
#'  \item{rs1}{Genotype code for rs1}
#'  \item{rs2}{Genotype code for rs2}
#'  \item{rs3}{Genotype code for rs3}
#' }
#' @references This data set was artificially created and modified for the ZIM4rv package.
#' @keywords datasets
#' @examples
#'
#' data(Ex1_genedata)
#' head(Ex1_genedata)
#'
"Ex1_genedata"
