#' Autism dataset
#'
#' The autism dataset from Voineagu et al. 2011, extracted from the Gene Expression Omnibus (accession number GSE28521). 
#' See ?sampleTable for metadata. 
#'
#' @docType data
#' @keywords datasets
#' @name datExpr
#' @usage data(datExpr)
#' @format An expression data frame with 9934 rows (genes) and 59 columns (probe names + 58 samples)
NULL

#' Autism dataset metadata
#'
#' The metadata for the autism dataset (Voineagu et al. 2011). See ?datExpr for expression data. 
#'
#'#' \itemize{
#'   \item Sample the sample name corresponding to columns of datExpr
#'   \item Status autism or control
#'   \item Tissue frontal cortex (FC) or temporal cortex (TC)
#' }
#' 
#' @docType data
#' @keywords datasets
#' @name sampleTable
#' @usage data(sampleTable)
#' @format A data frame with 58 rows (samples) and 3 columns (metadata)
NULL