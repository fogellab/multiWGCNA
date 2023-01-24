#' Autism dataset
#'
#' The autism dataset from Voineagu et al. 2011, extracted from the Gene Expression Omnibus (accession number GSE28521). 
#' See ?sampleTable for metadata. 
#'
#' @docType data
#' @keywords datasets
#' @name autism_data
#' @usage data(autism_data)
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
#' @name autism_metadata
#' @usage data(autism_metadata)
#' @format A data frame with 58 rows (samples) and 3 columns (metadata)
NULL
#' Astrocyte dataset
#'
#' The astrocyte dataset from Itoh et al. PNAS. 2018, extracted from the supplementary materials.  
#' See astrocyte_metadata for metadata. 
#'
#' @docType data
#' @keywords datasets
#' @name astrocyte_data
#' @usage data(astrocyte_data)
#' @format An expression data frame with 20352 rows (genes) and 37 columns (gene names + 36 samples)
NULL

#' Astrocyte dataset metadata
#'
#' The metadata for the autism dataset (Voineagu et al. 2011). See ?datExpr for expression data. 
#'
#'#' \itemize{
#'   \item Sample the sample name corresponding to columns of datExpr
#'   \item Status experimental autoimmune encephalitis (EAE) or wildtype (WT)
#'   \item Tissue cerebellum (Cbl), cortex (Ctx), hippocampus (Hippo), or spinal cord (Sc)
#' }
#' 
#' @docType data
#' @keywords datasets
#' @name astrocyte_data
#' @usage data(astrocyte_data)
#' @format A data frame with 36 rows (samples) and 3 columns (metadata)
NULL