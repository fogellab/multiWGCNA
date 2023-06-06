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

#' Autism Summarized Experiment
#'
#' The autism dataset from Voineagu et al. 2011, extracted from the Gene Expression Omnibus (accession number GSE28521). 
#' See ?sampleTable for metadata. 
#'
#' @docType data
#' @keywords datasets
#' @name autism_se
#' @usage data(autism_se)
#' @format An expression data frame with 9934 rows (genes) and 59 columns (probe names + 58 samples)
NULL

#' Autism dataset metadata
#'
#' The metadata for the autism dataset (Voineagu et al. 2011). See ?autism_data for expression data. 
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

#' Astrocyte Ribotag Summarized Experiment
#'
#' The astrocyte ribotag dataset from Itoh et al. PNAS. 2018, extracted from the supplementary materials, 
#' provided as a SummarizedExperiment object. See ?astrocyte_networks for pre-computed WGCNA networks.  
#'
#' @docType data
#' @keywords datasets
#' @name astrocyte_se
#' @usage data(astrocyte_se)
#' @format An SummarizedExperiment object with 20352 genes and 36 samples
NULL

#' Astrocyte dataset metadata
#'
#' The metadata for the astrocyte Ribotag dataset (Itoh et al. PNAS. 2018). See ?astrocyte_data for expression data. 
#'
#'#' \itemize{
#'   \item Sample the sample name corresponding to columns of datExpr
#'   \item Status experimental autoimmune encephalitis (EAE) or wildtype (WT)
#'   \item Tissue cerebellum (Cbl), cortex (Ctx), hippocampus (Hippo), or spinal cord (Sc)
#' }
#' 
#' @docType data
#' @keywords datasets
#' @name astrocyte_metadata
#' @usage data(astrocyte_metadata)
#' @format A data frame with 36 rows (samples) and 3 columns (metadata)
NULL

#' Astrocyte networks
#'
#' The astrocyte networks from Tommasini and Fogel, BMC Bioinformatics, 2023. 
#' derived from the astrocyte Ribotag data from Itoh et al. PNAS. 2018. 
#'
#' @docType data
#' @keywords datasets
#' @name astrocyte_networks
#' @usage data(astrocyte_networks)
#' @format A list of WGCNA objects (level 1: combined; level 2: EAE, WT; and level 3: Cbl, Ctx, Hippo, Sc)
NULL

#' Tau pathology (rTg4510) Mouse Summarized Experiment
#'
#' The tau pathology dataset from Castanho et al. Cell Rep. 2020, extracted from GSE125957, 
#' processed as described in methods (https://git.exeter.ac.uk/ic322/ad-mice-rna-seq-cell-reports), 
#' and provided as a SummarizedExperiment object. Values correspond to DESeq2 rlog-transformed 
#' counts. See ?tau_networks for pre-computed WGCNA networks.  
#'
#' @docType data
#' @keywords datasets
#' @name tau_se
#' @usage data(tau_se)
#' @format A SummarizedExperiment object with 18822 genes and 58 samples
NULL

# Generation of summarized experiment objects
# tau_data = read.csv("../combined_summary.csv")
# tau_trait = read.csv("../traitData.csv")
# tau_mat = tau_data[,-c((ncol(tau_data)-4):ncol(tau_data))]
# rownames(tau_mat) = tau_data$X
# tau_colData = tau_trait
# rownames(tau_colData) = tau_trait$Sample
# tau_rowRanges = tau_data["X"]
# rownames(tau_rowRanges) = tau_data$X
# SummarizedExperiment(assays=list(counts=as.matrix(tau_mat)),
#                      rowData=tau_rowRanges, colData=tau_colData)
