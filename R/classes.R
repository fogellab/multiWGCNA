#' The WGCNA Class
#'
#' The WGCNA class is the main class used in multiWGCNA to store results from 
#' a weighted gene co-expression nework analysis. These include the original 
#' unaltered expression data used as input, connectivity metrics, module 
#' assignment, input sample conditions, trait 
#'
#' @slot datExpr The expression data, connectivity data, and module assignment
#' @slot conditions A data.frame with integer conditions for WGCNA
#' @slot trait A data.frame showing pearson correlation values to traits
#' @slot moduleEigengenes A data.frame of module eigengenes for each module across samples
#' @slot outlierModules A vector of modules classified by our algorithm as being driven by sample outliers
#'
#' @return NA
#'
#' @import methods
setClass("WGCNA", slots=list(datExpr="data.frame", 
                             conditions="data.frame", 
                             trait="data.frame", 
                             moduleEigengenes="data.frame", 
                             outlierModules="vector"))

setMethod("show", "WGCNA", function(object) {
		cat("##### datExpr #####\n")
		print(head(object@datExpr[,1:5]))
		if(length(object@conditions)>0){
			cat("\n##### conditions #####\n")
			print(head(object@conditions))
		}
		if(length(object@conditions)>0){
			cat("\n##### module-trait correlation #####\n")
			print(head(object@trait))
		}
		if(length(object@moduleEigengenes)>0){
			cat("\n##### module eigengenes #####\n")
			print(head(object@moduleEigengenes[,1:5]))
		}
		if(length(object@outlierModules)>0){
			cat("\n##### outlier modules #####\n")
			print(object@outlierModules)
		}
	}
)

