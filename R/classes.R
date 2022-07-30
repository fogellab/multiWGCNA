setClass("WGCNA", slots=list(datExpr="data.frame", conditions="data.frame", trait="data.frame", moduleEigengenes="data.frame", outlierModules="vector"))

setMethod("show", "WGCNA", function(object) {
		cat("##### datExpr #####\n")
		print(head(object@datExpr))
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


