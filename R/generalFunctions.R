#' name: Name of WGCNAobject
#'
#' Returns the module name for a given gene
#'
#' @param WGCNAobject an object of class WGCNA
#' @param gene a gene symbol, should be present in WGCNA object
#' @param use which module label to use e.g. dynamicLabels or dynamicColors
#' 
#' @return A module name
#'
#' @export
#' 
#' @examples 
#' library(ExperimentHub)
#' eh = ExperimentHub()
#' eh_query = query(eh, c("multiWGCNAdata"))
#' astrocyte_networks = eh_query[["EH8222"]]
#' getModule(astrocyte_networks[[1]], 'Srgn')
getModule <- function(WGCNAobject, gene, use = 'dynamicLabels'){
  if(!gene %in% WGCNAobject@datExpr$X) stop('Gene not found!')
  WGCNAobject@datExpr[WGCNAobject@datExpr$X == gene, use]
}

printModules <- function(WGCNAobject){
	name=name(WGCNAobject)
	dir_name=paste0(name, "_Genes")
	dir.create(dir_name)
	datExpr=WGCNAobject@datExpr
	for(module in sort(unique(datExpr$dynamicLabels))){
		file_name=paste0(dir_name, "/", module, ".txt")
		message("### writing ", file_name, " ###\n")
		write.table(datExpr$X[datExpr$dynamicLabels==module], 
			file_name, 
			quote=FALSE, 
			row.names=FALSE, 
			col.names=FALSE)
	}
}

#' name: Name of WGCNAobject
#'
#' Returns the name of a WGCNAobject.
#'
#' @param WGCNAobject an object of class WGCNA
#' 
#' @return Returns the name of the WGCNA object, ie "EAE" for 
#' astrocyte_networks$EAE. 
#'
#' @import stringr
#' @export
#' 
#' @examples 
#' library(ExperimentHub)
#' eh = ExperimentHub()
#' eh_query = query(eh, c("multiWGCNAdata"))
#' astrocyte_networks = eh_query[["EH8222"]]
#' name(astrocyte_networks$EAE)
name <- function(WGCNAobject){
	return(str_split_fixed(WGCNAobject@datExpr$dynamicLabels[[1]], "_", 2)[[1]])
}

getLevel <- function(level, design){
	if(level==1){
		return("combined")
	}
	return(unique(design[,level]))
}

colors <- function(nColors, random=FALSE){
	colors=c("cyan", "grey", "blue", "brown", "darkgreen", "darkmagenta", "darkolivegreen", "darkorange", "darkred",
		"darkturquoise", "green", "darkgrey", "greenyellow", "lightgreen", "magenta", "midnightblue", "orange",
		"paleturquoise", "pink", "purple", "red", "black", "royalblue", "saddlebrown", "salmon", "sienna3", "skyblue",
		"tan", "turquoise", "violet", "yellowgreen", "gold", "maroon", "yellow")
	if(random){
		return(sample(colors, nColors))
	} else {
		return(colors[1:nColors])
	}
}

#' topNGenes: Top N genes of a module
#'
#' Returns the top N number of genes of a module.
#' All genes returned if no number is specified.
#' Genes are in order of intramodular connectivity.
#'
#' @param WGCNAobject an object of class WGCNA
#' @param module the name of the module in WGCNAobject
#' @param nGenes an integer from 1 to module size; returns all genes if left NULL
#'
#' @return a character vector of the genes/probes in the module
#' 
#' @import dplyr
#' @export
#' 
#' @examples 
#' library(ExperimentHub)
#' eh = ExperimentHub()
#' eh_query = query(eh, c("multiWGCNAdata"))
#' astrocyte_networks = eh_query[["EH8222"]]
#' topNGenes(astrocyte_networks$EAE, "EAE_015", nGenes = 10)
#' 
topNGenes <- function(WGCNAobject, module, nGenes=NULL){
	datExpr = WGCNAobject@datExpr
	subsetDatExpr = datExpr[datExpr$dynamicLabels==module,]
	orderedDatExpr = subsetDatExpr %>% arrange(-kWithin)
	if(is.null(nGenes)){
		orderedDatExpr$X
	} else{
		head(orderedDatExpr$X, nGenes)
	}
}

subsetDatExpr <- function(WGCNAobject, geneList){
  datExpr= WGCNAobject@datExpr
  return(cleanDatExpr(datExpr[datExpr$X %in% geneList,]))
}

prepareWGCNA <- function(WGCNAobject){
	newWGCNA= findBestTrait(WGCNAobject)
	newWGCNA = findOutlierModules(newWGCNA)
	newWGCNA
}

addWGCNA <- function(FUN){
	FUN <- match.fun(FUN)
	mixedWGCNA=FUN(mixedWGCNA)
	dimensionOneWGCNAs=lapply(dimensionOneWGCNAs, FUN)
	dimensionTwoWGCNAs= lapply(dimensionTwoWGCNAs, FUN)
}

removeOutlierModules <- function(WGCNAobject, outlierModules=NULL){
	if(missing(outlierModules)){
		outlierModules=outlierModules(WGCNAobject)
	}
	WGCNAobject@datExpr=WGCNAobject@datExpr[!WGCNAobject@datExpr$dynamicLabels %in% outlierModules,]
	WGCNAobject@trait=WGCNAobject@trait[!gsub("ME","",WGCNAobject@trait$Module) %in% outlierModules,]
	WGCNAobject@moduleEigengenes=WGCNAobject@moduleEigengenes[!gsub("ME", "", rownames(WGCNAobject@moduleEigengenes)) %in% outlierModules,]
	WGCNAobject
}

#' Generate a trait table from a sample table
#' 
#' Generates a WGCNA-compatible trait table from a sampleTable dataframe
#' 
#' @param inputTable the sampleTable data.frame
#' @param column the column from the sampleTable to use as traits
#' @param detectNumbers whether to consider traits with numbers as numerical rather than categorical variables
#' 
#' @return a data.frame with integer values denoting the categorical sample traits
#' 
#' @export
#' 
#' @examples 
#' sampleTable = data.frame(Sample = c(paste0("EAE", 1:10), paste0("WT", 1:10)), 
#'                          Disease = c(rep("EAE", 10), rep("WT", 10)), 
#'                          Region = c(rep(c("Cbl", "Sc"), 5))) 
#' makeTraitTable(sampleTable, 2)
#' 
makeTraitTable <- function(inputTable, column, detectNumbers=FALSE) {
	traits=unique(inputTable[,column])
	traitTable=list()
	traitTable[[1]]=inputTable[,1]
	element=2
	for(trait in traits){
		traitTable[[element]]=lapply(inputTable[, column], function(x) if(x==trait) {
			1
		} else {
			2
		}
		)
		element=element+1
	}
	if(detectNumbers){
	  if(is.na(parse_number(inputTable[1, column]))){
                	traitTable[[element]]=as.list(parse_number(inputTable[,column]))
                	traitTable=t(do.call(rbind,traitTable))
                	colnames(traitTable)=c("Sample", c(traits),"Cont")
	  }
	} else {
                traitTable=t(do.call(rbind,traitTable))
                colnames(traitTable)=c("Sample", c(traits))
	}
	traitTable=as.data.frame(apply(traitTable, 2, unlist))
	traitTable
}

#' cleanDatExpr
#'
#' A function that converts a data.frame where row 1 is gene symbols to a 
#' numeric matrix where columns are genes and rows are samples for compatibility 
#' with most WGCNA functions. 
#'
#' @param datExpr a data.frame were columns are samples and rows are samples and the gene symbols are in the first row
#' @param checkGenesSamples call the WGCNA function checkGenesSamples?
#'
#' @return Returns a datExpr with rows as samples and columns as genes
#' 
#' @author Dario Tommasini
#'
#' @export
#' 
#' @examples
#' library(ExperimentHub)
#' eh = ExperimentHub()
#' eh_query = query(eh, c("multiWGCNAdata"))
#' astrocyte_se = eh_query[["EH8223"]]
#' datExpr = data.frame(X = rownames(assays(astrocyte_se)[[1]]), assays(astrocyte_se)[[1]])
#' cleanDatExpr(datExpr)
#' 
cleanDatExpr <- function(datExpr, checkGenesSamples=FALSE) {
	cleanDatExpr <- t(datExpr[ ,!colnames(datExpr) %in% c("X","kTotal","kWithin","kOut","kDiff","dynamicColors","dynamicLabels")])
	colnames(cleanDatExpr) = datExpr$X
	if(checkGenesSamples){
		gsg=goodSamplesGenes(cleanDatExpr)
		if(!gsg$allOK) {
			cleanDatExpr = cleanDatExpr[gsg$goodSamples, gsg$goodGenes]
		}
	}
	return(cleanDatExpr)
}

keyModules <- function(WGCNAobject){
	keyModules=WGCNAobject@trait$Module[WGCNAobject@trait$trait!="None"]
	keyModules=keyModules[!keyModules %in% WGCNAobject@outlierModules]
	keyModules=keyModules[!endsWith(keyModules, "0")]
	return(keyModules)
}

#' summarizeResults: Summarize results from a results list object
#'
#' Prints (or writes) a summary of the results from a results list object
#'
#' @param myNetworks a list of WGCNAobjects
#' @param results results list
#' @param alphaLevel alpha level of significance
#' @param write write to file?
#' @param outputFile name of output file, defaults to results.txt
#' 
#' @return prints a summary of results from the multiWGCNA analysis
#' 
#' @export
summarizeResults <- function(myNetworks, results, alphaLevel=0.05, write=FALSE, outputFile="results.txt") {
	if(write) sink(outputFile)

	#identify interesting (trait-associated) modules across levels
	modulesOfInterest=vector()
	for(network in myNetworks){
		modulesOfInterest=append(modulesOfInterest, keyModules(network))
	}

	for(element in results$preservation){
		message("### Non-preserved modules ###")
		mod1Pres=element[[1]][which(element[[1]]$Zsummary.pres < 10 & rownames(element[[1]]) %in% modulesOfInterest),]
		if(nrow(mod1Pres>0))	message(paste0(capture.output(mod1Pres), collapse = "\n"))
		mod2Pres=element[[2]][which(element[[2]]$Zsummary.pres < 10 & rownames(element[[2]]) %in% modulesOfInterest),]
		if(nrow(mod2Pres>0)) message(paste0(capture.output(mod2Pres), collapse = "\n"))
	}

	diffModExp=results$diffModExp
	message("### Differentially expressed modules ###")
	message(paste0(capture.output(diffModExp[apply(diffModExp, 1, function(p) any(p<0.05)),]), collapse = "\n"))
	if(write) sink()
}

#' iterate: Iterate function across networks
#'
#' A high level function that iterates functions across
#' a list of WGCNA objects
#'
#' @param WGCNAlist a vector of objects of type WGCNAobject
#' @param FUN function to iterate, either overlapComparisons or preservationComparisons
#' @param ... argmuents to be passed on to overlapComparisons or preservationComparisons
#' 
#' @return a comparison list from overlapComparisons or preservationComparisons
#'
#' @author Dario Tommasini
#'
#' @export
#' 
#' @examples
#' library(ExperimentHub)
#' eh = ExperimentHub()
#' eh_query = query(eh, c("multiWGCNAdata"))
#' astrocyte_networks = eh_query[["EH8222"]]
#' results = list()
#' iterate(astrocyte_networks, overlapComparisons, plot=FALSE)
#' 
iterate <- function(WGCNAlist, FUN, ...){
  
  # Check arguments
  stopifnot(deparse(substitute(FUN)) %in% c("overlapComparisons", "preservationComparisons"))
  
	comparisonList=list()
	FUN <- match.fun(FUN)
	element=1
	for(first in 1:(length(WGCNAlist)-1)){
		for(second in (first+1):length(WGCNAlist)){
			comparisonList <- FUN(comparisonList,
							WGCNAlist,
							first,
							second,
							element,
							...)
			element=element+1
		}
	}
	return(comparisonList)
}

#' @importFrom grDevices dev.off pdf rgb
#' @importFrom stats IQR anova dist lm na.omit p.adjust phyper quantile runif sd t.test var
#' @importFrom utils capture.output head write.csv write.table
#' @importFrom graphics legend
NULL
