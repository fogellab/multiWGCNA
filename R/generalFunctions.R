printModules <- function(WGCNAobject){
	name=name(WGCNAobject)
	dir_name=paste0(name, "_Genes")
	dir.create(dir_name)
	datExpr=WGCNAobject@datExpr
	for(module in sort(unique(datExpr$dynamicLabels))){
		file_name=paste0(dir_name, "/", module, ".txt")
		cat("### writing ", file_name, " ###\n")
		write.table(datExpr$X[datExpr$dynamicLabels==module], 
			file_name, 
			quote=F, 
			row.names=F, 
			col.names=F)
	}
}

#' name: Name of WGCNAobject
#'
#' Returns the name of a WGCNAobject.
#'
#' @param WGCNAobject an object of type WGCNAobject
#'
#' @import stringr
#' @export
name <- function(WGCNAobject){
	return(str_split_fixed(WGCNAobject@datExpr$dynamicLabels[[1]], "_", 2)[[1]])
}

getLevel <- function(level, design=sampleTable){
	if(level==1){
		return("combined")
	}
	return(unique(design[,level]))
}

colors <- function(nColors, random=F){
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
#' @param WGCNAobject an object of type WGCNAobject
#' @param module the name of the module in WGCNAobject
#' @param nGenes an integer from 1 to module size; returns all genes if left NULL
#'
#' @import dplyr
#' @export
topNGenes <- function(WGCNAobject, module, nGenes=NULL){
	datExpr= WGCNAobject@datExpr
	subsetDatExpr=datExpr[datExpr$dynamicLabels==module,]
	orderedDatExpr=subsetDatExpr %>% arrange(-kWithin)
	if(missing(nGenes)){
		orderedDatExpr$X
	} else{
		head(orderedDatExpr$X, geneList)
	}
}

subsetDatExpr <- function(WGCNAobject, geneList){
  datExpr= WGCNAobject@datExpr
  return(cleanDatExpr(datExpr[datExpr$X %in% geneList,]))
}

prepareWGCNA <- function(WGCNAobject){
	newWGCNA= findBestTrait(WGCNAobject)
	#newWGCNA = findModuleEigengenes(newWGCNA)
	newWGCNA = findOutlierModules(newWGCNA)
	newWGCNA
}

addWGCNA <- function(FUN){
	FUN <- match.fun(FUN)
	mixedWGCNA=FUN(mixedWGCNA)
	dimensionOneWGCNAs=lapply(dimensionOneWGCNAs, FUN)
	dimensionTwoWGCNAs= lapply(dimensionTwoWGCNAs, FUN)
}

loadMultiWGCNA <- function(expressionFile, traitFile, softPower, alphaLevel){
	sourceFiles()
	loadDependencies()
	datExpr <<- read.csv(expressionFile, header=T)
	sampleTable <<- read.csv(traitFile, header=T)
	alphaLevel <<- alphaLevel
	firstDimension <<- unique(sampleTable[,2])
	secondDimension <<- unique(sampleTable[,3])
	mixedWGCNA=new("WGCNA",
				datExpr=read.csv("Mixed/ModuleSummary/Mixed_datExpr2_ksummary_dynamiccolors_dynamiclabels.csv"),
				trait=read.csv("Mixed/Traits/Mixed_heatmap_modulelabelsTraitCor.csv"))

	dimensionOneWGCNAs=list()
	for(element in 1:length(firstDimension)){
		dimensionOneWGCNAs[[element]]=new("WGCNA",
			datExpr=read.csv(paste0(firstDimension[[element]],
				"/ModuleSummary/",
				firstDimension[[element]],
				"_datExpr2_ksummary_dynamiccolors_dynamiclabels.csv"),
				header=T),
			trait=read.csv(paste0(firstDimension[[element]],
				"/Traits/",
				firstDimension[[element]],
				"_heatmap_modulelabelsTraitCor.csv"),
				header=T))
	}
	names(dimensionOneWGCNAs)=firstDimension

	dimensionTwoWGCNAs=list()
	for(element in 1:length(secondDimension)){
		dimensionTwoWGCNAs[[element]]=new("WGCNA",
			datExpr=read.csv(paste0(secondDimension[[element]],
				"/ModuleSummary/",
				secondDimension[[element]],
				"_datExpr2_ksummary_dynamiccolors_dynamiclabels.csv"),
				header=T),
			trait=read.csv(paste0(secondDimension[[element]],
				"/Traits/",
				secondDimension[[element]],
				"_heatmap_modulelabelsTraitCor.csv"),
				header=T))
	}
	names(dimensionTwoWGCNAs)=secondDimension

	mixedWGCNA <- prepareWGCNA(mixedWGCNA)
	dimensionOneWGCNAs <- lapply(dimensionOneWGCNAs, prepareWGCNA)
	dimensionTwoWGCNAs <- lapply(dimensionTwoWGCNAs, prepareWGCNA)

	mixedWGCNA <<- mixedWGCNA
	dimensionOneWGCNAs <<- dimensionOneWGCNAs
	dimensionTwoWGCNAs <<- dimensionTwoWGCNAs

	allWGCNAlist <<- c(mix=mixedWGCNA, dimensionOneWGCNAs, dimensionTwoWGCNAs)

	firstDimensionComparisons=list()
	firstDimensionComparisons <<- iterate(firstDimensionComparisons, dimensionOneWGCNAs, overlapComparisons, plot=F)
	secondDimensionComparisons=list()
	secondDimensionComparisons <<- iterate(secondDimensionComparisons, dimensionTwoWGCNAs, overlapComparisons, plot=F)
	allOverlapComparisons=list()
	allOverlapComparisons <<- iterate(allOverlapComparisons, allWGCNAlist, overlapComparisons, plot=F)
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


#' generate a trait table from a design table
#' @export
makeTraitTable <- function(inputTable, column, detectNumbers=get("detectNumbers", envir = parent.frame())) {
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
#' @author Dario Tommasini
#'
#' @examples
#'
#' datExpr=data.frame(c("A", "B", "C"), matrix(sample.int(100, size=12, replace=TRUE),nrow=3,ncol=4), c("mod_1", "mod_2", "mod_1"))
#' colnames(datExpr)=c("X", "sample.1", "sample.2", "sample.3", "sample.4", "dynamicLabels")
#' cleanDatExpr(datExpr)
#'
#' @export
cleanDatExpr <- function(datExpr, checkGenesSamples=F) {
	cleanDatExpr <- t(datExpr[ ,!colnames(datExpr) %in% c("X","kTotal","kWithin","kOut","kDiff","dynamicColors","dynamicLabels")])
	colnames(cleanDatExpr) = datExpr$X
	if(checkGenesSamples){
		gsg=goodSamplesGenes(cleanDatExpr)
		if(!gsg$allOK) {
			cleanDatExpr = cleanDatExpr[gsg$goodSamples, gsg$goodGenes]
		}
	}
	cleanDatExpr
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
#' @export
summarizeResults <- function(myNetworks, results, alphaLevel=get("alphaLevel", envir = parent.frame()), write=FALSE, outputFile="results.txt") {
	if(write) sink(outputFile)
	#name1=str_split_fixed(rownames(comparison$mod1Preservation[rownames(comparison$mod1Preservation) != "gold",]),"_",2)[,1][[1]]
	#name2=str_split_fixed(rownames(comparison$mod2Preservation[rownames(comparison$mod2Preservation) != "gold",]),"_",2)[,1][[1]]

	#identify interesting (trait-associated) modules across levels
	modulesOfInterest=vector()
	for(network in myNetworks){
		modulesOfInterest=append(modulesOfInterest, keyModules(network))
	}

	for(element in results$preservation){
		cat("\n### Non-preserved modules ###\n")
		mod1Pres=element[[1]][which(element[[1]]$Zsummary.pres < 10 & rownames(element[[1]]) %in% modulesOfInterest),]
		if(nrow(mod1Pres>0)) print(mod1Pres)
		mod2Pres=element[[2]][which(element[[2]]$Zsummary.pres < 10 & rownames(element[[2]]) %in% modulesOfInterest),]
		if(nrow(mod2Pres>0)) print(mod2Pres)
	}

	diffModExp=results$diffModExp
	#diffModExp$isInterestingModule = input1$log10Pvalue[match(paste0("ME", diffModExp$mod1), input1$Module)] > -log10(alpha) |
	#							input2$log10Pvalue[match(paste0("ME", diffModExp$mod2),input2$Module)] > -log10(alpha)
	cat(paste0("\n### Differentially expressed modules ###\n"))
	#print(diffModExp[rownames(diffModExp) %in% modulesOfInterest & apply(diffModExp, 1, function(p) any(p<0.05)), ])
	#print(modulesOfInterest)
	print(diffModExp[apply(diffModExp, 1, function(p) any(p<0.05)),])
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
#'
#' @author Dario Tommasini
#'
#' @export
iterate <- function(WGCNAlist, FUN, ...){
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
	comparisonList
}


