name <- function(WGCNAobject){
	return(str_split_fixed(WGCNAobject@datExpr$dynamicLabels[[1]], "_", 2)[[1]])
}

loadWGCNAfromDirectory <- function(directory){

	datExpr_file=paste0(directory, grep("datExpr.csv", list.files(directory, recursive=T), value=T))
	conditions_file=paste0(directory, grep("conditions.csv", list.files(directory, recursive=T), value=T))
	ME_file=paste0(directory, grep("moduleEigengenes.csv", list.files(directory, recursive=T), value=T))
	trait_file=paste0(directory, grep("TraitCor", list.files(directory, recursive=T), value=T))


	cat("datExpr file: ", datExpr_file, "\n")
	cat("conditions file: ", conditions_file, "\n")
	cat("module eigengenes file: ", ME_file, "\n")
	cat("trait file: ", trait_file, "\n")

	return(new("WGCNA", datExpr=read.csv(datExpr_file),
		conditions=read.csv(conditions_file),
		moduleEigengenes=read.csv(ME_file),
		trait=read.csv(trait_file)))
}

loadWGCNAfromDirectoryOld <- function(directory){

	datExpr_file=paste0(directory, grep("datExpr2_ksummary_dynamiccolors_dynamiclabels.csv", list.files(directory, recursive=T), value=T))
	conditions_file=paste0(directory, grep("Trait_edited.csv", list.files(directory, recursive=T), value=T))
	trait_file=paste0(directory, grep("modulelabelsTraitCor.csv", list.files(directory, recursive=T), value=T))

	cat("datExpr file: ", datExpr_file, "\n")
	cat("conditions file: ",conditions_file, "\n")
	cat("trait file: ",trait_file, "\n")

	return(new("WGCNA", datExpr=read.csv(datExpr_file),
		conditions=read.csv(conditions_file),
		trait=read.csv(trait_file)))
}

getLevel <- function(level, design=sampleTable){
	if(level==1){
		return("combined")
	}
	return(unique(design[,level]))
}

colors <- function(nColors, random=F){
	colors=c("grey", "cyan", "blue", "brown", "darkgreen", "darkmagenta", "darkolivegreen", "darkorange", "darkred",
		"darkturquoise", "green", "darkgrey", "greenyellow", "lightgreen", "magenta", "midnightblue", "orange",
		"paleturquoise", "pink", "purple", "red", "black", "royalblue", "saddlebrown", "salmon", "sienna3", "skyblue",
		"tan", "turquoise", "violet", "yellowgreen", "gold", "maroon", "yellow")
	if(random){
		return(sample(colors, nColors))
	} else {
		return(colors[1:nColors])
	}
}

#' Top N genes of a module
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
		head(orderedDatExpr$X, nGenes)
	}
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


loadDependencies <- function(){
	library(WGCNA)
	library(ggrepel)
	library(ggpubr)
	library(stringr)
	library(dplyr)
	library(reshape2)
	library(RColorBrewer)
	library(tidyverse)
	library(data.table)
	library(ggalluvial)
	library(readr)
	library(psych)
	library(patchwork)
	library(doParallel)
	library(scales)
	library(igraph)
	library(flashClust)
	library(filesstrings)
	library(biomaRt)
	library(dcanr)
	library(limma)
}

sourceFiles <- function(){
	source("~/projects/novelNetworkApproaches/generalFunctions.R")
	source("~/projects/novelNetworkApproaches/runWGCNA.R")
	source("~/projects/novelNetworkApproaches/overlaps.R")
	source("~/projects/novelNetworkApproaches/preservation.R")
	source("~/projects/novelNetworkApproaches/diffModuleExpAnalysis.R")
	source("~/projects/novelNetworkApproaches/networks.R")
	source("~/projects/novelNetworkApproaches/classes.R")
	source("~/projects/novelNetworkApproaches/geneOntology.R")

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

#' Summarize results from a results list object
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

#' Iterate function across networks
#'
#' A high level function that iterates functions across
#' a list of WGCNA objects
#'
#' @param WGCNAlist a vector of objects of type WGCNAobject
#' @param FUN function to iterate, either overlapComparisons or preservationComparisons
#' @param plot generate plots?
#'
#' @author Dario Tommasini
#'
#' @examples
#'
#' iterate(myNetworks, overlapComparisons, plot=TRUE)
#'
#' @export
iterate <- function(WGCNAlist, FUN, plot=FALSE){
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
							plot=plot)
			element=element+1
		}
	}
	comparisonList
}


