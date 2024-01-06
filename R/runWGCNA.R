performWGCNA <- function(datExpr, traitData, identifier, alphaLevel=0.05, write=FALSE, plot=FALSE, ...){
  arguments=list(...)
  
  # Set network type to WGCNA's default unsigned networkType if not defined in constructNetworks function
  if(is.null(arguments$networkType)) arguments$networkType = "unsigned"
  
  datExpr = t(cleanDatExpr(datExpr, checkGenesSamples = TRUE))
  my_net = blockwiseModules(t(datExpr), ...)
  degrees1=intramodularConnectivity.fromExpr(t(datExpr), my_net$colors,
                                            networkType=arguments$networkType, power=arguments$power)
  dynamicLabels=paste(identifier, "_", str_pad(my_net$colors, 3, pad="0"), sep="")
  summary = cbind(data.frame(X = rownames(datExpr), datExpr), degrees1, dynamicLabels)
  if(write) write.csv(summary, file=paste0(identifier, "_summary.csv"), row.names=FALSE)
  myWGCNA <- new("WGCNA", datExpr=summary, conditions=traitData)
	myWGCNA=findModuleEigengenes(myWGCNA, write=write)
	myWGCNA=traitCor(myWGCNA, write=write)
	myWGCNA=findBestTrait(myWGCNA, alphaLevel = alphaLevel)
	myWGCNA=findOutlierModules(myWGCNA)
	if(plot) plotModules(myWGCNA, mode="PC1")
	return(myWGCNA)
}

#returns the indices of modules that are driven by single samples (outliers)
findOutlierModules <- function(WGCNAobject, byName=TRUE, method="Var", varCutoff=0.01, IQRcutoff=4, zScoreCutoff=4){
        modules=gsub("ME", "", rownames(WGCNAobject@moduleEigengenes))
        if(method=="zScore"){
                ME.zScores=zScoreMatrix(WGCNAobject@moduleEigengenes)
                if(byName){
                        outlierModules=modules[apply(ME.zScores, 1, function(x) any(abs(x) > zScoreCutoff))]
                } else {
                        outlierModules=which(apply(ME.zScores, 1, function(x) any(abs(x) > zScoreCutoff)))
                }
        }
        if(method=="IQR"){
                if(byName){
                        outlierModules=modules[apply(WGCNAobject@moduleEigengenes, 1, function(MEs){
                                any(MEs < quantile(MEs)[[2]]-IQRcutoff*IQR(MEs)) | any(MEs > quantile(MEs)[[4]]+IQRcutoff*IQR(MEs))
                        })]
                } else {

                }
        }
        if(method=="Var"){
                varWithMax=apply(WGCNAobject@moduleEigengenes, 1, range)
                varWithoutMax=apply(WGCNAobject@moduleEigengenes, 1, function(x) {
                        x=x[-which.max(abs(x))]; var(x)
                        })
                outlierModules=modules[varWithoutMax<varCutoff]
        }
        WGCNAobject@outlierModules=outlierModules
        WGCNAobject
}

findBestTrait <- function(WGCNAobject, alphaLevel=0.05, p.adjust=FALSE, write=FALSE) {
        traitTable=WGCNAobject@trait
        group=apply(traitTable[,which(startsWith(colnames(traitTable),"p.value"))], 1, which.min)
        traitTable$trait=gsub("p.value.", "", colnames(traitTable)[which(startsWith(colnames(traitTable), "p.value"))])[group]
        bestPvalues=apply(traitTable[,which(startsWith(colnames(traitTable),"p.value"))], 1, function(x) x[[which.min(x)]])
	      bestPvalues=as.numeric(bestPvalues)
	      traitTable$log10Pvalue= -log10(bestPvalues)
        traitTable$trait[traitTable$log10Pvalue<(-log10(as.numeric(alphaLevel)))]="None"
        traitTable$Module=gsub("ME", "", traitTable$Module)
        WGCNAobject@trait=traitTable
        return(WGCNAobject)
}

traitCor <- function(WGCNAobject, write=FALSE){
	moduleEigengenes=t(WGCNAobject@moduleEigengenes)
	datExpr2=cleanDatExpr(WGCNAobject@datExpr)
	traitData=WGCNAobject@conditions
        identifier=name(WGCNAobject)
	nSamples=nrow(datExpr2)
	Traits=traitData[match(rownames(datExpr2), traitData$Sample),-1]
        rownames(Traits)=traitData[match(rownames(datExpr2), traitData$Sample),c("Sample")]
        datTraits=as.data.frame(Traits)
        moduleTraitCorL = cor(moduleEigengenes, datTraits, use = "p");
        moduleTraitPvalueL = corPvalueStudent(moduleTraitCorL, nSamples);
        colnames(moduleTraitPvalueL) = paste0("p.value.", colnames(moduleTraitCorL));
        traitCor=cbind(Module=gsub("ME", "", rownames(moduleTraitCorL)), moduleTraitCorL, moduleTraitPvalueL)
        rownames(traitCor)=seq_len(nrow(traitCor))
        if(write) write.csv(traitData, paste0(identifier,"_conditions.csv"), row.names=FALSE)
	      if(write) write.csv(traitCor, paste0(identifier,"_TraitCor.csv"), row.names=FALSE)
        WGCNAobject@trait=as.data.frame(traitCor)
        return(WGCNAobject)
}

findModuleEigengenes <- function(WGCNAobject, write=FALSE){
        datExpr2=WGCNAobject@datExpr
        identifier=name(WGCNAobject)
	moduleEigengenes=moduleEigengenes(cleanDatExpr(datExpr2), colors = datExpr2$dynamicLabels, nPC=1)$eigengenes
        moduleEigengenes=as.data.frame(t(moduleEigengenes))
	rownames(moduleEigengenes)=gsub("ME", "", rownames(moduleEigengenes))
	if(write) write.csv(moduleEigengenes, paste0(identifier,"_moduleEigengenes.csv"), row.names=TRUE)
	WGCNAobject@moduleEigengenes=moduleEigengenes
	return(WGCNAobject)
}

plotModules <- function(WGCNAobject, mode="PC1"){
	datExpr2=WGCNAobject@datExpr
	traitData=WGCNAobject@conditions
	identifier=name(WGCNAobject)
	pdf(paste0(identifier, "_moduleExpression.pdf"))
        for(module in sort(unique(datExpr2$dynamicLabels))){
                moduleGenes=datExpr2$X[datExpr2$dynamicLabels==module]
                print(moduleExpressionPlot(WGCNAobject, moduleGenes, mode=mode, title=module))
        }
        dev.off()
}

#' constructNetworks: Construct all the weighted gene correlation networks
#'
#' A high level function that returns all networks
#' possible for a given experimental design
#'
#' @param datExpr either a SummarizedExperiment object or data.frame with genes are rows and samples as columns
#' @param sampleTable data.frame with sample names in first column and sample traits in the second and third column. First column should be called "Sample"
#' @param conditions1 first design conditions, ie healthy/disease
#' @param conditions2 second design conditions, ie frontal lobe/temporal lobe
#' @param write write results out to files? 
#' @param alphaLevel significance value passed to findBestTrait function, default is 0.05 
#' @param plot plot modules? Default is false
#' @param ... Arguments to pass to blockwiseModules function 
#'
#' @return A list of WGCNA objects, ie level one, two, and three networks. 
#' 
#' @author Dario Tommasini
#'
#' @import stringr
#' @import readr
#' @import WGCNA
#' @import flashClust
#' @import SummarizedExperiment
#' @export
#' 
#' @examples
#' library(ExperimentHub)
#' eh = ExperimentHub()
#' eh_query = query(eh, c("multiWGCNAdata"))
#' autism_se = eh_query[["EH8219"]]
#' set.seed(1)
#' autism_se = autism_se[sample(rownames(autism_se), 500),]
#' sampleTable = colData(autism_se)
#' conditions1 = unique(sampleTable[,2])
#' conditions2 = unique(sampleTable[,3])
#' autism_networks = constructNetworks(autism_se, sampleTable, conditions1[[1]], conditions2[[1]], 
#'   networkType = "signed", TOMType = "unsigned", 
#'   power = 10, minModuleSize = 100, maxBlockSize = 25000,
#'   reassignThreshold = 0, minKMEtoStay = 0, mergeCutHeight = 0,
#'   numericLabels = TRUE, pamRespectsDendro = FALSE, 
#'   deepSplit = 4, verbose = 3)
#' autism_networks[["combined"]]
#' 
constructNetworks <- function(datExpr, sampleTable, conditions1, conditions2, write=FALSE, alphaLevel=0.05, plot=FALSE, ...){

  # Check input data format
  stopifnot(inherits(datExpr, "SummarizedExperiment") | inherits(datExpr, "data.frame"))
  
  # Put data in expected format, with genes in first column names "X"
  if(inherits(datExpr, "SummarizedExperiment")){
    datExpr = data.frame(X = rownames(assays(datExpr)[[1]]), assays(datExpr)[[1]])
  } else if(inherits(datExpr, "data.frame")){
    datExpr = data.frame(X = rownames(datExpr), datExpr)
  }
  
	conditions1TraitTable=makeTraitTable(sampleTable, 3) #subset by conditions1, resolve conditions2
	conditions2TraitTable=makeTraitTable(sampleTable, 2) #subset by conditions2, resolve conditions1
	combinedTraitTable=cbind(conditions1TraitTable, conditions2TraitTable[,-1])

	myNetworks=list()
	if(write) {
	  dir.create("combined")
	  setwd("combined")
	}
	myNetworks = append(myNetworks, performWGCNA(datExpr, combinedTraitTable, "combined", alphaLevel = alphaLevel, plot = plot, ...))
	if(write) setwd("..")
	
	# first dimension
	for(trait in unique(conditions1)){
	  if(write) {
	    dir.create(trait)
		  setwd(trait)
	  }
		myNetworks=append(myNetworks, performWGCNA(datExpr[,c("X", sampleTable$Sample[sampleTable[,2]==trait])], conditions1TraitTable[sampleTable[,2]==trait,], trait, plot=plot, ...))
		if(write) setwd("..")
	}

	# second dimension
	for(trait in unique(conditions2)){
	  if(write) {
	    dir.create(trait)
		  setwd(trait)
	  }
		myNetworks=append(myNetworks, performWGCNA(datExpr[,c("X", sampleTable$Sample[sampleTable[,3]==trait])], conditions2TraitTable[sampleTable[,3]==trait,], trait, plot=plot, ...))
		if(write) setwd("..")
	}

	names(myNetworks)=c("combined", conditions1, conditions2)
	return(myNetworks)
}
