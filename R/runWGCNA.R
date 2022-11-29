performWGCNA <- function(datExpr, traitData, identifier, alphaLevel=0.05, write=FALSE, plot=FALSE, ...){
  
  arguments=list(...)
  my_net=blockwiseModules(cleanDatExpr(datExpr), ...)
  degrees1=intramodularConnectivity.fromExpr(cleanDatExpr(datExpr), my_net$colors,
                                            networkType=arguments$TOMType, power=arguments$power)
  dynamicLabels=paste(identifier, "_", str_pad(my_net$colors, 3, pad="0"), sep="")
  summary = cbind(datExpr, degrees1, dynamicLabels)
  if(write) write.csv(summary, file=paste0(identifier, "_summary.csv"), row.names=F)
  myWGCNA <- new("WGCNA", datExpr=summary, conditions=traitData)
	myWGCNA=findModuleEigengenes(myWGCNA, write=write)
	myWGCNA=traitCor(myWGCNA, write=write)
	myWGCNA=findBestTrait(myWGCNA)
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
                        #rownames(WGCNAobject@moduleEigengenes)[outlierModules]
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
        #identifier=name(WGCNAobject)
        group=apply(traitTable[,which(startsWith(colnames(traitTable),"p.value"))], 1, which.min)
        traitTable$trait=gsub("p.value.", "", colnames(traitTable)[which(startsWith(colnames(traitTable), "p.value"))])[group]
        bestPvalues=apply(traitTable[,which(startsWith(colnames(traitTable),"p.value"))], 1, function(x) x[[which.min(x)]])
	      bestPvalues=as.numeric(bestPvalues)
	      traitTable$log10Pvalue= -log10(bestPvalues)
        traitTable$trait[traitTable$log10Pvalue<(-log10(as.numeric(alphaLevel)))]="None"
        traitTable$Module=gsub("ME", "", traitTable$Module)
        #if(write) write.csv(traitTable, paste0(identifier, "_modTraitCor.csv", row.names=F))
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
        #rownames(traitCor)=rownames(moduleTraitCorL)
        rownames(traitCor)=seq_len(nrow(traitCor))
        if(write) write.csv(traitData, paste0(identifier,"_conditions.csv"), row.names=F)
	      if(write) write.csv(traitCor, paste0(identifier,"_TraitCor.csv"), row.names=F)
        WGCNAobject@trait=as.data.frame(traitCor)
        return(WGCNAobject)
}

findModuleEigengenes <- function(WGCNAobject, write=FALSE){
        datExpr2=WGCNAobject@datExpr
        identifier=name(WGCNAobject)
	moduleEigengenes=moduleEigengenes(cleanDatExpr(datExpr2), colors = datExpr2$dynamicLabels, nPC=1)$eigengenes
        moduleEigengenes=as.data.frame(t(moduleEigengenes))
	rownames(moduleEigengenes)=gsub("ME", "", rownames(moduleEigengenes))
	if(write) write.csv(moduleEigengenes, paste0(identifier,"_moduleEigengenes.csv"), row.names=T)
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

#' Construct all the weighted gene correlation networks
#'
#' A high level function that returns all networks
#' possible for a given experimental design
#'
#' @param datExpr data.frame where column 1 has samples and row 1 has genes
#' @param sampleTable data.frame with sample traits
#' @param conditions1 first design conditions, ie healthy/disease
#' @param conditions2 second design conditions, ie frontal lobe/temporal lobe
#' @param softPower thresholding power used for network construction
#' @param originalDirectory directory to output files, defaults to getwd()
#' @param blockWise use the blockwiseModules function? default is TRUE
#' @param nCores number of cores to use for multithreaded network construction
#'
#' @author Dario Tommasini
#'
#' @import stringr
#' @import readr
#' @import WGCNA
#' @import flashClust
#' @export
constructNetworks <- function(datExpr, sampleTable, conditions1, conditions2, write=FALSE, alphaLevel=0.05, plot=FALSE, ...){

	conditions1TraitTable=makeTraitTable(sampleTable, 3) #subset by conditions1, resolve conditions2
	conditions2TraitTable=makeTraitTable(sampleTable, 2) #subset by conditions2, resolve conditions1
	combinedTraitTable=cbind(conditions1TraitTable, conditions2TraitTable[,-1])

	myNetworks=list()
	if(write) {
	  dir.create("combined")
	  setwd("combined")
	}
	myNetworks=append(myNetworks, performWGCNA(datExpr, combinedTraitTable, "combined", alphaLevel=alphaLevel, plot=plot, ...))
	if(write) setwd("..")

	#first dimension
	for(trait in unique(conditions1)){
	  if(write) {
	    dir.create(trait)
		  setwd(trait)
	  }
		myNetworks=append(myNetworks, performWGCNA(datExpr[,c("X", sampleTable$Sample[sampleTable[,2]==trait])], conditions1TraitTable[sampleTable[,2]==trait,], trait, plot=plot, ...))
		if(write) setwd("..")
	}

	#second dimension
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


runWGCNA_minimal <- function(expressionMatrix, identifier, traitData=as.data.frame(NULL), minModuleSize=100, MEDissThres=0, softPower=12, softPowerCalculation=T, saveTOM=F) {

	options(stringsAsFactors=FALSE)

	data = expressionMatrix
	traitData
	if(length(traitData)>0) data = data[,c("X", traitData$Sample)]
	data_1=data[,-1]
	rownames(data_1)=data$X
	data=data_1
	datExpr=as.data.frame(t(data))
	colnames(datExpr)=rownames(data)
	rownames(datExpr)=colnames(data)

	gsg = goodSamplesGenes(datExpr, verbose = 3)
	gsg$allOK

	dim(datExpr)
	if (!gsg$allOK){
		datExpr = datExpr[gsg$goodSamples, gsg$goodGenes]
	}

	datExpr2 = datExpr
	nGenes = ncol(datExpr2);
	nSamples = nrow(datExpr2);

	if(softPowerCalculation){
	pdf(file=paste0("2-", identifier, "_power_dynamic.pdf"),height=10,width=18)
	powers = c(c(1:10), seq(from = 12, to=40, by=2))
	sft = pickSoftThreshold(datExpr2, RsquaredCut = 0.92, powerVector = powers, verbose = 5)
	par(mfrow = c(1,2));
	cex1 = 0.9;
	plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
		xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
		main = paste("Scale independence"));
	text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
		labels=powers,cex=cex1,col="red");
	abline(h=0.90,col="red")
	plot(sft$fitIndices[,1], sft$fitIndices[,5],
		xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
		main = paste("Mean connectivity"))
	text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
	dev.off()
	}
	softPower=softPower

	adjacency = adjacency(datExpr2, power=softPower, type="signed");
	TOM = TOMsimilarity(adjacency);
	dissTOM = 1-TOM
	if(saveTOM) write.csv(dissTOM, paste0("../", identifier, "_dissTOM.csv"))

	geneTree = flashClust(as.dist(dissTOM), method="average");

	dynamicMods= cutreeDynamic(dendro=geneTree, distM=dissTOM, deepSplit=4, pamRespectsDendro=FALSE, minClusterSize=minModuleSize, verbose=2);
	table(dynamicMods)
	dynamicColors = labels2colors(dynamicMods)
	table(dynamicColors)
	library(stringr)
	dynamicLabels <- paste(identifier, "_", str_pad(dynamicMods, 3, pad="0"), sep="")

	mergeC = mergeCloseModules(datExpr2, dynamicColors, cutHeight = MEDissThres, verbose = 3)
	mergeddynamicColors = mergeC$colors
	dynamicColors = mergeddynamicColors

	mergeL = mergeCloseModules(datExpr2, dynamicLabels, cutHeight = MEDissThres, verbose = 3)
	mergeddynamicLabels = mergeL$colors
	dynamicLabels = mergeddynamicLabels

	module_names = setdiff(unique(dynamicColors), "grey")
	gene_names = rownames(t(datExpr2))
	for(label in module_names){
	  module=gene_names[which(dynamicColors==label)]
	}

	module_names = setdiff(unique(dynamicLabels), "grey")
	gene_names = rownames(t(datExpr2))
	for(label in module_names){
	  module=gene_names[which(dynamicLabels==label)]
	}

	degrees1=intramodularConnectivity(adjacency, dynamicColors)
	mydata = t(datExpr2)
	mydata2 = cbind(X=rownames(mydata), mydata, degrees1, dynamicColors, dynamicLabels)
	write.csv(mydata2, file=paste0(identifier, "_summary.csv"), row.names=F)

	myWGCNA <- new("WGCNA", datExpr=mydata2, conditions=traitData)

	return(myWGCNA)
}

runWGCNA_noTraitCor <- function(expressionMatrix, identifier, minModuleSize=100, MEDissThres=0, softPower=12, saveTOM=F) {

	#Step7 Load required libraries
	options(stringsAsFactors=FALSE)

	#Step8 Input normalized data for RNA-seq or expression data normalized for microarray
	data = expressionMatrix
	print(head(data))

	#skip this step if did not import csv file
	data_1=data[,-1]
	rownames(data_1)=data$X
	data=data_1
	head(data[1:3,1:3])

	#Step9 Cluster Samples
	#For rows=samples, and columns=genes do transpose
	datExpr=as.data.frame(t(data))
	#check the data, looks good!
	datExpr[,1:3]

	names(datExpr)=rownames(data)
	rownames(datExpr)=names(data)

	#check for missing data or outlier samples
	gsg = goodSamplesGenes(datExpr, verbose = 3)
	gsg$allOK

	#returned FALSE for this dataset, i.e. there is an outlier gene or sample. Use if statement to remove #outliers
	dim(datExpr)
	if (!gsg$allOK)
	{
	# Remove the offending genes and samples from the data:
	datExpr = datExpr[gsg$goodSamples, gsg$goodGenes]
	}
	#Output = Removed a bunch of genes that had zero variance, missing values, etc
	dim(datExpr)

	#show cluster of sample and visually check for outlier samples
	pdf(file=paste0("1-",identifier,"_sampleClustering_dynamic.pdf"), width = 12, height = 9)
	sampleTree = hclust(dist(datExpr), method = "average");
	par(mar = c(0,4,2,0))
	par(oma=c(1,1,1,1))
	plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,cex.main=2)
	dev.off()

	GSGGrow.no=dim(datExpr)[1]
	GSGGcol.no=dim(datExpr)[2]

	var = apply(datExpr, 2, var)
	dat = rbind(datExpr, var)
	rownames(dat) = c(rownames(datExpr), "variance")
	dat2 = dat[,order(dat["variance",], decreasing=T)] #order columns by variance

	dat2[,1:3]
	dim(dat2)
	dat1 = dat2[1:GSGGrow.no,1:GSGGcol.no] #option2: keep ALL goodsamplegoodgenes
	dim(dat1)
	datExpr2 = dat1

	pdf(file=paste0("2-",identifier,"_power_dynamic.pdf"),height=10,width=18)
	powers = c(c(1:10), seq(from = 12, to=40, by=2))
	sft = pickSoftThreshold(datExpr2, RsquaredCut = 0.92, powerVector = powers, verbose = 5)
	#pickSoftThreshold takes about 15-20min to finish
	par(mfrow = c(1,2));
	cex1 = 0.9;
	#scale-free topology fit index as function of soft-thresholding power
	plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
	xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
	main = paste("Scale independence"));
	text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
	labels=powers,cex=cex1,col="red");
	#this line corresponds to using R^2 cut-off of h
	abline(h=0.90,col="red")
	# Mean connectivity as a function of the soft-thresholding power
	plot(sft$fitIndices[,1], sft$fitIndices[,5],
	xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
	main = paste("Mean connectivity"))
	text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
	#dev.print(pdf, "2-identifier_power_dynamic.pdf",height=10,width=18)
	dev.off()

	#sft$powerEstimate is systems best guess
	if(!is.na(sft$powerEstimate)) {
		softPower=sft$powerEstimate
	} else {
		softPower=sft$fitIndices$Power[which.max(sft$fitIndices$SFT.R.sq)]
	}
	softPower=softPower
	print(softPower)
	write.table(softPower, paste0(identifier, "_softPowerEstimate.txt"))

	#### Variables labeled are called 'dynamic' instead of 'block' ###

	#Step12 Adjacency and TOM

	library(flashClust)
	#calculate adjacency
	adjacency = adjacency(datExpr2, power=softPower, type="signed");
	#Topological overlap matrix (TOM), turn adjacency into topological overlap
	TOM = TOMsimilarity(adjacency);
	dissTOM = 1-TOM
	if(saveTOM) write.csv(dissTOM, paste0("../", identifier, "_dissTOM.csv"))

	#For clustering using TOM,
	#call the hierarchical clustering function
	geneTree = flashClust(as.dist(dissTOM), method="average");

	#plot the resulting tree (dendrogram)
	pdf(file=paste0("3-",identifier,"_Gene clustering_dynamic.pdf"), height=10, width=15)
	plot(geneTree, xlab="", sub="", main="Gene clustering on TOM-based dissimilarity", labels=FALSE, hang=0.04);
	#dev.print(pdf,"3-identifier_Gene clustering_dynamic.pdf", height=10, width=15)
	dev.off()

	#module identification using dynamic tree cut
	dynamicMods= cutreeDynamic(dendro=geneTree, distM=dissTOM, deepSplit=4, pamRespectsDendro=FALSE, minClusterSize=minModuleSize, verbose=2);
	table(dynamicMods)
	#convert numberic labels into colors
	dynamicColors = labels2colors(dynamicMods)
	table(dynamicColors)
	#create labels for modules too
	library(stringr)
	dynamicLabels <- paste(identifier,"_", str_pad(dynamicMods, 3, pad="0"), sep="")
	#the pad number ex 1 becomes 01

	###################################Here onwards the code is again identical to the blockwise module detection, variables labeled with called 'dynamic' instead of 'block'#################

	dim(table(dynamicColors))
	dim(table(dynamicMods))
	dim(table(dynamicLabels))

	head(table(dynamicColors))
	head(table(dynamicMods))
	head(table(dynamicLabels))

	#Plot the dendrogram and colors underneath
	pdf(file=paste0("4-",identifier,"_Gene dendrogram and module colors_dynamic.pdf"), height=10, width=15)
	plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", dendroLabels=FALSE, hang=0.03, addGuide=TRUE, guideHang=0.05, main="Gene dendrogram and module colors")
	#dev.print(pdf,"4-identifier_Gene dendrogram and module colors_dynamic.pdf", height=10, width=15)
	dev.off()

	pdf(file=paste0("5-",identifier,"_Gene dendrogram and module labels_dynamic.pdf"), height=10, width=15)
	plotDendroAndColors(geneTree, dynamicColors,
	c("Dynamic Tree Cut"), dendroLabels = FALSE, hang = 0.03,
	addGuide = TRUE, guideHang = 0.05)
	#dev.print(pdf,"5-identifier_Clustering of module eigengene_dynamic.pdf", height=10, width=15)
	dev.off()

	# Call an automatic merging function for colors
	mergeC = mergeCloseModules(datExpr2, dynamicColors, cutHeight = MEDissThres, verbose = 3)
	mergeddynamicColors = mergeC$colors
	#relabel dynamicColors to the new merged dynamicColors
	dynamicColors = mergeddynamicColors
	pdf(file="6-Feng_Nature_Gene dendrogram and module colors_merged_dynamic.pdf", height=10, width=15)
	plotDendroAndColors(geneTree, dynamicColors,
	c("Dynamic Tree Cut"),dendroLabels = FALSE, hang = 0.03,
	addGuide = TRUE, guideHang = 0.05)
	dev.off()

	# Call an automatic merging function for dynamicLabels
	mergeL = mergeCloseModules(datExpr2, dynamicLabels, cutHeight = MEDissThres, verbose = 3)
	mergeddynamicLabels = mergeL$colors
	#relabel dynamicLabels to the new merged dynamicLabels
	dynamicLabels = mergeddynamicLabels
	pdf(file="6L-Feng_Nature_Gene dendrogram and module labels_merged_dynamic.pdf", height=10, width=15)
	plotDendroAndColors(geneTree, dynamicColors,
	c("Dynamic Tree Cut"),dendroLabels = FALSE, hang = 0.03,
	addGuide = TRUE, guideHang = 0.05)
	dev.off()

	#Step13 Save module size, module gene list, intramodular connectivity and signedkME data

	#WITH COLORS
	#save merged gene list txt file with module colors
	module_names = setdiff(unique(dynamicColors), "grey")
	gene_names = rownames(t(datExpr2))
	for(label in module_names){
	  module=gene_names[which(dynamicColors==label)]
	  write.table(module, paste(identifier,"_", label, ".txt", sep=""), sep="\t", row.names=F, col.names=F, quote=F)
	}

	#WITH LABELS
	#save merged gene list txt file with module labels
	module_names = setdiff(unique(dynamicLabels), "grey")
	gene_names = rownames(t(datExpr2))
	for(label in module_names){
	  module=gene_names[which(dynamicLabels==label)]
	  write.table(module, paste(label, ".txt", sep=""), sep="\t", row.names=F, col.names=F, quote=F)
	}

	#adjacency already calulated in stepwise module detection
	#Summary output of network analysis results for colors in moduleColors
	degrees1=intramodularConnectivity(adjacency, dynamicColors)
	head(degrees1)

	#want to make summary file with genes, genes intramodular connectivity, dynamicColors and dynamicLabels
	mydata = t(datExpr2)
	mydata2 = cbind(mydata, degrees1, dynamicColors, dynamicLabels)
	write.csv(mydata2,file=paste0(identifier,"_datExpr2_ksummary_dynamiccolors_dynamiclabels.csv"))

	#with colors
	# Calculate eigengenes
	MEList = moduleEigengenes(datExpr2, colors = dynamicColors, softPower= softPower, nPC=1)
	MEs = MEList$eigengenes

	pdf(file=paste0("7-",identifier,"_eigengenes-b_dynamic.pdf"),height=9,width=9)
	par(cex = 0.6)
	plotEigengeneNetworks(MEs, "", marDendro = c(0,8,1,2), marHeatmap = c(5,8,1,2), cex.adjacency = 0.3 , cex.preservation = 0.3, plotPreservation = "standard")
	#dev.print(pdf,"7-identifier_eigengenes-b_dynamic.pdf",height=9,width=9)
	dev.off()

	#with labels
	# Calculate eigengenes
	MEListL = moduleEigengenes(datExpr2, colors = dynamicLabels, softPower= softPower, nPC=1)
	MEsL = MEListL$eigengenes

	pdf(file=paste0("7L-",identifier,"_eigengenes-b_dynamic.pdf"),height=9,width=9)
	par(cex = 0.6)
	plotEigengeneNetworks(MEsL, "", marDendro = c(0,8,1,2), marHeatmap = c(5,8,1,2), cex.adjacency = 0.3 , cex.preservation = 0.3, plotPreservation = "standard")
	#dev.print(pdf,"7L-identifier_eigengenes-b_dynamic.pdf",height=9,width=9)
	dev.off()

	#Calculate signedkME. Higher the signed kME for a given gene higher its correlation with eigen gene i.e. its a hub gene for the given module Ref: Miller 2010, Gandal 2018
	geneMM=signedKME(datExpr2,MEs)
	colnames(geneMM)=paste(colnames(geneMM),".cor",sep="")
	MMPvalue=corPvalueStudent(as.matrix(geneMM),dim(datExpr2)[[1]])
	colnames(MMPvalue)=paste(colnames(MMPvalue),".pval",sep="")

	geneMML=signedKME(datExpr2,MEsL)
	colnames(geneMML)=paste(colnames(geneMML),".cor",sep="")
	MMPvalueL=corPvalueStudent(as.matrix(geneMML),dim(datExpr2)[[1]])
	colnames(MMPvalueL)=paste(colnames(MMPvalueL),".pval",sep="")

	#want to make summary file with genes, signedkME, blockColors and blockLabels
	mydata3 = cbind(geneMM, MMPvalue, dynamicColors)
	mydata4 = cbind(geneMML, MMPvalueL, dynamicLabels)
	write.csv(mydata3,file=paste0(identifier,"_datExpr2_signedKME_dynamiccolors.csv"))
	write.csv(mydata4,file=paste0(identifier,"_datExpr2_signedKME_dynamiclabels.csv"))

	write.csv(table(dynamicColors), file=paste0(identifier,"_final_dynamiccolors.csv"))
	write.csv(table(dynamicLabels), file=paste0(identifier,"_final_dynamiclabels.csv"))
	#these make an csv file with the module colors or my labels and how many members are in each

	#Step14 Visualization of final modules, PC1 and eigengenes
	#heatmaps for module colors
	aa = read.csv(paste0(identifier,"_final_dynamiccolors.csv"), header=T, row.names=1)
	MEss = MEs[,sort(colnames(MEs))]
	dim(aa)

	pdf(file=paste0("8-",identifier,"_heatmaps_dynamic.pdf"))
	for (k in 1:length(rownames(aa))) {
		whichmodule = aa$dynamicColors [k]
		par(mfrow=c(2,1))
		par(mar=c(3.5,3,2.5,3))
		par(oma=c(4,0,2,0))
		datcombined=datExpr2[,dynamicColors==whichmodule]
		datcombined=datcombined[order(row.names(datcombined)),]
		#sort so samples in order
		scaledModExpr =scale(datcombined)
		myEx = max(abs(scaledModExpr))
		plotMat(t(scaledModExpr), main=paste(aa$dynamicColors [k], "(", aa$Freq[k], ")", sep=""))
		names=sort(dimnames(datExpr2)[[1]])
		#sort labels
		numSamp = length(rownames(datExpr2))
		datsv=MEss[,k]
		whichmodule=aa$dynamicColors[k]
		space=0.174
		width=(10.65-(space*numSamp))/numSamp
		xpos=seq(0.35+(width/2),11,width+space)
		barplot(datsv,col="grey",xlim=c(0.35,11), width=width, space=space)
		axis(side=1,at=xpos,labels=F,tick=T)
		box()
		mtext(names,1,at=xpos,cex=1,col="black",adj=1.2,las=3)
		print(k)
	}
	dev.off()

	#heatmaps for module LABELS
	aaL = read.csv(paste0(identifier,"_final_dynamiclabels.csv"), header=T, row.names=1)
	MEssL = MEsL[,sort(colnames(MEsL))]
	dim(aaL)

	pdf(file=paste0("8L-",identifier,"_heatmaps_dynamic.pdf"))
	for (k in 1:length(rownames(aaL))) {
	  	whichmodule = aaL$dynamicLabels[k]
		par(mfrow=c(2,1))
		par(mar=c(3.5,3,2.5,3))
		par(oma=c(4,0,2,0))
		datcombined=datExpr2[,dynamicLabels==whichmodule]
		#datcombined=datcombined[order(row.names(datcombined)),] #sort so samples in order
		scaledModExpr =scale(datcombined)
		#myEx = max(abs(scaledModExpr))
		#plotMat(t(scaledModExpr), main=paste(aaL$dynamicModsLabels[k], "(", aaL$Freq[k], ")", sep=""))
		plotMat(t(scaledModExpr), main=paste(aaL$dynamicLabels[k], "(", aaL$Freq[k], ")", sep=""))
		names=rownames(datcombined)
		numSamp = length(rownames(datExpr2))
		datsv=MEssL[,k]
		whichmodule=aaL$dynamicLabels[k]
		space=0.05
		width=(11-(space*numSamp))/numSamp
		#edited on 20201215 by DT
		#xpos=seq(0.35+(width/2),11,width+space)
		#barplot(datsv,col="grey",xlim=c(0.35,11), width=width, space=space)
		#space in barplot takes that fraction of width variable, so we need to divide by width in order to make space equal space
		barplot(datsv,col="grey",xlim=c(0.41,10.6), space=space/width, width=width, axisnames=TRUE, axis.lty=1, names.arg=names, las=2)
		#axis(side=1,at=xpos,labels=F,tick=T)
		#we include the x axis labels in the barplot using names.arg for the names and axis.lty for the tick marks
		box()
		#mtext(names,1,at=xpos,cex=1,col="black",adj=1.2,las=3)

		print(k)
	}
	dev.off()

	#Note if goodsamplegenes removes samples then Count or data variable sample definition won't hold true so we need to redefine the samples based on datExpr2 here
	samples=rownames(datExpr2)
	samples

	#single plots module COLORS (moduleColors)
	PCvalues=MEs
	PCvalues= PCvalues[,order(names(PCvalues),decreasing=F)]
	modules=colnames(PCvalues)
	colorsMod=gsub("ME","",modules)
	library(stringr)
	mains=paste0(identifier,"_", colorsMod, sep="")
	namesSamp=samples

	pdf(file=paste0("9-",identifier,"_singlePC_dynamic.pdf"),width=10,height=10)
	#this makes 5x5 graphs per page, with margins so labels are cut off the page
	par(mfrow=c(5,5))
	par(oma=c(2,1,2,1))
	for (mod in 1:length(PCvalues)) {
	plot(smooth.spline(PCvalues[,mod]), xlab="", type="n", ylim=c(-0.5, 0.5), ylab="First PC", axes=F)
	abline(h=0)
	axis(1, at=1:length(namesSamp), labels=namesSamp, cex.axis=1, las=3)
	axis(2, at=seq(-0.5, 0.5, 0.5))
	#par here lets you plot a second graph on the same plot
	par(new=TRUE)
	barplot(PCvalues[,mod],col="white", ylim=c(-0.5, 0.5), axes=F,main=mains[mod])
	SecNumber=length(samples)+0.5 ##1.5 to Secondnumber or SecNumber just centers the line alittle between around the samples
	lines(smooth.spline(x=c(1.5:SecNumber),PCvalues[,mod], spar=0.3), col=colorsMod[mod],  lwd=4)
	}
	dev.off()

	#single plots module LABELS (moduledynamicLabels)
	PCvaluesL=MEsL
	PCvaluesL= PCvaluesL[,order(names(PCvaluesL),decreasing=F)]
	modulesL=colnames(PCvaluesL)
	colorsModL=gsub("ME","",modulesL)
	colorsModL2=gsub(paste0(identifier,"_"),"",colorsModL)
	colorsModcolors = labels2colors(colorsModL2)
	mainsL=paste0(identifier,"_", colorsModL2, sep="")
	namesSampL =samples

	pdf(file=paste0("9L-",identifier,"_singlePC_dynamic.pdf"),width=10,height=10)
	#this makes 5x5 graphs per page, with margins so labels are cut off the page
	par(mfrow=c(5,5))
	par(oma=c(2,1,2,1))
	for (mod in 1:length(PCvaluesL)) {
	plot(smooth.spline(PCvaluesL[,mod]), xlab="", type="n", ylim=c(-0.5, 0.5), ylab="First PC", axes=F)
	abline(h=0)
	axis(1, at=1:length(namesSampL), labels=namesSampL, cex.axis=1, las=3)
	axis(2, at=seq(-0.5, 0.5, 0.5))
	#par here lets you plot a second graph on the same plot
	par(new=TRUE)
	barplot(PCvaluesL[,mod],col="white", ylim=c(-0.5, 0.5), axes=F,main=mainsL[mod])
	#keeping the module colors for the line, but using numbers for the main labels
	SecNumber=length(samples)+0.5 ##1.5 to Secondnumber or SecNumber just centers the line alittle between around the samples
	lines(smooth.spline(x=c(1.5:SecNumber),PCvaluesL[,mod], spar=0.3), col=colorsModcolors[mod],  lwd=4)
	}
	dev.off()

	###check if number of modules are too many or to few and if biggest module has too many or too few genes
	#Import new module data we got from our analysis containing genes, VST counts, module and module information summary file from this analysis folder identifier_datExpr2_ksummary_dynamiccolors_dynamiclabels.csv
	dat1=read.csv(paste0(identifier,"_datExpr2_ksummary_dynamiccolors_dynamiclabels.csv"), sep=',')
	head(dat1)

	dat1GeneMod=dat1[,c("X","dynamicColors","dynamicLabels")]
	colnames(dat1GeneMod)=c("Gene","dynamicColors","dynamicLabels")
	dat1GeneMod[1:3,]
	dat1GeneMod_table<-as.data.frame(table(dat1GeneMod$dynamicLabels))
	dat1GeneMod_table<-dat1GeneMod_table[order(-dat1GeneMod_table$Freq),]
	write.csv(dat1GeneMod_table,paste0(identifier,"_final_dynamiccolors_sorted.csv"))
	#This result is same as "identifier_final_dynamiccolors.csv" but sorted

	dat1GeneMod_table2<-as.data.frame(table(dat1GeneMod$dynamicColors))
	dat1GeneMod_table2<-dat1GeneMod_table2[order(-dat1GeneMod_table2$Freq),]
	write.csv(dat1GeneMod_table2,paste0(identifier,"_final_dynamiclabels_sorted.csv"))
	#This result is same as "identifier_final_dynamiclabels.csv" but sorted

	#For dynamicColors
	#This is same as the dat1 file imported above
	dat1=dat1
	#Export genes ranked by net connectivity weight or kIM or IM for each module
	for (i in unique(dat1GeneMod$dynamicColors)){

	  #Select modules
	  modules= i
	  modules

	  #net gene Connectivity or weights is a.k.a. kIM or IM or kWithin
	  dat1_sel=dat1[which(dynamicColors==modules),] #select those rows which match module name
	  dat1_sort=dat1_sel[order(-dat1_sel$kWithin),] #order by stronger to weaker connectivity

	  #net gene Connectivity ranks
	  dat1_sort$rank<- rank(-(dat1_sort$kWithin)) #This gives the highest kWithin the highest rank
	  #export output
	  write.table(dat1_sort, paste0("ForGenesInModule-", modules, "-allGenesInModuleRankedByNetKWithin.txt", sep=""),row.names=F, quote=F)
	}

	#For dynamicLabels
	#This is same as the dat1 file imported above
	dat1=dat1
	#Export genes ranked by net connectivity weight or kIM or IM for each module
	for (i in unique(dat1GeneMod$dynamicLabels)){

	  #Select modules
	  modules= i
	  modules

	  #net gene Connectivity or weights is a.k.a. kIM or IM or kWithin
	  dat1_sel=dat1[which(dynamicLabels==modules),] #select those rows which match module name
	  dat1_sort=dat1_sel[order(-dat1_sel$kWithin),] #order by stronger to weaker connectivity

	  #net gene Connectivity ranks
	  dat1_sort$rank<- rank(-(dat1_sort$kWithin)) #This gives the highest kWithin the highest rank
	  #export output
	  write.table(dat1_sort, paste0("ForGenesInModule-", modules, "-allGenesInModuleRankedByNetKWithin.txt", sep=""),row.names=F, quote=F)
	}

	#Import module names and corrosponding gene names from WGCNA summary file 'identifier_datExpr2_ksummary_dynamiccolors_dynamiclabels.csv' from analysis folder
	dat1=read.csv(paste0(identifier,"_datExpr2_ksummary_dynamiccolors_dynamiclabels.csv"), sep=',')
	dat1GeneMod=dat1[,c("X","dynamicColors","dynamicLabels")]
	colnames(dat1GeneMod)=c("Gene","dynamicColors","dynamicLabels")
	dat1GeneMod[1:3,]
	#For MEDissThresh >0.05 for smaller merged modules, multiple colors can be assigned to the same moduledynamicLabels. So its better to use the final moduledynamicLabels for exporting for cytoscape and VisANT

	#Step17 Organization of files
	library(filesstrings)

	dir.create("GeneLists")
	file.move(list.files(pattern='ForGenesInModule-.*'), "GeneLists")
	file.move(list.files(pattern = paste0(identifier,"_.*.txt")), "GeneLists")

	dir.create("Plots")
	file.move(list.files(pattern = '*.pdf'), "Plots")


	#remove extra files
	file.remove(list.files(pattern = "final_modules.csv"))
	file.remove(list.files(pattern = "final_modules_labels.csv"))

	dir.create("ModuleSummary")
	file.move(list.files(pattern = "_datExpr2*"), "ModuleSummary")
	file.move(list.files(pattern = "_final_dynamiccolors.csv"), "ModuleSummary")
	file.move(list.files(pattern = "_final_dynamiclabels.csv"), "ModuleSummary")
	file.move(list.files(pattern = "_final_dynamiccolors_sorted.csv"), "ModuleSummary")
	file.move(list.files(pattern = "_final_dynamiclabels_sorted.csv"), "ModuleSummary")

	#remove extra files
	file.remove(list.files(pattern="dat1GeneLabel.csv"))
	file.remove(list.files(pattern="tdatExpr.csv"))
	file.remove(list.files(pattern="PaperTablewithExpr.csv"))
	file.remove(list.files(pattern="*.RData")) #optional deletion of .RData

	#Step18 Save version numbers for reproducibility and version control
	sessionInfo()
	toLatex(sessionInfo())
	save.image(file=paste0(identifier,'.RData'))
}



runWGCNA <- function(expressionMatrix, designMatrix, identifier, minModuleSize=100, MEDissThres=0, softPower=12, saveTOM=F, saveRData=F) {
	#define the variables
	traitData=designMatrix

	#Step7 Load required libraries
	options(stringsAsFactors=FALSE)

	#Step8 Input normalized data for RNA-seq or expression data normalized for microarray
	data = expressionMatrix
	data = data[,c("X", traitData$Sample)]
	print(head(data))

	#skip this step if did not import csv file
	data_1=data[,-1]
	rownames(data_1)=data$X
	data=data_1
	head(data[1:3,1:3])

	#Step9 Cluster Samples
	#For rows=samples, and columns=genes do transpose
	datExpr=as.data.frame(t(data))
	#check the data, looks good!
	datExpr[,1:3]

	names(datExpr)=rownames(data)
	rownames(datExpr)=names(data)

	#check for missing data or outlier samples
	gsg = goodSamplesGenes(datExpr, verbose = 3)
	gsg$allOK

	#returned FALSE for this dataset, i.e. there is an outlier gene or sample. Use if statement to remove #outliers
	dim(datExpr)
	if (!gsg$allOK)
	{
	# Remove the offending genes and samples from the data:
	datExpr = datExpr[gsg$goodSamples, gsg$goodGenes]
	}
	#Output = Removed a bunch of genes that had zero variance, missing values, etc
	dim(datExpr)

	#show cluster of sample and visually check for outlier samples
	pdf(file=paste0("1-",identifier,"_sampleClustering_dynamic.pdf"), width = 12, height = 9)
	sampleTree = hclust(dist(datExpr), method = "average");
	par(mar = c(0,4,2,0))
	par(oma=c(1,1,1,1))
	plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,cex.main=2)
	dev.off()

	#Step10 OPTIONAL SELECT ALL GENES OR ALTERNATIVELY SELECT TOP VARIABLE GENES USING ONE OF THE CODES BELOW

	#store row number or number of genes after goodsamplegoodgenes (GSGG)
	GSGGrow.no=dim(datExpr)[1]
	#store column number or number of samples after goodsamplegoodgenes
	GSGGcol.no=dim(datExpr)[2]

	###Select All or top5000 genes###
	var = apply(datExpr, 2, var)
	#head(var)
	dat = rbind(datExpr,var)
	#rownames(dat) = c(samples, "variance")
	rownames(dat) = c(rownames(datExpr), "variance")
	dat2 = dat[,order(dat["variance",], decreasing=T)] #order columns by variance

	dat2[,1:3] #check that variance is in order
	dim(dat2)
	#Use one of the two lines below to select top5000 or all genes respectively
	#dat1 = dat2[1:GSGGrow.no,1:5000] #option1: keep top5000 goodsamplegoodgenes
	dat1 = dat2[1:GSGGrow.no,1:GSGGcol.no] #option2: keep ALL goodsamplegoodgenes
	dim(dat1)
	datExpr2 = dat1

	pdf(file=paste0("2-",identifier,"_power_dynamic.pdf"),height=10,width=18)
	powers = c(c(1:10), seq(from = 12, to=40, by=2))
	sft = pickSoftThreshold(datExpr2, RsquaredCut = 0.92, powerVector = powers, verbose = 5)
	#pickSoftThreshold takes about 15-20min to finish
	par(mfrow = c(1,2));
	cex1 = 0.9;
	#scale-free topology fit index as function of soft-thresholding power
	plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
	xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
	main = paste("Scale independence"));
	text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
	labels=powers,cex=cex1,col="red");
	#this line corresponds to using R^2 cut-off of h
	abline(h=0.90,col="red")
	# Mean connectivity as a function of the soft-thresholding power
	plot(sft$fitIndices[,1], sft$fitIndices[,5],
	xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
	main = paste("Mean connectivity"))
	text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
	#dev.print(pdf, "2-identifier_power_dynamic.pdf",height=10,width=18)
	dev.off()

	#sft$powerEstimate is systems best guess
	if(!is.na(sft$powerEstimate)) {
		softPower=sft$powerEstimate
	} else {
		softPower=sft$fitIndices$Power[which.max(sft$fitIndices$SFT.R.sq)]
	}
	softPower=softPower
	print(softPower)
	write.table(softPower, "softPowerEstimate.txt")

	#### Variables labeled are called 'dynamic' instead of 'block' ###

	#Step12 Adjacency and TOM

	library(flashClust)
	#calculate adjacency
	adjacency = adjacency(datExpr2, power=softPower, type="signed");
	#Topological overlap matrix (TOM), turn adjacency into topological overlap
	TOM = TOMsimilarity(adjacency);
	dissTOM = 1-TOM
	if(saveTOM) write.csv(dissTOM, paste0("../", identifier, "_dissTOM.csv"))

	#For clustering using TOM,
	#call the hierarchical clustering function
	geneTree = flashClust(as.dist(dissTOM), method="average");

	#plot the resulting tree (dendrogram)
	pdf(file=paste0("3-",identifier,"_Gene clustering_dynamic.pdf"), height=10, width=15)
	plot(geneTree, xlab="", sub="", main="Gene clustering on TOM-based dissimilarity", labels=FALSE, hang=0.04);
	#dev.print(pdf,"3-identifier_Gene clustering_dynamic.pdf", height=10, width=15)
	dev.off()

	#module identification using dynamic tree cut
	dynamicMods= cutreeDynamic(dendro=geneTree, distM=dissTOM, deepSplit=4, pamRespectsDendro=FALSE, minClusterSize=minModuleSize, verbose=2);
	table(dynamicMods)
	#convert numberic labels into colors
	dynamicColors = labels2colors(dynamicMods)
	table(dynamicColors)
	#create labels for modules too
	library(stringr)
	dynamicLabels <- paste(identifier,"_", str_pad(dynamicMods, 2, pad="0"), sep="")
	#the pad number ex 1 becomes 01

	###################################Here onwards the code is again identical to the blockwise module detection, variables labeled with called 'dynamic' instead of 'block'#################


	dim(table(dynamicColors))
	dim(table(dynamicMods))
	dim(table(dynamicLabels))

	head(table(dynamicColors))
	head(table(dynamicMods))
	head(table(dynamicLabels))

	#Plot the dendrogram and colors underneath
	pdf(file=paste0("4-",identifier,"_Gene dendrogram and module colors_dynamic.pdf"), height=10, width=15)
	plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", dendroLabels=FALSE, hang=0.03, addGuide=TRUE, guideHang=0.05, main="Gene dendrogram and module colors")
	#dev.print(pdf,"4-identifier_Gene dendrogram and module colors_dynamic.pdf", height=10, width=15)
	dev.off()

	pdf(file=paste0("5-",identifier,"_Gene dendrogram and module labels_dynamic.pdf"), height=10, width=15)
	plotDendroAndColors(geneTree, dynamicColors,
	c("Dynamic Tree Cut"),dendroLabels = FALSE, hang = 0.03,
	addGuide = TRUE, guideHang = 0.05)
	#dev.print(pdf,"5-identifier_Clustering of module eigengene_dynamic.pdf", height=10, width=15)
	dev.off()

	# Call an automatic merging function for colors
	mergeC = mergeCloseModules(datExpr2, dynamicColors, cutHeight = MEDissThres, verbose = 3)
	mergeddynamicColors = mergeC$colors
	#relabel dynamicColors to the new merged dynamicColors
	dynamicColors = mergeddynamicColors
	pdf(file="6-Feng_Nature_Gene dendrogram and module colors_merged_dynamic.pdf", height=10, width=15)
	plotDendroAndColors(geneTree, dynamicColors,
	c("Dynamic Tree Cut"),dendroLabels = FALSE, hang = 0.03,
	addGuide = TRUE, guideHang = 0.05)
	dev.off()

	# Call an automatic merging function for dynamicLabels
	mergeL = mergeCloseModules(datExpr2, dynamicLabels, cutHeight = MEDissThres, verbose = 3)
	mergeddynamicLabels = mergeL$colors
	#relabel dynamicLabels to the new merged dynamicLabels
	dynamicLabels = mergeddynamicLabels
	pdf(file="6L-Feng_Nature_Gene dendrogram and module labels_merged_dynamic.pdf", height=10, width=15)
	plotDendroAndColors(geneTree, dynamicColors,
	c("Dynamic Tree Cut"),dendroLabels = FALSE, hang = 0.03,
	addGuide = TRUE, guideHang = 0.05)
	dev.off()

	#Step13 Save module size, module gene list, intramodular connectivity and signedkME data

	#WITH COLORS
	#save merged gene list txt file with module colors
	module_names = setdiff(unique(dynamicColors), "grey")
	gene_names = rownames(t(datExpr2))
	for(label in module_names){
	  module=gene_names[which(dynamicColors==label)]
	  write.table(module, paste(identifier,"_", label, ".txt", sep=""), sep="\t", row.names=F, col.names=F, quote=F)
	}

	#WITH LABELS
	#save merged gene list txt file with module labels
	module_names = setdiff(unique(dynamicLabels), "grey")
	gene_names = rownames(t(datExpr2))
	for(label in module_names){
	  module=gene_names[which(dynamicLabels==label)]
	  write.table(module, paste(label, ".txt", sep=""), sep="\t", row.names=F, col.names=F, quote=F)
	}

	#adjacency already calulated in stepwise module detection
	#Summary output of network analysis results for colors in moduleColors
	degrees1=intramodularConnectivity(adjacency, dynamicColors)
	head(degrees1)

	#want to make summary file with genes, genes intramodular connectivity, dynamicColors and dynamicLabels
	mydata = t(datExpr2)
	mydata2 = cbind(mydata, degrees1, dynamicColors, dynamicLabels)
	write.csv(mydata2,file=paste0(identifier,"_datExpr2_ksummary_dynamiccolors_dynamiclabels.csv"))

	#with colors
	# Calculate eigengenes
	MEList = moduleEigengenes(datExpr2, colors = dynamicColors, softPower= softPower, nPC=1)
	MEs = MEList$eigengenes

	pdf(file=paste0("7-",identifier,"_eigengenes-b_dynamic.pdf"),height=9,width=9)
	par(cex = 0.6)
	plotEigengeneNetworks(MEs, "", marDendro = c(0,8,1,2), marHeatmap = c(5,8,1,2), cex.adjacency = 0.3 , cex.preservation = 0.3, plotPreservation = "standard")
	#dev.print(pdf,"7-identifier_eigengenes-b_dynamic.pdf",height=9,width=9)
	dev.off()

	#with labels
	# Calculate eigengenes
	MEListL = moduleEigengenes(datExpr2, colors = dynamicLabels, softPower= softPower, nPC=1)
	MEsL = MEListL$eigengenes

	pdf(file=paste0("7L-",identifier,"_eigengenes-b_dynamic.pdf"),height=9,width=9)
	par(cex = 0.6)
	plotEigengeneNetworks(MEsL, "", marDendro = c(0,8,1,2), marHeatmap = c(5,8,1,2), cex.adjacency = 0.3 , cex.preservation = 0.3, plotPreservation = "standard")
	#dev.print(pdf,"7L-identifier_eigengenes-b_dynamic.pdf",height=9,width=9)
	dev.off()

	#Calculate signedkME. Higher the signed kME for a given gene higher its correlation with eigen gene i.e. its a hub gene for the given module Ref: Miller 2010, Gandal 2018
	geneMM=signedKME(datExpr2,MEs)
	colnames(geneMM)=paste(colnames(geneMM),".cor",sep="")
	MMPvalue=corPvalueStudent(as.matrix(geneMM),dim(datExpr2)[[1]])
	colnames(MMPvalue)=paste(colnames(MMPvalue),".pval",sep="")

	geneMML=signedKME(datExpr2,MEsL)
	colnames(geneMML)=paste(colnames(geneMML),".cor",sep="")
	MMPvalueL=corPvalueStudent(as.matrix(geneMML),dim(datExpr2)[[1]])
	colnames(MMPvalueL)=paste(colnames(MMPvalueL),".pval",sep="")

	#want to make summary file with genes, signedkME, blockColors and blockLabels
	mydata3 = cbind(geneMM, MMPvalue, dynamicColors)
	mydata4 = cbind(geneMML, MMPvalueL, dynamicLabels)
	write.csv(mydata3,file=paste0(identifier,"_datExpr2_signedKME_dynamiccolors.csv"))
	write.csv(mydata4,file=paste0(identifier,"_datExpr2_signedKME_dynamiclabels.csv"))

	write.csv(table(dynamicColors), file=paste0(identifier,"_final_dynamiccolors.csv"))
	write.csv(table(dynamicLabels), file=paste0(identifier,"_final_dynamiclabels.csv"))
	#these make an csv file with the module colors or my labels and how many members are in each

	#Step14 Visualization of final modules, PC1 and eigengenes
	#heatmaps for module colors
	aa = read.csv(paste0(identifier,"_final_dynamiccolors.csv"), header=T, row.names=1)
	MEss = MEs[,sort(colnames(MEs))]
	dim(aa)

	pdf(file=paste0("8-",identifier,"_heatmaps_dynamic.pdf"))
	for (k in 1:length(rownames(aa))) {
		whichmodule = aa$dynamicColors [k]
		par(mfrow=c(2,1))
		par(mar=c(3.5,3,2.5,3))
		par(oma=c(4,0,2,0))
		datcombined=datExpr2[,dynamicColors==whichmodule]
		datcombined=datcombined[order(row.names(datcombined)),]
		#sort so samples in order
		scaledModExpr =scale(datcombined)
		myEx = max(abs(scaledModExpr))
		plotMat(t(scaledModExpr), main=paste(aa$dynamicColors [k], "(", aa$Freq[k], ")", sep=""))
		names=sort(dimnames(datExpr2)[[1]])
		#sort labels
		numSamp = length(rownames(datExpr2))
		datsv=MEss[,k]
		whichmodule=aa$dynamicColors[k]
		space=0.174
		width=(10.65-(space*numSamp))/numSamp
		xpos=seq(0.35+(width/2),11,width+space)
		barplot(datsv,col="grey",xlim=c(0.35,11), width=width, space=space)
		axis(side=1,at=xpos,labels=F,tick=T)
		box()
		mtext(names,1,at=xpos,cex=1,col="black",adj=1.2,las=3)
		print(k)
	}
	dev.off()

	#heatmaps for module LABELS
	aaL = read.csv(paste0(identifier,"_final_dynamiclabels.csv"), header=T, row.names=1)
	MEssL = MEsL[,sort(colnames(MEsL))]
	dim(aaL)

	pdf(file=paste0("8L-",identifier,"_heatmaps_dynamic.pdf"))
	for (k in 1:length(rownames(aaL))) {
	  	whichmodule = aaL$dynamicLabels[k]
		par(mfrow=c(2,1))
		par(mar=c(3.5,3,2.5,3))
		par(oma=c(4,0,2,0))
		datcombined=datExpr2[,dynamicLabels==whichmodule]
		#datcombined=datcombined[order(row.names(datcombined)),] #sort so samples in order
		scaledModExpr =scale(datcombined)
		#myEx = max(abs(scaledModExpr))
		#plotMat(t(scaledModExpr), main=paste(aaL$dynamicModsLabels[k], "(", aaL$Freq[k], ")", sep=""))
		plotMat(t(scaledModExpr), main=paste(aaL$dynamicLabels[k], "(", aaL$Freq[k], ")", sep=""))
		names=rownames(datcombined)
		numSamp = length(rownames(datExpr2))
		datsv=MEssL[,k]
		whichmodule=aaL$dynamicLabels[k]
		space=0.05
		width=(11-(space*numSamp))/numSamp
		#edited on 20201215 by DT
		#xpos=seq(0.35+(width/2),11,width+space)
		#barplot(datsv,col="grey",xlim=c(0.35,11), width=width, space=space)
		#space in barplot takes that fraction of width variable, so we need to divide by width in order to make space equal space
		barplot(datsv,col="grey",xlim=c(0.41,10.6), space=space/width, width=width, axisnames=TRUE, axis.lty=1, names.arg=names, las=2)
		#axis(side=1,at=xpos,labels=F,tick=T)
		#we include the x axis labels in the barplot using names.arg for the names and axis.lty for the tick marks
		box()
		#mtext(names,1,at=xpos,cex=1,col="black",adj=1.2,las=3)

		print(k)
	}
	dev.off()

	#Note if goodsamplegenes removes samples then Count or data variable sample definition won't hold true so we need to redefine the samples based on datExpr2 here
	samples=rownames(datExpr2)
	samples

	#single plots module COLORS (moduleColors)
	PCvalues=MEs
	PCvalues= PCvalues[,order(names(PCvalues),decreasing=F)]
	modules=colnames(PCvalues)
	colorsMod=gsub("ME","",modules)
	library(stringr)
	mains=paste0(identifier,"_", colorsMod, sep="")
	namesSamp=samples

	pdf(file=paste0("9-",identifier,"_singlePC_dynamic.pdf"),width=10,height=10)
	#this makes 5x5 graphs per page, with margins so labels are cut off the page
	par(mfrow=c(5,5))
	par(oma=c(2,1,2,1))
	for (mod in 1:length(PCvalues)) {
	plot(smooth.spline(PCvalues[,mod]), xlab="", type="n", ylim=c(-0.5, 0.5), ylab="First PC", axes=F)
	abline(h=0)
	axis(1, at=1:length(namesSamp), labels=namesSamp, cex.axis=1, las=3)
	axis(2, at=seq(-0.5, 0.5, 0.5))
	#par here lets you plot a second graph on the same plot
	par(new=TRUE)
	barplot(PCvalues[,mod],col="white", ylim=c(-0.5, 0.5), axes=F,main=mains[mod])
	SecNumber=length(samples)+0.5 ##1.5 to Secondnumber or SecNumber just centers the line alittle between around the samples
	lines(smooth.spline(x=c(1.5:SecNumber),PCvalues[,mod], spar=0.3), col=colorsMod[mod],  lwd=4)
	}
	dev.off()

	#single plots module LABELS (moduledynamicLabels)
	PCvaluesL=MEsL
	PCvaluesL= PCvaluesL[,order(names(PCvaluesL),decreasing=F)]
	modulesL=colnames(PCvaluesL)
	colorsModL=gsub("ME","",modulesL)
	colorsModL2=gsub(paste0(identifier,"_"),"",colorsModL)
	colorsModcolors = labels2colors(colorsModL2)
	mainsL=paste0(identifier,"_", colorsModL2, sep="")
	namesSampL =samples

	pdf(file=paste0("9L-",identifier,"_singlePC_dynamic.pdf"),width=10,height=10)
	#this makes 5x5 graphs per page, with margins so labels are cut off the page
	par(mfrow=c(5,5))
	par(oma=c(2,1,2,1))
	for (mod in 1:length(PCvaluesL)) {
	plot(smooth.spline(PCvaluesL[,mod]), xlab="", type="n", ylim=c(-0.5, 0.5), ylab="First PC", axes=F)
	abline(h=0)
	axis(1, at=1:length(namesSampL), labels=namesSampL, cex.axis=1, las=3)
	axis(2, at=seq(-0.5, 0.5, 0.5))
	#par here lets you plot a second graph on the same plot
	par(new=TRUE)
	barplot(PCvaluesL[,mod],col="white", ylim=c(-0.5, 0.5), axes=F,main=mainsL[mod])
	#keeping the module colors for the line, but using numbers for the main labels
	SecNumber=length(samples)+0.5 ##1.5 to Secondnumber or SecNumber just centers the line alittle between around the samples
	lines(smooth.spline(x=c(1.5:SecNumber),PCvaluesL[,mod], spar=0.3), col=colorsModcolors[mod],  lwd=4)
	}
	dev.off()

	#Step15 plot module-traits relation
	#IMPORTANT CAUTION: In the Trait file the samples have to be the same order as the output in the column order in the count file, so make the Trait.csv after the Count file has been generated
	dim(datExpr2)
	head(datExpr2[,1:3])
	#datExpr2 contains the samples as rows and most variable probes 1 to 10000 as columns

	nGenes = ncol(datExpr2);
	nSamples = nrow(datExpr2);
	nGenes
	nSamples

	#Import trait data
	head(traitData)

	#Normally may have different traits so here select whatever traits you need
	#Note if goodsamplegenes removes samples then Count or data variable sample definition won't hold true so we need to redefine the samples based on datExpr2 here
	save.image(file='4.32_identifier_settings4.25t4.26e_MultipleTraitColumns.RData')
	Traits=traitData[match(rownames(datExpr2), traitData$Sample),-1]
	rownames(Traits)=traitData[match(rownames(datExpr2), traitData$Sample),c("Sample")]
	#datTraits=Traits[,-which(names(Traits) %in% c("sample"))]
	datTraits=as.data.frame(Traits)
	head(datTraits)

	#This is for color labeled modules
	moduleTraitCor = cor(MEs, datTraits, use = "p");
	moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
	colnames(moduleTraitPvalue) = paste("p.value.", colnames(moduleTraitCor), sep="");
	out3<-cbind(Module=rownames(moduleTraitCor ), moduleTraitCor, moduleTraitPvalue)
	dim(out3)
	write.table(out3, paste0(identifier,"_heatmap_moduleTraitCor.csv"), sep=",",row.names=F)

	pdf(file=paste0("10-",identifier,"_heatmap_moduleTraitCor_dynamic.pdf"), width=5,height=15)
	textMatrix = paste( signif(moduleTraitCor, 2), '\n(',
						   signif(moduleTraitPvalue, 1), ')', sep = ''
						   );
	dim(textMatrix) = dim(moduleTraitCor);
	par( mar = c(8, 9, 3, 3) );
	labeledHeatmap( Matrix = moduleTraitCor,
					xLabels = names(datTraits),
					yLabels = names(MEs),
					ySymbols = names(MEs),
					colorLabels = F,
					colors = blueWhiteRed(50),
					textMatrix = textMatrix,
					setStdMargins = F,
					cex.text = 0.25,
					zlim = c(-1, 1),
	        cex.lab.y = 0.25,
					main = paste("module-trait relationships")
					);
	sig = moduleTraitPvalue < (.05 / (ncol(moduleTraitPvalue) * nrow(moduleTraitPvalue)));
	#dev.print(pdf,"10-identifier_heatmap_moduleTraitCor_dynamic.pdf", width=5,height=10)
	dev.off()
	#In the exported pdf this plot looks better

	#This is for dynamic labeled named modules
	#Note datTraits file and nSamples remains the same
	moduleTraitCorL = cor(MEsL, datTraits, use = "p");
	moduleTraitPvalueL = corPvalueStudent(moduleTraitCorL, nSamples);
	colnames(moduleTraitPvalueL) = paste("p.value.", colnames(moduleTraitCorL), sep="");
	out3L<-cbind(Module=rownames(moduleTraitCorL ), moduleTraitCorL, moduleTraitPvalueL)
	dim(out3L)
	write.table(out3L, paste0(identifier,"_heatmap_modulelabelsTraitCor.csv"), sep=",",row.names=F)

	pdf(file=paste0("10L-",identifier,"_heatmap_moduleTraitCor_dynamic.pdf"), width=5,height=15)
	textMatrixL = paste( signif(moduleTraitCorL, 2), '\n(',
	                     signif(moduleTraitPvalueL, 1), ')', sep = ''
	);
	dim(textMatrixL) = dim(moduleTraitCorL);
	par( mar = c(8, 9, 3, 3) );
	labeledHeatmap( Matrix = moduleTraitCorL,
	                xLabels = names(datTraits),
	                yLabels = names(MEsL),
	                ySymbols = names(MEsL),
	                colorLabels = F,
	                colors = blueWhiteRed(50),
	                textMatrix = textMatrixL,
	                setStdMargins = F,
	                cex.text = 0.25,
	                zlim = c(-1, 1),
	                cex.lab.y = 0.25,
	                main = paste("module-trait relationships")
	);
	sig = moduleTraitPvalueL < (.05 / (ncol(moduleTraitPvalueL) * nrow(moduleTraitPvalueL)));
	#dev.print(pdf,"10L-identifier_heatmap_moduleTraitCor_dynamic.pdf", width=5,height=10)
	dev.off()

	###Single Principal Component Plots with sample Trait###
	#Here we generate plots same as 9 and 9L, but here the sample numbers are replaced by the trait of our interest
	Trait_chr=traitData[match(rownames(datExpr2), traitData$X),c("interest_trait_chr")]
	Trait_chr

	#single plots module COLORS (moduleColors)
	PCvalues=MEs
	PCvalues= PCvalues[,order(names(PCvalues),decreasing=F)]
	modules=colnames(PCvalues)
	colorsMod=gsub("ME","",modules)
	library(stringr)
	mains=paste0(identifier,"_", colorsMod, sep="")
	namesSamp=Trait_chr

	pdf(file=paste0("9.1-",identifier,"_singlePC_dynamic.pdf"),width=10,height=10)
	#this makes 5x5 graphs per page, with margins so labels are cut off the page
	par(mfrow=c(5,5))
	par(oma=c(2,1,2,1))
	for (mod in 1:length(PCvalues)) {
	plot(smooth.spline(PCvalues[,mod]), xlab="", type="n", ylim=c(-0.5, 0.5), ylab="First PC", axes=F)
	abline(h=0)
	axis(1, at=1:length(namesSamp), labels=namesSamp, cex.axis=0.5, las=3)
	axis(2, at=seq(-0.5, 0.5, 0.5))
	#par here lets you plot a second graph on the same plot
	par(new=TRUE)
	barplot(PCvalues[,mod],col="white", ylim=c(-0.5, 0.5), axes=F,main=mains[mod])
	SecNumber=length(samples)+0.5 ##1.5 to Secondnumber or SecNumber just centers the line alittle between around the samples
	lines(smooth.spline(x=c(1.5:SecNumber),PCvalues[,mod], spar=0.3), col=colorsMod[mod],  lwd=4)
	}
	dev.off()

	#single plots module LABELS (moduledynamicLabels)
	PCvaluesL=MEsL
	PCvaluesL= PCvaluesL[,order(names(PCvaluesL),decreasing=F)]
	modulesL=colnames(PCvaluesL)
	colorsModL=gsub("ME","",modulesL)
	colorsModL2=gsub(paste0(identifier,"_"),"",colorsModL)
	colorsModcolors = labels2colors(colorsModL2)
	mainsL=paste(identifier,"_", colorsModL2, sep="")
	namesSampL=Trait_chr

	pdf(file=paste0("9L.1-",identifier,"_singlePC_dynamic.pdf"),width=10,height=10)
	#this makes 5x5 graphs per page, with margins so labels are cut off the page
	par(mfrow=c(5,5))
	par(oma=c(2,1,2,1))
	for (mod in 1:length(PCvaluesL)) {
	plot(smooth.spline(PCvaluesL[,mod]), xlab="", type="n", ylim=c(-0.5, 0.5), ylab="First PC", axes=F)
	abline(h=0)
	axis(1, at=1:length(namesSampL), labels=namesSampL, cex.axis=0.5, las=3)
	axis(2, at=seq(-0.5, 0.5, 0.5))
	#par here lets you plot a second graph on the same plot
	par(new=TRUE)
	barplot(PCvaluesL[,mod],col="white", ylim=c(-0.5, 0.5), axes=F,main=mainsL[mod])
	#keeping the module colors for the line, but using numbers for the main labels
	SecNumber=length(samples)+0.5 ##1.5 to Secondnumber or SecNumber just centers the line alittle between around the samples
	lines(smooth.spline(x=c(1.5:SecNumber),PCvaluesL[,mod], spar=0.3), col=colorsModcolors[mod],  lwd=4)
	}
	dev.off()


	###check if modules have significant correlation and p-value to trait of interest###
	#Mod_corrTrait=read.csv(paste0(identifier,"_heatmap_modulelabelsTraitCor.csv"), header=TRUE, sep=',')
	#head(Mod_corrTrait)

	#Mod_corrTrait_interest_trait_code_sig=Mod_corrTrait[Mod_corrTrait$p.value.interest_trait_code<=0.05,]
	#Mod_corrTrait_interest_trait_code_sig=Mod_corrTrait_interest_trait_code_sig[order(Mod_corrTrait_interest_trait_code_sig$p.value.interest_trait_code),]
	#By default, sorting is ASCENDING
	#Mod_corrTrait_interest_trait_code_sig
	#write.csv(Mod_corrTrait_interest_trait_code_sig,"Mod_corrTrait_interest_trait_code_sig.csv")

	###check if number of modules are too many or to few and if biggest module has too many or too few genes
	#Import new module data we got from our analysis containing genes, VST counts, module and module information summary file from this analysis folder identifier_datExpr2_ksummary_dynamiccolors_dynamiclabels.csv
	dat1=read.csv(paste0(identifier,"_datExpr2_ksummary_dynamiccolors_dynamiclabels.csv"), sep=',')
	head(dat1)

	dat1GeneMod=dat1[,c("X","dynamicColors","dynamicLabels")]
	colnames(dat1GeneMod)=c("Gene","dynamicColors","dynamicLabels")
	dat1GeneMod[1:3,]
	dat1GeneMod_table<-as.data.frame(table(dat1GeneMod$dynamicLabels))
	dat1GeneMod_table<-dat1GeneMod_table[order(-dat1GeneMod_table$Freq),]
	write.csv(dat1GeneMod_table,paste0(identifier,"_final_dynamiccolors_sorted.csv"))
	#This result is same as "identifier_final_dynamiccolors.csv" but sorted

	dat1GeneMod_table2<-as.data.frame(table(dat1GeneMod$dynamicColors))
	dat1GeneMod_table2<-dat1GeneMod_table2[order(-dat1GeneMod_table2$Freq),]
	write.csv(dat1GeneMod_table2,paste0(identifier,"_final_dynamiclabels_sorted.csv"))
	#This result is same as "identifier_final_dynamiclabels.csv" but sorted

	#For dynamicColors
	#This is same as the dat1 file imported above
	dat1=dat1
	#Export genes ranked by net connectivity weight or kIM or IM for each module
	for (i in unique(dat1GeneMod$dynamicColors)){

	  #Select modules
	  modules= i
	  modules

	  #net gene Connectivity or weights is a.k.a. kIM or IM or kWithin
	  dat1_sel=dat1[which(dynamicColors==modules),] #select those rows which match module name
	  dat1_sort=dat1_sel[order(-dat1_sel$kWithin),] #order by stronger to weaker connectivity

	  #net gene Connectivity ranks
	  dat1_sort$rank<- rank(-(dat1_sort$kWithin)) #This gives the highest kWithin the highest rank
	  #export output
	  write.table(dat1_sort, paste0("ForGenesInModule-", modules, "-allGenesInModuleRankedByNetKWithin.txt", sep=""),row.names=F, quote=F)
	}

	#For dynamicLabels
	#This is same as the dat1 file imported above
	dat1=dat1
	#Export genes ranked by net connectivity weight or kIM or IM for each module
	for (i in unique(dat1GeneMod$dynamicLabels)){

	  #Select modules
	  modules= i
	  modules

	  #net gene Connectivity or weights is a.k.a. kIM or IM or kWithin
	  dat1_sel=dat1[which(dynamicLabels==modules),] #select those rows which match module name
	  dat1_sort=dat1_sel[order(-dat1_sel$kWithin),] #order by stronger to weaker connectivity

	  #net gene Connectivity ranks
	  dat1_sort$rank<- rank(-(dat1_sort$kWithin)) #This gives the highest kWithin the highest rank
	  #export output
	  write.table(dat1_sort, paste0("ForGenesInModule-", modules, "-allGenesInModuleRankedByNetKWithin.txt", sep=""),row.names=F, quote=F)
	}

	#Import module names and corrosponding gene names from WGCNA summary file 'identifier_datExpr2_ksummary_dynamiccolors_dynamiclabels.csv' from analysis folder
	dat1=read.csv(paste0(identifier,"_datExpr2_ksummary_dynamiccolors_dynamiclabels.csv"), sep=',')
	dat1GeneMod=dat1[,c("X","dynamicColors","dynamicLabels")]
	colnames(dat1GeneMod)=c("Gene","dynamicColors","dynamicLabels")
	dat1GeneMod[1:3,]
	#For MEDissThresh >0.05 for smaller merged modules, multiple colors can be assigned to the same moduledynamicLabels. So its better to use the final moduledynamicLabels for exporting for cytoscape and VisANT

	#blockwise module TOM is "dist" matrix object so we need to convert it to what we need that is a regular matrix
	#our TOM is already loaded in the environment
	#TOM was calculated from 'adjacency'
	#Checking format of TOM, it needs to be a regular matrix
	typeof(TOM)
	is.matrix(TOM)
	str(TOM)
	class(TOM)

	#Export using already generated TOM and power
	TOM = TOM
	for (i in unique(dat1GeneMod$dynamicLabels)){
	  #We use original TOM that was used to build the modules. This is a matrix of matrices, with rows=columns=number of genes input into WGCNA.
	  TOM = TOM

	  #Select modules
	  modules= i
	  modules

	  #We have many diffent moduledynalicLabels this function only keeps those moduledynamicLabels that match with our selected module name
	  inModule = is.finite(match(dynamicLabels,modules))
	  #Select only those genes that are in the selected module. In keeping with microarrya nomenclature modProbes is the default name of the genes in modules    for the WGCNA package.
	  probes = names(datExpr2)
	  modProbes = probes[inModule]
	  #modProbes = substring(modProbes,1,25)
	  length(modProbes)
	  #this number should match with the number of genes in the selected module name

	  #Select the corresponding Topological Overlap for the selected module i.e. only those genes are selected that belong to the selected module
	  modTOM = TOM[inModule, inModule];
	  dimnames(modTOM) = list(modProbes, modProbes)

	  #Export the network into files for VisANT for the selected module.
	  vis = exportNetworkToVisANT(modTOM,file = paste("VisANTInput-", modules, "-allConnected.txt", sep=""),weighted = TRUE, threshold = 0)
	  vis[1:3,]

	  #Using this files as input for visANT will be too many just look like a green blob
	  #Instead us this file as input for VisANT which is a selection of top400 and top 800 most connected genes.
	  #Connectivity or weights is a.k.a. kIM or IM or kWithin
	  vis_sort=vis[order(-vis$weight),] #order by stronger to weaker connectivity
	  #export output
	  write.table(vis_sort, paste0("VisANTInput-", modules, "-allConnectedSorted.txt", sep=""),row.names=F, col.names=F, quote=F)
	  #Connectivity ranks
	  vis_rank=vis_sort
	  colnames(vis_rank)=c("FromGene", "ToGene", "thresold", "VisANTID", "Weights_kIM_IM")
	  vis_rank$rank<- rank(-(vis_rank$Weights_kIM_IM)) #This gives the highest weight the highest rank
	  #export output
	  write.table(vis_rank, paste0("VisANTRanks-", modules, "-allConnectedRankedByWeight.txt", sep=""),row.names=F, quote=F)

	  vis_top400wts=vis_sort[1:400,] #select only the top 400
	  #export output
	  write.table(vis_top400wts, paste0("VisANTInput-", modules, "-top400ConnectedSorted.txt", sep=""),row.names=F, col.names=F, quote=F)
	  #Connectivity ranks
	  vis_rank_top400wts=vis_top400wts
	  colnames(vis_rank_top400wts)=c("FromGene", "ToGene", "thresold", "VisANTID", "Weights_kIM_IM")
	  vis_rank_top400wts$rank<- rank(-(vis_rank_top400wts$Weights_kIM_IM)) #This gives the highest weight the highest rank
	  #export output
	  write.table(vis_rank_top400wts, paste0("VisANTRanks-", modules, "-top400ConnectedRankedByWeight.txt", sep=""),row.names=F, quote=F)

	  vis_top800wts=vis_sort[1:800,] #select only the top 800
	  #export output
	  write.table(vis_top800wts, paste0("VisANTInput-", modules, "-top800ConnectedSorted.txt", sep=""),row.names=F, col.names=F, quote=F)
	  #Connectivity ranks
	  vis_rank_top800wts=vis_top800wts
	  colnames(vis_rank_top800wts)=c("FromGene", "ToGene", "thresold", "VisANTID", "Weights_kIM_IM")
	  vis_rank_top800wts$rank<- rank(-(vis_rank_top800wts$Weights_kIM_IM)) #This gives the highest weight the highest rank
	  #export output
	  write.table(vis_rank_top800wts, paste0("VisANTRanks-", modules, "-top800ConnectedRankedByWeight.txt", sep=""),row.names=F, quote=F)

	}

	#Cytoscape Export 1

	#Export using already generated TOM and power
	TOM = TOM
	for (i in unique(dat1GeneMod$dynamicLabels)){
	  #We use original TOM that was used to build the modules. This is a matrix of matrices, with rows=columns=number of genes input into WGCNA.
	  TOM = TOM

	  #Select modules
	  modules= i
	  modules

	  #We have many diffent moduledynalicLabels this function only keeps those moduledynamicLabels that match with our selected module name
	  inModule=is.finite(match(dynamicLabels,modules))
	  #Select only those genes that are in the selected module
	  genes=names(datExpr2)
	  modGenes=genes[inModule]
	  length(modGenes)
	  #this number should match with the number of genes in the selected module name

	  #Select the corresponding Topological Overlap for the selected module i.e. only those genes are selected that belong to the selected module
	  modTOM = TOM[inModule, inModule]
	  dimnames(modTOM) = list(modGenes, modGenes)

	  #Export the network into edge and node list files for Cytoscape for the selected module and use as input for cytoscape
	  cyt = exportNetworkToCytoscape(modTOM,
	                                 edgeFile=paste("CytoEdge",paste(modules,collapse="-"),"allConnected.txt",sep=""),
	                                 nodeFile=paste("CytoNode",paste(modules,collapse="-"),"allConnected.txt",sep=""),
	                                 weighted = TRUE, threshold = 0.02,nodeNames=modGenes,
	                                 altNodeNames = modGenes, nodeAttr = dynamicLabels[inModule])
	}

	#Step17 Organization of files
	library(filesstrings)
	dir.create("ForCytoscapeInput")
	file.move(list.files(pattern='CytoEdge*'), "ForCytoscapeInput")
	file.move(list.files(pattern='CytoNode*'), "ForCytoscapeInput")

	dir.create("ForVisANTInput")
	file.move(list.files(pattern='VisANTInput-.*'), "ForVisANTInput")
	file.move(list.files(pattern='VisANTRanks-.*'), "ForVisANTInput")

	dir.create("GeneLists")
	file.move(list.files(pattern='ForGenesInModule-.*'), "GeneLists")
	file.move(list.files(pattern = paste0(identifier,"_.*.txt")), "GeneLists")

	dir.create("Plots")
	file.move(list.files(pattern = '*.pdf'), "Plots")

	dir.create("Traits")
	file.move(list.files(pattern = "Mod_.*.sig.csv"), "Traits")
	file.move(list.files(pattern = "*TraitCor.csv"), "Traits")

	dir.create("Normalizations")
	file.move(list.files(pattern = "_log*"), "Normalizations")
	file.move(list.files(pattern = "_CPM.csv"), "Normalizations")
	file.move(list.files(pattern = "_RPKM.csv"), "Normalizations")
	file.move(list.files(pattern = "_TPM.csv"), "Normalizations")
	file.move(list.files(pattern = "_VST.csv"), "Normalizations")
	file.move(list.files(pattern = "_Count_merged_matchedToGtf.csv"), "Normalizations")

	#remove extra files
	file.remove(list.files(pattern = "final_modules.csv"))
	file.remove(list.files(pattern = "final_modules_labels.csv"))

	dir.create("ModuleSummary")
	file.move(list.files(pattern = "_datExpr2*"), "ModuleSummary")
	file.move(list.files(pattern = "_final_dynamiccolors.csv"), "ModuleSummary")
	file.move(list.files(pattern = "_final_dynamiclabels.csv"), "ModuleSummary")
	file.move(list.files(pattern = "_final_dynamiccolors_sorted.csv"), "ModuleSummary")
	file.move(list.files(pattern = "_final_dynamiclabels_sorted.csv"), "ModuleSummary")

	#remove extra files
	file.remove(list.files(pattern="dat1GeneLabel.csv"))
	file.remove(list.files(pattern="tdatExpr.csv"))
	file.remove(list.files(pattern="PaperTablewithExpr.csv"))
	file.remove(list.files(pattern="*.RData")) #optional deletion of .RData

	#Step18 Save version numbers for reproducibility and version control
	sessionInfo()
	toLatex(sessionInfo())
	if(saveRData) save.image(file=paste0(identifier,'.RData'))

	return(new("WGCNA", datExpr=mydata2, conditions=traitData, trait=out3L))
}

