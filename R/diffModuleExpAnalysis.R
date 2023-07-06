#' Plots an expression profile for a module
#'
#' A plotting function that returns a heatmap and barplot for a module
#'
#' @param WGCNAobject an object of class WGCNAobject
#' @param geneList a vector of gene names to be extracted from WGCNAobject
#' @param mode use first principal component or averageZscore?
#' @param legend plot legend?
#' @param title title of the plot
#' @param clusterGenes cluster heatmap genes by hierarchical clustering?
#'
#' @return a patchworked ggplot object
#' 
#' @author Dario Tommasini
#'
#' @import ggplot2
#' @import patchwork
#' @import WGCNA
#' @import dplyr
#' @export
#' 
#' @examples
#' library(ExperimentHub)
#' eh = ExperimentHub()
#' eh_query = query(eh, c("multiWGCNAdata"))
#' astrocyte_networks = eh_query[["EH8222"]]
#' moduleExpressionPlot(astrocyte_networks[["combined"]], 
#'   geneList = topNGenes(astrocyte_networks$combined, "combined_013"))
#' 
moduleExpressionPlot <- function(WGCNAobject, geneList, mode=c("PC1", "averageZscore"), legend=FALSE, title=NULL, clusterGenes=FALSE){

  # Check that mode is either PC1 or averageZscore
  mode = match.arg(mode)
  
	datExpr=WGCNAobject@datExpr
	design=WGCNAobject@conditions
	heatmap <- expressionHeatmap(datExpr, 
	                             geneList, 
	                             column.labels=FALSE, 
	                             axis.titles=FALSE, 
	                             legend=legend, 
	                             clusterGenes=clusterGenes, 
	                             title=title)
	cleanDatExpr=t(cleanDatExpr(datExpr))
	subset=cleanDatExpr[toupper(rownames(cleanDatExpr)) %in% toupper(geneList),]

	if(mode=="averageZscore"){
		mean=rowMeans(subset)
		stdev=apply(subset, 1, sd)
		zscoreMatrix=(subset-mean)/stdev
		zscoreMatrix=na.omit(zscoreMatrix)
		moduleExpression=reshape2::melt(zscoreMatrix)
		colnames(moduleExpression)=c("Gene", "Sample", "Expression")
		mergedData=cbind(moduleExpression, design[match(moduleExpression$Sample, design$X),-1])
		boxplot <- ggplot(data=mergedData, aes(x=Sample, y=Expression, color=eval(parse(text=colnames(design)[[2]])))) +
    				geom_boxplot() +
    				guides(color=guide_legend(title="Condition")) +
    				{if(!legend) theme(legend.position = "none")} +
    				theme(panel.background = element_blank(), axis.line = element_line(colour = "black"),
    					panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    					axis.text.x=element_text(angle=45, hjust=1), plot.title=element_text(hjust=0.5))

	return(heatmap / boxplot)
	}

	if(mode=="PC1"){
		PC1=moduleEigengenes(t(subset), colors = rep("Module", nrow(subset)), nPC=1)$eigengenes
		moduleExpression=data.frame(Sample=rownames(PC1), moduleExpression=PC1)
		colnames(moduleExpression)=c("Sample","moduleExpression")
		mergedData=cbind(moduleExpression, design[match(moduleExpression$Sample, design$X),-1])
		bargraph <- ggplot(data=moduleExpression, aes(x=factor(Sample, levels=Sample), y=moduleExpression))+
					ylab(mode) +
					geom_bar(aes(fill= moduleExpression), stat="identity", color="black", position=position_dodge(9)) +
	       			scale_fill_gradient2(name=mode, low = "blue", mid="white", high = "red", midpoint=0) +
	     		  	theme(axis.ticks.x=element_blank(), axis.text.x=element_text(angle=90, hjust=1),
	     		  	panel.background = element_blank(), axis.line.y = element_line(colour = "black"),
	       			plot.title=element_text(hjust=0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
	       			axis.title.x=element_blank()) +
	       			{if(!legend) theme(legend.position = "none")} +
	      			geom_hline(yintercept=0, linetype="solid", color = "black", size=0.5)
	return(heatmap / bargraph)
	}
}

# heatmap of z-scored expression for a given expression data.frame and list of genes
expressionHeatmap <- function(datExpr, 
                              geneList, 
                              lower=-2, 
                              upper=2, 
                              design=NULL, 
                              column.labels=TRUE, 
                              row.labels=FALSE, 
                              legend=FALSE, 
                              axis.titles=TRUE, 
                              plotTitle=NULL, 
                              clusterGenes=FALSE, 
                              title=NULL){
	datExpr=t(cleanDatExpr(datExpr))
	subset=datExpr[toupper(rownames(datExpr)) %in% toupper(geneList),]
	subset=subset[match(toupper(geneList), toupper(rownames(subset))),]
	subset=na.omit(subset)
	mean=rowMeans(subset)
	stdev=apply(subset, 1, sd)
	zscoreMatrix=(subset-mean)/stdev
	zscoreMatrix=na.omit(zscoreMatrix)
	zscoreMatrix[zscoreMatrix<lower]=lower
	zscoreMatrix[zscoreMatrix>upper]=upper
	collapsedMatrix=reshape2::melt(zscoreMatrix)
	colnames(collapsedMatrix)=c("Gene", "Sample", "Zscore")
	if(clusterGenes) {
		tree=hclust(dist(zscoreMatrix), method = "average")
		collapsedMatrix$Gene=factor(collapsedMatrix$Gene, levels=rev(rownames(zscoreMatrix)[tree$order]))
	}
	plt = ggplot(collapsedMatrix, aes(x = Sample, y = Gene, fill = Zscore)) +
				geom_tile() +
				ylab("Gene") +
				{if(!is.null(plotTitle)) ggtitle(plotTitle)}+
				theme_classic()+
				scale_fill_gradient2(low = "blue", mid="white", high = "red", na.value="white", limits=c(lower, upper)) +
				{if(!is.null(title)) ggtitle(title)} +
				{if(!axis.titles) theme(axis.title=element_blank())} +
				{if(!column.labels) theme(axis.text.x = element_blank()) else theme(axis.text.x = element_text(angle = 90, vjust=(0.5)))} +
				{if(!row.labels) theme(axis.text.y = element_blank())} +
				{if(!legend) theme(legend.position = "none")} +
				theme(axis.ticks=element_blank(), axis.line=element_blank(), plot.title=element_text(hjust=0.5))
	
	return(plt)
}

#' Run differential module expression
#'
#' A wrapper to run diffModuleExpression on all the modules in a network
#'
#' @param WGCNAobject object of class WGCNA with the modules to run DME on
#' @param design the sampleTable
#' @param alphaLevel level of significance
#' @param testCondition the column of the sampleTable to be resolved
#' @param refCondition the column of the sampleTable to be used as biological variation
#' @param test statistical test to perform, either "ANOVA" or "PERMANOVA"
#' @param p.adjust adjust for multiple comparisons, argument to pass to p.adjust function
#' @param plot generate a plot?
#' @param write write results to a file?
#' @param out file name for DME plots, only used if write is TRUE
#'
#' @return a data.frame summarizing the results of the analysis
#'
#' @author Dario Tommasini
#'
#' @importFrom data.table fwrite
#'
#' @export
#' 
#' @examples
#' library(ExperimentHub)
#' eh = ExperimentHub()
#' eh_query = query(eh, c("multiWGCNAdata"))
#' astrocyte_se = eh_query[["EH8223"]]
#' sampleTable = colData(astrocyte_se)
#' astrocyte_networks = eh_query[["EH8222"]]
#' runDME(astrocyte_networks[["combined"]], 
#'   design = sampleTable,
#'   p.adjust = "fdr", 
#'   refCondition = "Region", 
#'   testCondition = "Disease") 
#' 
runDME <- function(WGCNAobject, design, alphaLevel=0.05, testCondition=NULL, refCondition=NULL, p.adjust="fdr", plot=FALSE, test=c("ANOVA", "PERMANOVA"), write=FALSE, out=NULL){
	
  # Check argument
  test = match.arg(test)
  
  datExpr=WGCNAobject@datExpr
	if(is.null(refCondition)) refCondition=colnames(design)[[3]]
	if(is.null(testCondition)) testCondition=colnames(design)[[2]]
	modules=sort(unique(datExpr$dynamicLabels))
	modulePrefix=name(WGCNAobject)
	pval.dfs=list()
	if(plot) {
	  if(is.null(out)) pdf(paste0(modulePrefix,"_DME.pdf")) else pdf(paste0(out))
	}
	element=1
	for(module in modules){
		moduleGenes=datExpr$X[datExpr$dynamicLabels==module]
		pval.dfs[[element]]=diffModuleExpression(WGCNAobject, moduleGenes, design, plotTitle=module, plot=plot, test=test)
		colnames(pval.dfs[[element]])=c("Factors", module)
		element=element+1
	}
	if(plot) dev.off()
	mergedPval=Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "Factors", all.x = TRUE), pval.dfs)
	mergedPadj=as.matrix(rbind(p.adjust(mergedPval[1,-1], method=p.adjust),
			  	p.adjust(mergedPval[2,-1], method=p.adjust),
			  	p.adjust(mergedPval[3,-1], method=p.adjust)))
	rownames(mergedPadj)=mergedPval$Factors
	colnames(mergedPadj)=colnames(mergedPval[,-1])
	results.df=as.data.frame(t(mergedPadj))
	if(write) data.table::fwrite(results.df, paste0(modulePrefix, "/", modulePrefix, "_diffModExp.txt"), sep="\t", row.names=TRUE)

	return(results.df)
}

#' Differential module expression
#'
#' Runs (and plots) the differential module expression analysis
#'
#' @param WGCNAobject WGCNA object
#' @param geneList vector of genes in WGCNAobject
#' @param design the sampleTable
#' @param plotTitle title for the plot
#' @param mode either PC1 or Zscore, default is PC1
#' @param testColumn the column of the sampleTable to be resolved
#' @param refColumn the column of the sampleTable to be used as biological variation
#' @param test statistical test to perform, either "ANOVA" or "PERMANOVA"
#' @param plot generate a plot?
#' 
#' @return a data.frame with the resulting p-values
#'
#' @import patchwork
#' @import ggplot2
#' @import WGCNA
#' @import dplyr
#' @export
#' 
#' @examples
#' library(ExperimentHub)
#' eh = ExperimentHub()
#' eh_query = query(eh, c("multiWGCNAdata"))
#' astrocyte_se = eh_query[["EH8223"]]
#' sampleTable = colData(astrocyte_se)
#' astrocyte_networks = eh_query[["EH8222"]]
#' diffModuleExpression(astrocyte_networks[["combined"]], 
#'   topNGenes(astrocyte_networks$combined, "combined_013"), 
#'   sampleTable,
#'   test = "ANOVA",
#'   plotTitle = "combined_013",
#'   plot = TRUE)
#' 
diffModuleExpression <- function(WGCNAobject, 
                                 geneList, 
                                 design, 
                                 plotTitle=NULL, 
                                 mode=c("PC1", "Zscore"), 
                                 testColumn=2, 
                                 refColumn=3, 
                                 test=c("ANOVA", "PERMANOVA"), 
                                 plot=TRUE){
  
  # Check arguments
  mode = match.arg(mode)
  test = match.arg(test)

  datExpr=WGCNAobject@datExpr
  
	refCondition=colnames(design)[[refColumn]]
	testCondition=colnames(design)[[testColumn]]

	cleanDatExpr=t(cleanDatExpr(datExpr))
	geneList=geneList[geneList %in% rownames(cleanDatExpr)]
	subset=cleanDatExpr[rownames(cleanDatExpr) %in% geneList,]
	if(mode=="Zscore"){
		mean=rowMeans(subset)
		stdev=apply(subset,1,sd)
		zscoreMatrix=(subset-mean)/stdev
		zscoreMatrix=na.omit(zscoreMatrix)
		averageExpression=apply(zscoreMatrix, 2, mean)
		moduleExpression=data.frame(Sample=names(averageExpression), 
		                            moduleExpression=averageExpression)
	}

	if(mode=="PC1"){
		PC1=moduleEigengenes(t(subset), colors = rep("Module", length(geneList)), nPC=1)$eigengenes
		moduleExpression=data.frame(Sample=rownames(PC1), moduleExpression=PC1)
		colnames(moduleExpression)=c("Sample", "moduleExpression")
	}

	mergedData=cbind(moduleExpression, design[match(moduleExpression$Sample, design$Sample),-1])
	mergedData=mergedData %>% arrange(eval(parse(text=testCondition)), eval(parse(text=refCondition)))
	bargraph <- ggplot(data=mergedData, aes(x=factor(Sample, levels=Sample), y=moduleExpression))+
				ylab(mode) +
				geom_bar(aes(fill= moduleExpression), stat="identity", color="black", position=position_dodge(9)) +
	       			scale_fill_gradient2(name=mode, low = "blue", mid="white", high = "red", midpoint=0) +
	       			scale_x_discrete(labels=tolower(paste0(substr(mergedData[,3],0,3),"_",substr(mergedData[,4],0,3)))) +
	     		  	theme(axis.ticks.x=element_blank(), axis.text.x=element_text(angle=90, hjust=1),
	     		  	panel.background = element_blank(), axis.line.y = element_line(colour = "black"),
	       			plot.title=element_text(hjust=0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
	       			axis.title.x=element_blank()) +
	      			geom_hline(yintercept=0, linetype="solid", color = "black", size=0.5)

	datExpr=datExpr[,c("X", mergedData$Sample)]
	heatmap <- expressionHeatmap(datExpr, geneList, plotTitle=plotTitle, column.labels=FALSE, axis.titles=FALSE, legend=TRUE)

	if(test=="ANOVA"){
		pval.df <- performANOVA(moduleExpression, design, testCondition, refCondition) #category1=test, category2=ref)
	}
	
	if(test=="PERMANOVA"){
	  requireNamespace("vegan", quietly = TRUE)
	  permanova=vegan::adonis(t(subset) ~ design[[testCondition]] + design[[refCondition]] + design[[testCondition]]*design[[refCondition]], method = "euclidean", permutations = 9999)
	  Factors=c(testCondition, refCondition, paste0(testCondition, "*", refCondition))
	  p.value=c(permanova$aov.tab$`Pr(>F)`[1:3])
	  pval.df=data.frame(Factors, p.value)
	}

	boxplot <- ggplot(data=mergedData, aes(x=eval(parse(text=refCondition)), y=moduleExpression, color=eval(parse(text=testCondition)))) +
    				labs(title=paste0("p: ", paste(pval.df$Factors, signif(pval.df$p.value,2), sep = '=', collapse = ', ')),
   					y="Expression", x=refCondition) +
    				geom_boxplot(width=1/length(unique(mergedData[,refCondition]))) +
    				scale_color_manual(values=c("red", "blue")) +
    				guides(color=guide_legend(title=testCondition)) +
    				theme(panel.background = element_blank(), axis.line = element_line(colour = "black"),
    					panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title=element_text(hjust=0.5))

	if(plot){
		print(heatmap / bargraph / boxplot)
		message(paste0("#### plotting ", plotTitle, " ####\n"))
	}

	return(pval.df)
}

#' Perform ANOVA
#'
#' Test association between module expression to traits using ANOVA
#'
#' @param datExpr expression data.frame
#' @param testCondition test column in sampleTable
#' @param refCondition reference column in sampleTable
#' @param design the sampleTable
#' @param alphaLevel the significance level
#' 
#' @return a data.frame with p-values for each association
#'
#' @export
performANOVA <- function(datExpr, design, testCondition, refCondition, alphaLevel=0.05){ #category1, category2
	mergedData = cbind(datExpr, 
	                   test=design[match(datExpr$Sample, design$Sample), testCondition],
	                   ref=design[match(datExpr$Sample, design$Sample), refCondition])

	# Bug found 20230122: should use matched categorical variable from mergedData
	# H1=lm(as.numeric(datExpr$moduleExpression) ~ sampleTable[[testCondition]] + sampleTable[[refCondition]] +
	#       sampleTable[[testCondition]]*sampleTable[[refCondition]])
	# H0=lm(as.numeric(datExpr$moduleExpression) ~ sampleTable[[testCondition]] + sampleTable[[refCondition]])
	# H0_2=lm(as.numeric(datExpr$moduleExpression) ~ sampleTable[[testCondition]])
	
	H1=lm(as.numeric(mergedData$moduleExpression) ~ mergedData$test + mergedData$ref +
	        mergedData$test*mergedData$ref)
	H0=lm(as.numeric(mergedData$moduleExpression) ~ mergedData$test + mergedData$ref)
	H0_2=lm(as.numeric(mergedData$moduleExpression) ~ mergedData$test)
	
	Factors=c(testCondition, refCondition, paste0(testCondition, "*", refCondition))
	p.value=c(summary(H0)$coefficients[2,4], anova(H0, H0_2)[2,6], anova(H0, H1)[2,6])
	pval.df=data.frame(Factors, p.value)
  
	return(pval.df)
}

