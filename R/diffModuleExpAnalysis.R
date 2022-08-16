#' Plots an expression profile for a module
#'
#' A plotting function that returns a heatmap and barplot for a module
#'
#' @param WGCNAobject an object of class WGCNAobject
#' @param geneList a vector of gene names to be extracted from WGCNAobject
#' @param mode="PC1" use first principal component or averageZscore?
#' @param legend=FALSE plot legend?
#' @param title=NULL title of the plot
#' @param clusterGenes=F cluster heatmap genes by hierarchical clustering?
#'
#' @author Dario Tommasini
#'
#' @examples
#'
#' moduleExpressionPlot(WGCNAobject, moduleGenes, title=myModule)
#'
#' @import ggplot2
#' @import patchwork
#' @import WGCNA
#' @import dplyr
#' @export
moduleExpressionPlot <- function(WGCNAobject, geneList, mode="PC1", legend=FALSE, title=NULL, clusterGenes=F){

	datExpr=WGCNAobject@datExpr
	design=WGCNAobject@conditions
	heatmap <- expressionHeatmap(datExpr, geneList, column.labels=F, axis.titles=F, legend=legend, clusterGenes=clusterGenes, title=title)
	cleanDatExpr=t(cleanDatExpr(datExpr))
	subset=cleanDatExpr[toupper(rownames(cleanDatExpr)) %in% toupper(geneList),]

	if(mode=="averageZscore"){
		mean=rowMeans(subset)
		stdev=apply(subset, 1, sd)
		zscoreMatrix=(subset-mean)/stdev
		zscoreMatrix=na.omit(zscoreMatrix)
		#averageExpression=apply(zscoreMatrix, 2, mean)
		#moduleExpression=data.frame(Sample=names(averageExpression), moduleExpression=averageExpression)
		moduleExpression=reshape2::melt(zscoreMatrix)
		colnames(moduleExpression)=c("Gene", "Sample", "Expression")
		mergedData=cbind(moduleExpression, design[match(moduleExpression$Sample, design$X),-1])
		boxplot <- ggplot(data=mergedData, aes(x=Sample, y=Expression, color=eval(parse(text=colnames(design)[[2]])))) +
    				geom_boxplot() +
    				#scale_color_manual(values=c("red", "blue")) +
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
					#ggtitle(title) +
					geom_bar(aes(fill= moduleExpression), stat="identity", color="black", position=position_dodge(9)) +
	       			scale_fill_gradient2(name=mode, low = "blue", mid="white", high = "red", midpoint=0) +
	       			#scale_x_discrete(labels=tolower(paste0(substr(mergedData[,3],0,3), "_", substr(mergedData[,4],0,3)))) +
	     		  	theme(axis.ticks.x=element_blank(), axis.text.x=element_text(angle=90, hjust=1),
	     		  	panel.background = element_blank(), axis.line.y = element_line(colour = "black"),
	       			plot.title=element_text(hjust=0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
	       			axis.title.x=element_blank()) +
	       			{if(!legend) theme(legend.position = "none")} +
	      			geom_hline(yintercept=0, linetype="solid", color = "black", size=0.5)
	return(heatmap / bargraph)
	}

}


targetCumulativeDistribution <- function(DEobject, geneList, alternative="two.sided"){

	DEobject=na.omit(DEobject)
	DEobject$inList=FALSE
	DEobject$inList[toupper(rownames(DEobject)) %in% toupper(geneList)]=TRUE
	#cumulativeDat=data.frame(Gene=bothDat$X, kWithin=bothDat$kWithin, inTreatment=bothDat$X %in% treatDat$X)
	ksTestGreater=ks.test(DEobject$log2FoldChange[DEobject$inList], DEobject$log2FoldChange[!DEobject$inList], alternative="greater")
	ksTestLess=ks.test(DEobject$log2FoldChange[DEobject$inList], DEobject$log2FoldChange[!DEobject$inList], alternative="less")
	ksTest=ks.test(DEobject$log2FoldChange[DEobject$inList], DEobject$log2FoldChange[!DEobject$inList], alternative=alternative)
	print(ksTest)

	cumulative <- ggplot(DEobject, aes(log2FoldChange, colour=inList))+
                stat_ecdf(geom = "step") +
                ylab("Cumulative Fraction") +
                xlab("Log2 fold change") +
                xlim(-1,1)+
                scale_color_manual(values=c("black", "red")) +
                theme(panel.grid.major = element_blank(),legend.position="none",
                	panel.grid.minor = element_blank(),
                	panel.background=element_blank(),panel.border = element_rect(colour = "black", fill=NA)) +
                #theme(plot.title = element_text(hjust = 0.5)) +
                annotate("text", x = 0, y = 0.9, hjust=1, parse=T,
                	label = paste0("italic(P) == ", pretty10exp(signif(ksTestGreater$p.value,1)))) +
                annotate("text", x = 0, y = 0.2, hjust=0, parse=T,
                	label = paste0("italic(P) == ", pretty10exp(signif(ksTestLess$p.value,1))))

	cumulative
}



geneExpBarPlot <- function(datExpr, gene, clean=F, group=NULL){
	if(clean) {
		datExpr=cleanDatExpr(datExpr)
	} else {
		datExpr=t(datExpr)
	}
	expression=datExpr[,gene]
	df=data.frame(Sample=names(expression), Expression=expression, Condition=str_split_fixed(names(expression), "_", 2)[, 1])
	if(!is.null(group)){
		df$Group=as.character(group)
		df$Sample=factor(df$Sample, levels=df$Sample[order(df$Group)])
	}
	ggplot(df, aes(x= Sample, y=Expression, fill=Group, group=Group)) +
		geom_bar(stat="identity", color="black", position=position_dodge()) +
		#scale_fill_manual(values=colors(length(unique(df$Condition))))+
		xlab("Sample") +
		ylab("Expression") +
		theme_classic()+
		theme(legend.position="none", axis.text.x=element_text(angle=45, hjust=1))

}

#PCA of expression dataframe
PCAplot <- function(datExpr, WGCNA=FALSE, groups="none", labels="none"){
	if(WGCNA){
		data=cleanDatExpr(datExpr)
		design=sampleTable
		groups=as.factor(design[,3][match(rownames(data), design$Sample)])
	} else {
		data=t(datExpr)

	}

	#groups <- as.factor(c(rep('Cerebellum',6),rep('Frontal cortex',7),rep('Spinal cord',5),rep('Corpus callosum',6)))

	#remove genes with no variance
	data=data[, which(apply(data, 2, var) != 0)]
	res.pca <- prcomp(data, scale = T)
	fviz_pca_ind(res.pca,
             axes=c(1,2),
             axes.linetype=NA,
             #col.ind = groups, # color by groups
             #palette = c("cyan3","chartreuse2", "red3", "blue"),
             habillage=groups,
             invisible="quali",
             legend.title = ""#,
             #label=labels
  			 )+
             theme(axis.title=element_text(size=16), axis.text=element_text(size=12), panel.background = element_blank(),
             	axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
             	panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
             labs(x = "PC1", y = "PC2", title=NULL)
}

zScoreMatrix <- function(matrix){
	mean=rowMeans(matrix)
	stdev=apply(matrix, 1, sd)
	zScoreMatrix=(matrix-mean)/stdev
	na.omit(zScoreMatrix)
}

#heatmap of z-scored expression for a given expression dataframe and list of genes
expressionHeatmap <- function(datExpr, geneList, lower=-2, upper=2, design=NULL, column.labels=TRUE, row.labels=FALSE, legend=FALSE, axis.titles=TRUE, moduleName=NULL, clusterGenes=F, title=NULL){
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
	ggplot(collapsedMatrix, aes(x = Sample, y = Gene, fill = Zscore)) +
				geom_tile() +
				ylab("Gene") +
				{if(!is.null(moduleName)) ggtitle(moduleName)}+
				theme_classic()+
				scale_fill_gradient2(low = "blue", mid="white", high = "red", na.value="white", limits=c(lower, upper)) +
				{if(!is.null(title)) ggtitle(title)} +
				{if(!axis.titles) theme(axis.title=element_blank())} +
				{if(!column.labels) theme(axis.text.x = element_blank()) else theme(axis.text.x = element_text(angle = 90, vjust=(0.5)))} +
				{if(!row.labels) theme(axis.text.y = element_blank())} +
				{if(!legend) theme(legend.position = "none")} +
				theme(axis.ticks=element_blank(), axis.line=element_blank(), plot.title=element_text(hjust=0.5))

}

#' Best matching modules
#'
#' Find all the modules from dataset1 that have a best match to a module in dataset2
#' if that module in dataset2 is also a best match to the module in dataset1
#'
#'
#' @param overlapDf a data.frame resulting from a call to computeOverlapsFromWGCNA
#' @param plot generate a heatmap of best matching modules?
#'
#' @author Dario Tommasini
#'
#' @examples
#'
#' bidirectionalBestMatches(comparisonList[[element]]$overlap, WGCNAlist[[first]], WGCNAlist[[second]])
#'
#' @export
runDME <- function(WGCNAobject, alpha=get("alphaLevel", envir = parent.frame()), design=sampleTable, testCondition=NULL, refCondition=NULL, p.adjust="fdr", plot=FALSE, write=FALSE){
	datExpr=WGCNAobject@datExpr
	if(is.null(refCondition)) refCondition=colnames(design)[[3]]
	if(is.null(testCondition)) testCondition=colnames(design)[[2]]
	modules=sort(unique(datExpr$dynamicLabels))
	modulePrefix=name(WGCNAobject)
	pval.dfs=list()
	if(plot) pdf(paste0(modulePrefix,"_DME.pdf"))
	element=1
	for(module in modules){
		moduleGenes=datExpr$X[datExpr$dynamicLabels==module]
		pval.dfs[[element]]=diffModuleExpression(datExpr, moduleGenes, moduleName=module, plot=plot)
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
	if(write) fwrite(results.df, paste0(modulePrefix, "/", modulePrefix, "_diffModExp.txt"), sep="\t", row.names=T)

	return(results.df)
}

#' Differential module expression
#'
#' Runs (and plots if turned on) the differential module expression analysis
#'
#' @param datExpr expression data.frame
#' @param geneList vector of genes in datExpr
#' @param moduleName name of the module for title
#' @param mode either PC1 or Zscore, default is PC1
#' @param design the sampleTable
#' @param testColumn the column of the sampleTable to be resolved
#' @param refColumn the column of the sampleTable to be used as biological variation
#' @param test currently only ANOVA is supported
#' @param plot generate a plot?
#'
#' @import patchwork
#' @import ggplot2
#' @export
diffModuleExpression <- function(datExpr, geneList, moduleName=NULL, mode="PC1", design=sampleTable, testColumn=2, refColumn=3, test="ANOVA", plot=TRUE){
	#test=design[1, testColumn], ref=design[1, refColumn]

	refCondition=colnames(design)[[refColumn]]
	testCondition=colnames(design)[[testColumn]]
	#ref=refCondition
	#test=testCondition

	cleanDatExpr=t(cleanDatExpr(datExpr))
	geneList=geneList[geneList %in% rownames(cleanDatExpr)]
	subset=cleanDatExpr[rownames(cleanDatExpr) %in% geneList,]
	if(mode=="Zscore"){
		mean=rowMeans(subset)
		stdev=apply(subset,1,sd)
		zscoreMatrix=(subset-mean)/stdev
		zscoreMatrix=na.omit(zscoreMatrix)
		averageExpression=apply(zscoreMatrix, 2, mean)
		moduleExpression=data.frame(Sample=names(averageExpression), moduleExpression=averageExpression)
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
	heatmap <- expressionHeatmap(datExpr, geneList, moduleName=moduleName, column.labels=F, axis.titles=F, legend=TRUE)

	if(test=="ANOVA"){
		pval.df <- performANOVA(moduleExpression, testCondition, refCondition) #category1=test, category2=ref)
	}

	boxplot <- ggplot(data=mergedData, aes(x=eval(parse(text=refCondition)), y=moduleExpression, color=eval(parse(text=testCondition)))) +
    				labs(title=paste0("p: ", paste(pval.df$Factors, signif(pval.df$p.value,2), sep = '=', collapse = ', ')),
   					y="Expression", x=refCondition) +
    				#geom_hline(yintercept=0, linetype="solid", color = "black", size=0.5) +
    				geom_boxplot(width=1/length(unique(mergedData[,refCondition]))) +
    				scale_color_manual(values=c("red", "blue")) +
    				guides(color=guide_legend(title=testCondition)) +
    				theme(panel.background = element_blank(), axis.line = element_line(colour = "black"),
    					panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title=element_text(hjust=0.5))

	if(plot){
		print(heatmap / bargraph / boxplot)
		print(paste0("#### plotting ", moduleName, " ####"))
	}

	return(pval.df)
}

performANOVA <- function(expressionDf, testCondition, refCondition, design=get("design", envir = parent.frame()), alphaLevel=get("alphaLevel", envir = parent.frame())){ #category1, category2
	mergedData=cbind(expressionDf, test=design[match(expressionDf$Sample, design$Sample), testCondition],
		ref=design[match(expressionDf$Sample, design$Sample), refCondition])


	H1=lm(as.numeric(expressionDf$moduleExpression) ~ sampleTable[[testCondition]] + sampleTable[[refCondition]] +
	      sampleTable[[testCondition]]*sampleTable[[refCondition]])
	H0=lm(as.numeric(expressionDf$moduleExpression) ~ sampleTable[[testCondition]] + sampleTable[[refCondition]])
	H0_2=lm(as.numeric(expressionDf$moduleExpression) ~ sampleTable[[testCondition]])

	Factors=c(testCondition, refCondition, paste0(testCondition, "*", refCondition))
	p.value=c(summary(H0)$coefficients[2,4], anova(H0,H0_2)[2,6], anova(H0,H1)[2,6])
	pval.df=data.frame(Factors, p.value)

	pval.df
}


comparativeDME <- function(comparison, dataset1, dataset2, design=sampleTable, refCondition=3, testCondition=2, alphaLevel=get("alphaLevel", envir = parent.frame()), mode="PC1", test="zTestOfCorrelations", padjThreshold=10, overlapThreshold=30, useBestMatches=FALSE) {

	name1=str_split_fixed(comparison$overlap$mod1,"_",2)[,1][[1]]
	name2=str_split_fixed(comparison$overlap$mod2,"_",2)[,1][[1]]
	conditions=unique(design[,refCondition])

	top_overlaps=comparison$overlap[-log10(comparison$overlap$p.adj)>padjThreshold & comparison$overlap$overlap>overlapThreshold,]

	cat("#### Top Overlaps ####\n")
	print(head(top_overlaps))

	#calculate module eigengenes for each module
	ME1 = moduleEigengenes(cleanDatExpr(dataset1@datExpr), colors = dataset1@datExpr$dynamicLabels, nPC=1)
	ME2 = moduleEigengenes(cleanDatExpr(dataset2@datExpr), colors = dataset2@datExpr$dynamicLabels, nPC=1)

	if(mode=="averageExpr"){
		#moduleEigengenes1 = t(calculateAverageExpression(dataset1@datExpr))
		#moduleEigengenes2 = t(calculateAverageExpression(dataset2@datExpr))
		moduleEigengenes1=t(ME1$averageExpr)
		moduleEigengenes2=t(ME2$averageExpr)
		prefix="AE"
	}
	if(mode=="PC1"){
		moduleEigengenes1=t(ME1$eigengenes)
		moduleEigengenes2=t(ME2$eigengenes)
		prefix="ME"
	}

	collapseMeans <- function(moduleEigengeneTable, conditions){
		collapsedMeans=list()
		for(element in 1:length(conditions)) {
			collapsedMeans[[element]]=rowMeans(moduleEigengeneTable[, colnames(moduleEigengeneTable) %in%
				design$Sample[ design[,refCondition]==conditions[[element]] ] ])
		}
		expr1 = do.call(cbind,collapsedMeans)
		colnames(expr1)=conditions
		expr1
	}

	collapseStdDev <- function(moduleEigengeneTable,conditions) {
		collapsedDevs=list()
		for(element in 1:length(conditions)) {
			collapsedDevs[[element]]=apply(moduleEigengeneTable[,colnames(moduleEigengeneTable) %in%
				design$Sample[design[,refCondition]==conditions[[element]] ] ],1,sd)
		}
		sd1 = do.call(cbind, collapsedDevs)
		colnames(sd1)=conditions
		sd1
	}

	makeTable <- function(data, conditions){
		mean=collapseMeans(data, conditions)
		sd=collapseStdDev(data, conditions)
		table1=cbind(reshape2::melt(mean) %>% arrange(Var1),(reshape2::melt(sd) %>% arrange(Var1))$value)
		colnames(table1)=c("Module", "Condition", "Mean", "Stdev")
		table1
	}

	#prep tables for differential expression analysis
	table1=reshape2::melt(moduleEigengenes1) %>% arrange(Var1)
	colnames(table1)=c("Module","Sample","Expression")
	table1$Condition=design[match(table1$Sample,design$Sample),refCondition]

	table2=reshape2::melt(moduleEigengenes2) %>% arrange(Var1)
	colnames(table2)=c("Module","Sample","Expression")
	table2$Condition=design[match(table2$Sample,design$Sample),refCondition]

	#combinedTable=rbind(table1[table1$Module==paste0(prefix,"TG_01"),],
	#						table2[table2$Module==paste0(prefix,"WT_03"),])
	#combinedTable$Cont=parse_number(combinedTable$Condition)
	#combinedTable$RefCondition=(gsub(prefix,"",str_split_fixed(combinedTable$Module,"_",2)[,1]))
	#combinedTable$RefCondition[combinedTable$RefCondition=="TG"]=1
	#combinedTable$RefCondition[combinedTable$RefCondition=="WT"]=2

	#perform z-test of correlation coefficients using Fisher's z-transformation
	if(test=="zTestOfCorrelations"){
		#subset to only comparisons with at least one trait-associated module
		top_overlaps=top_overlaps[apply(cbind(lapply(top_overlaps$mod1,function(x)
										dataset1@trait$trait[gsub("ME","",dataset1@trait$Module)==x]!="None"),
									lapply(top_overlaps$mod2,function(x)
										dataset2@trait$trait[gsub("ME","",dataset2@trait$Module)==x]!="None")),1,any),]
		selectedTraits=list()
		cor1List=list()
		cor2List=list()
		pval=list()
		element=1
		for(row in 1:nrow(top_overlaps)){
			module1=top_overlaps$mod1[[row]]
			module2=top_overlaps$mod2[[row]]
			bestTrait1=dataset1@trait$trait[gsub("ME","",dataset1@trait$Module)==module1]
			bestTrait2=dataset2@trait$trait[gsub("ME","",dataset2@trait$Module)==module2]
			#if(bestTrait1!="None" | bestTrait2!="None"){
			bestP1=dataset1@trait[gsub("ME","",dataset1@trait$Module)==module1, paste0("p.value.", bestTrait1)]
			bestP2=dataset2@trait[gsub("ME","",dataset2@trait$Module)==module2, paste0("p.value.", bestTrait2)]
			if(is.null(bestP1)) bestP1=1
			if(is.null(bestP2)) bestP2=1
			if(bestP1<bestP2){
				selectedTrait=bestTrait1
				cor1=dataset1@trait[gsub("ME","",dataset1@trait$Module)==module1,
					if(grepl("^[[:digit:]]+", selectedTrait)) paste0("X", selectedTrait) else selectedTrait]
				cor2=dataset2@trait[gsub("ME","",dataset2@trait$Module)==module2,
					if(grepl("^[[:digit:]]+", selectedTrait)) paste0("X", selectedTrait) else selectedTrait]
			}
			if(bestP1>bestP2){
				selectedTrait=bestTrait2
				cor1=dataset1@trait[gsub("ME","",dataset1@trait$Module)==module1,
					if(grepl("^[[:digit:]]+", selectedTrait)) paste0("X", selectedTrait) else selectedTrait]
				cor2=dataset2@trait[gsub("ME","",dataset2@trait$Module)==module2,
					if(grepl("^[[:digit:]]+", selectedTrait)) paste0("X", selectedTrait) else selectedTrait]
			}
			selectedTraits[[element]]=selectedTrait
			cor1List[[element]]=cor1
			cor2List[[element]]=cor2
			n1=nrow(design[design[,testCondition]==name1,])
			n2=nrow(design[design[,testCondition]==name2,])
			pval[[element]]=paired.r(cor1, cor2, NULL, n1, n2, FALSE)$p
			element=element+1
		}
		p.adjust=p.adjust(pval)
		top_overlaps=cbind(top_overlaps, unlist(selectedTraits), unlist(cor1List), unlist(cor2List), unlist(pval), p.adjust)
		colnames(top_overlaps)[6:10]=c("selected.trait","cor1","cor2","zTest.pvalue", "zTest.padj")
		top_overlaps$Significance=FALSE
		top_overlaps$Significance[top_overlaps$zTest.padj<alpha]=TRUE
		top_overlaps=top_overlaps %>% arrange(zTest.padj)
	}

	#perform T-test comparing module eigengenes across conditions
	if(test=="t.test") {
		pval=list()
		element=1
		for(row in 1:nrow(top_overlaps)){
			for(condition in conditions){
				number1=sub("^0+","",gsub("\\D+", "", top_overlaps$mod1[[row]]))
				number2=sub("^0+","",gsub("\\D+", "", top_overlaps$mod2[[row]]))
				pval[[element]]=t.test(table1$value[sub("^0+","",sub(paste0(prefix,name1,"_"),"",table1$Var1))==number1 &
					table1$Var2 %in% design$Sample[design[,refCondition]==condition]],
					table2$value[sub("^0+","",sub(paste0(prefix,name2,"_"),"",table2$Var1))==number2 &
					table2$Var2 %in% design$Sample[design[,refCondition]==condition]])$p.value
				element=element+1
			}
		}
		pval=p.adjust(pval)
		top_overlaps=cbind(top_overlaps,t(matrix(pval,ncol=nrow(top_overlaps))))
		colnames(top_overlaps)[6:6+length(conditions)-1]=conditions
		top_overlaps$mostSignif=apply(top_overlaps[conditions],1,function(x) x[[which.min(x)]])
		top_overlaps$Significance=FALSE
		top_overlaps$Significance[top_overlaps$mostSignif<alpha]=TRUE
		top_overlaps=top_overlaps %>% arrange(mostSignif)
	}

	#plot all the interesting module expression changes as barplots
	meanAndSd1=makeTable(moduleEigengenes1,conditions)
	print(head(meanAndSd1))
	meanAndSd2=makeTable(moduleEigengenes2,conditions)
	print(head(meanAndSd2))
	dir.create(paste0(name1,"_vs_",name2))
	pdf(paste0(name1,"_vs_",name2,"/diffModuleExpressionBarplots.pdf"))
	for(row in 1:nrow(top_overlaps)) {
		number1=sub("^0+","",gsub("\\D+", "", top_overlaps$mod1[[row]]))
	    number2=sub("^0+","",gsub("\\D+", "", top_overlaps$mod2[[row]]))
		df=rbind(meanAndSd1[sub("^0+","",sub(paste0(prefix,name1,"_"),"", meanAndSd1$Module))==number1,],
			meanAndSd2[sub("^0+","",sub(paste0(prefix,name2,"_"),"", meanAndSd2$Module))==number2,])
		print(ggplot(data=df, aes(x=Condition, y=Mean, fill=Module)) +
			ggtitle(paste0("Overlap FDR = ",
			signif(top_overlaps$p.adj[[row]],2),
				" (",
				top_overlaps$overlap[[row]],
				"), z-test padj=",
			signif(top_overlaps$zTest.padj[[row]],2),
				" (",
			signif(top_overlaps$cor1[[row]],2),
				" vs ",
			signif(top_overlaps$cor2[[row]],2),")")) +
			geom_bar(stat="identity", color="black", position=position_dodge())+
			theme_linedraw()+
			theme(plot.title=element_text(hjust=0.5),panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
			geom_hline(yintercept=0, linetype="solid", color = "black", size=0.5)+
			geom_errorbar(aes(ymin=Mean-Stdev, ymax=Mean+Stdev), width=.2, position=position_dodge(.9)))
	}
	dev.off()

	#plot all the interesting module expression changes as scatterplots
	pdf(paste0(name1,"_vs_",name2,"/diffModuleExpressionScatterplots.pdf"))
	for(row in 1:nrow(top_overlaps)) {
		module1=top_overlaps$mod1[[row]]
		module2=top_overlaps$mod2[[row]]
		traitTable=makeTraitTable(sampleTable, refCondition)
		df=rbind(table1[table1$Module==paste0("ME",module1),], table2[table2$Module==paste0("ME",module2),])
		df$Condition=as.numeric(traitTable[match(df$Sample,traitTable$Sample), gsub("\\."," ", top_overlaps$selected.trait[[row]]) ])
		print(ggplot(data=df,
				aes(x=Condition, y=Expression, color=Module)) +
				geom_point() +
				geom_smooth(method='lm', se=FALSE) +
				theme_minimal() +
				scale_color_manual(values=c("red","blue")) +
				ggtitle(paste0("Overlap FDR = ",
				signif(top_overlaps$p.adj[[row]],2), " (", top_overlaps$overlap[[row]], "), z-test padj=",
				signif(top_overlaps$zTest.padj[[row]],2), " (",
				signif(top_overlaps$cor1[[row]],2), " vs ",
				signif(top_overlaps$cor2[[row]],2), ")")) +
				theme(plot.title=element_text(hjust=0.5)))
	}
	dev.off()

	fwrite(top_overlaps,paste0(name1,"_vs_",name2,"/top_overlaps_diffME.txt"),sep="\t")
	top_overlaps
}

calculateAverageExpression <- function(datExpr){
	cleanDatExpr=t(cleanDatExpr(datExpr))
	modules=sort(unique(datExpr$dynamicLabels))
	samples=(colnames(cleanDatExpr))
	averageExpressionOfModInSample=list()
	element=1
	for(module in modules){
		for(sample in samples){
			moduleGenes=datExpr$X[datExpr$dynamicLabels==module]
			averageExpressionOfModInSample[[element]]=mean(cleanDatExpr[rownames(cleanDatExpr) %in% moduleGenes,
														colnames(cleanDatExpr)==sample])
			element=element+1
		}
	}
	averageExpressionDf=as.data.frame(matrix(unlist(averageExpressionOfModInSample), ncol=length(modules), nrow=length(samples)))
	colnames(averageExpressionDf)=paste0("AE", modules)
	rownames(averageExpressionDf)=samples
	averageExpressionDf
}

geneExpBarPlot <- function(datExpr, gene){

  datExpr=cleanDatExpr(datExpr)
  expression=datExpr[,gene]
  df=data.frame(Sample=names(expression), Expression=expression, Condition=str_split_fixed(names(expression), "_", 2)[, 1])
  ggplot(df, aes(x= Sample, y=Expression, fill=Condition)) +
    geom_bar(stat="identity", color="black", position=position_dodge()) +
    scale_fill_manual(values=colors(length(unique(df$Condition))))+
    xlab("Sample") +
    ylab("Expression") +
    theme_classic()+
    theme(legend.position="none", axis.text.x=element_text(angle=45, hjust=1))

}
