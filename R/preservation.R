#' getPreservation
#'
#' Performs a network preservation analysis
#'
#' @param reference reference network of class WGCNA
#' @param test test network of class WGCNA
#' @param nPermutations number of permutations to perform
#' @param write write to file?
#' 
#' @author Dario Tommasini
#'
#' @import WGCNA
#' @import stringr
#' @export
getPreservation <- function(reference, test, nPermutations=100, write=FALSE) {

	reference=reference@datExpr
	test=test@datExpr

	#clean up reference datExpr
	name1=str_split_fixed(reference$dynamicLabels, "_", 2)[1,1]
	name2=str_split_fixed(test$dynamicLabels, "_", 2)[1,1]
	referenceExpr = cleanDatExpr(reference, checkGenesSamples=T)

	#clean up test datExpr
	testExpr = cleanDatExpr(test, checkGenesSamples=T)

	referenceModules=reference$dynamicLabels[match(colnames(referenceExpr), reference$X)]

	#testModules=test$dynamicLabels[match(colnames(testExpr), test$X)]

	#calculate perservation Z-scores
	multiExpr2 = list(A = list(data = referenceExpr), B = list(data = testExpr));
	modLabels2 = list(A = referenceModules) #, B = testModules)
	system.time( {
		mp2 = modulePreservation(multiExpr2, modLabels2,
              					referenceNetworks = 1,
                        nPermutations = nPermutations,
                        #maxGoldModuleSize = left as default
                        #maxModuleSize = left as default,
                        randomSeed = 1,
                        quickCor = 0,
                        verbose = 3,
              					parallelCalculation=TRUE)
		} );
	ref = 1
	test = 2

	if(write) write.csv(mp2$preservation$Z[[1]][[2]], paste0(name1, "_in_", name2,".csv"))

	as.data.frame(mp2$preservation$Z[[1]][[2]])
}

#' Preservation scatterplot
#'
#' A plotting  function that draws a scatterplot of preservation scores between
#' two WGCNAobjects
#'
#' @param comparison a data.frame resulting from a call to computeOverlapsFromWGCNA
#' @param dataset1 an object of class WGCNAobject to compare with dataset2
#' @param dataset2 an object of class WGCNAobject to compare with dataset1
#' @param alphaLevel alpha level of significance
#' @param outliers leave outlier modules? By default these are removed
#'
#' @author Dario Tommasini
#'
#' @import ggplot2
#' @import stringr
#' @import ggrepel
#' @export
preservationComparisonPlot <- function(comparison, dataset1, dataset2, alphaLevel, outliers=FALSE){
	name1=str_split_fixed(rownames(comparison$mod1Preservation[rownames(comparison$mod1Preservation) != "gold",]),"_",2)[,1][[1]]
	name2=str_split_fixed(rownames(comparison$mod2Preservation[rownames(comparison$mod2Preservation) != "gold",]),"_",2)[,1][[1]]

	#merge perservation statistic and trait correlation of each module into one df, then make scatterplot
	input1=dataset1@trait
	input1$Zsum=comparison$mod1Preservation$Zsummary.pres[match(input1$Module, rownames(comparison$mod1Preservation))]
	if(!outliers) input1=input1[! input1$Module %in% dataset1@outlierModules,]
	input2=dataset2@trait
        input2$Zsum=comparison$mod2Preservation$Zsummary.pres[match(input2$Module, rownames(comparison$mod2Preservation))]
        if(!outliers) input2=input2[! input2$Module %in% dataset2@outlierModules,]
	myColors=colors(length(sort(unique(c(input1$trait, input2$trait)))))
	names(myColors)=sort(unique(c(input1$trait, input2$trait)))
	myColors1=myColors[match(sort(unique(input1$trait)), names(myColors))]
	myColors2=myColors[match(sort(unique(input2$trait)), names(myColors))]

	dataset1Plot <- ggplot(input1, aes(log10Pvalue, Zsum, fill=trait)) +
			geom_point(shape=21, size=2)+
			theme_classic()+
			scale_fill_manual(values=myColors1)+
			labs(title=paste0(name1, " in ", name2), y = "Preservation (z-summary)", x="Trait association (-log10 p-value)")+
		        theme(plot.title=element_text(hjust=0.5))+
			geom_text_repel(data=input1,
				aes(label=paste0(gsub("^0+","",substr(Module, nchar(Module[[1]])-1, nchar(Module[[1]]))))))+
			geom_hline(yintercept=(2), linetype="dashed", color = "red", size=1)+
			geom_hline(yintercept=(10), linetype="dashed", color = "black", size=1)+
			geom_vline(xintercept=(-log10(alphaLevel)), linetype="dashed", color = "red", size=1)

	dataset2Plot <-ggplot(input2, aes(log10Pvalue, Zsum, fill=trait)) +
		geom_point(shape=21, size=2)+
			# theme(panel.grid.major = element_blank(),
				# panel.grid.minor = element_blank(),
				# panel.background = element_blank(),
				# axis.line = element_line(colour = "black"),
				# plot.title=element_text(hjust=0.5))+
			theme_classic()+
			scale_fill_manual(values=myColors2)+
			labs(title=paste0(name2, " in ", name1), y = "Preservation (z-summary)", x="Trait association (-log10 p-value)")+
	        	theme(plot.title=element_text(hjust=0.5))+
			geom_text_repel(data=input2,
				aes(label=paste0(gsub("^0+","", substr(Module, nchar(Module[[1]])-1, nchar(Module[[1]]))))))+
			geom_hline(yintercept=(2), linetype="dashed", color = "red", size=1)+
			geom_hline(yintercept=(10), linetype="dashed", color = "black", size=1)+
			geom_vline(xintercept=(-log10(alphaLevel)), linetype="dashed", color = "red", size=1)

	print(dataset1Plot | dataset2Plot)
}

#' Preservation comparisons
#'
#' A high level function that performs a perservation comparison between two
#' WGCNAobjects in a WGCNAlist, usually supplied by iterate function
#'
#' @param comparisonList a list passed by the iterate function
#' @param WGCNAlist list of objects of type WGCNAobject
#' @param first index of first WGCNAobject
#' @param second index of second WGCNAobject
#' @param element element position in the comparison list (passed by iterate function)
#' @param plot generate plots?
#' @param write write results to file? 
#' @param alphaLevel alpha level of significance for module-trait correlation
#' @param nPermutations number of permutations, defaults to 100
#'
#' @author Dario Tommasini
#'
#' @export
preservationComparisons <- function(comparisonList, WGCNAlist, first, second, element, plot=FALSE, write=FALSE, alphaLevel=get("alphaLevel", envir = parent.frame()), nPermutations=100){
	comparisonList[[element]]=list()
	comparisonList[[element]]=append(comparisonList[[element]],
								list(getPreservation(WGCNAlist[[first]],
									WGCNAlist[[second]], write=write, nPermutations=nPermutations)))
	comparisonList[[element]]=append(comparisonList[[element]],
								list(getPreservation(WGCNAlist[[second]],
									WGCNAlist[[first]], write=write, nPermutations=nPermutations)))
	names(comparisonList)[[element]]=paste0(names(WGCNAlist)[[first]], "_vs_", names(WGCNAlist)[[second]])
	names(comparisonList[[element]])[[1]]="mod1Preservation"
	names(comparisonList[[element]])[[2]]="mod2Preservation"
	if(plot){
		preservationComparisonPlot(comparisonList[[element]],
								WGCNAlist[[first]],
								WGCNAlist[[second]],
								alphaLevel)
	}
	comparisonList
}

coexpressionLineGraph <- function(datExpr, nDiseaseSamples, nWTSamples, splitBy=1, fontSize=2.15, colors=NULL){

	nGenes=ncol(datExpr)
	scaled=scale(datExpr)

	if(splitBy>0){
		for(column in 1:ncol(scaled)){
			scaled[, column]= scaled[,column]+splitBy*column
		}
	}

	ggplot(reshape2::melt(as.matrix(scaled)), aes(x = Var1, y = value, group= Var2, color=Var2, label=Var2)) +
				#geom_line(aes(color=factor(rep(c(rep("healthy", nWTSamples), rep("disease", nDiseaseSamples)), nGenes), levels=c("healthy", "disease")))) +
				geom_line()+
				#scale_color_hue(name = "Gene", l=70, c=30) +
				{if(is.null(colors)) scale_colour_manual(values=colors(nGenes))}+
				{if(!is.null(colors)) scale_colour_manual(values=colors)}+
				labs(y="Scaled expression", x="Samples") +
				theme_classic() +
				coord_cartesian(clip="off")+
				#scale_x_discrete(limits=c(0, 40))+
				annotate("text", x=nWTSamples+nDiseaseSamples+0.5, y=-0.5, hjust=0, vjust=0,
					label=paste0(rev(colnames(datExpr)), collapse='', sep='\n'), size=fontSize)+
				#geom_text_repel()+
				geom_vline(xintercept= nDiseaseSamples+1, linetype='dashed') +
				theme(legend.position="none", axis.text.x=element_text(angle=90, vjust=0.25, hjust=1),
					axis.text.y=element_blank(), axis.ticks=element_blank(),
					legend.key.width= unit(0.001, 'mm'), plot.margin = unit(c(1,1,1,1), "cm"))
}

correlationComparisonBoxplot <- function(diseaseDatExpr, healthyDatExpr, geneList, label=F, method="pearson"){
	name1=str_split_fixed(colnames(diseaseDatExpr)[[2]], "_", 2)[,1]
	name2=str_split_fixed(colnames(healthyDatExpr)[[2]], "_", 2)[,1]

	geneList=geneList[geneList %in% diseaseDatExpr$X & geneList %in% healthyDatExpr$X ]
	diseaseMatrix= cleanDatExpr(diseaseDatExpr[match(geneList, diseaseDatExpr$X),])
	healthyMatrix= cleanDatExpr(healthyDatExpr[match(geneList, healthyDatExpr$X),])

	diseaseCorMatrix=cor(diseaseMatrix, method=method)
	diseaseCor=diseaseCorMatrix[lower.tri(diseaseCorMatrix)]

	healthyCorMatrix=cor(healthyMatrix, method=method)
	healthyCor=healthyCorMatrix[lower.tri(healthyCorMatrix)]

	#boxplot(diseaseCor, healthyCor)
	indices <- which(lower.tri(diseaseCorMatrix) == TRUE, arr.ind=T)
	names1 <- geneList[indices[,1]]
	names2 <- geneList[indices[,2]]
	pair=paste0(names1,"/", names2)

	melted=as.data.frame(rbind(cbind(pair, diseaseCor, rep(name1, length(diseaseCor))),
						cbind(pair, healthyCor, rep(name2, length(healthyCor)))))
	colnames(melted)=c("Pair", "Correlation", "Status")
	melted$Correlation=as.numeric(melted$Correlation)

	print(t.test(diseaseCor))
	print(t.test(healthyCor))
	print(t.test(diseaseCor, healthyCor))

	ggplot(melted, aes(x=Status, y=Correlation, fill=Status)) +
		{if(label) geom_text_repel()} +
		theme_classic() +
		scale_y_continuous(position="right")+
		theme(legend.position="none") +
		scale_fill_manual(values=c("magenta", "cyan")) +
		stat_boxplot(geom = 'errorbar', lwd=1, width = 0.3, coef = 3) +
		geom_boxplot(notch = TRUE, outlier.shape=NA)
		#geom_dotplot(binaxis='y', stackdir='center', stackratio=1.5,
		#			binwidth=0.01, dotsize=1)
}

correlationComparisonHeatmaps <- function(diseaseDatExpr, healthyDatExpr, geneList, label=F, method="pearson", alphaLevel=0.05, p.adjust=T, z.score.limit=3){
	geneList=geneList[geneList %in% diseaseDatExpr$X & geneList %in% healthyDatExpr$X ]
	diseaseMatrix= cleanDatExpr(diseaseDatExpr[match(geneList, diseaseDatExpr$X),])
	healthyMatrix= cleanDatExpr(healthyDatExpr[match(geneList, healthyDatExpr$X),])

	diseaseCorMatrix=cor(diseaseMatrix, method=method)
	diseaseCor=diseaseCorMatrix[upper.tri(diseaseCorMatrix)]

	healthyCorMatrix=cor(healthyMatrix, method=method)
	healthyCor=healthyCorMatrix[upper.tri(healthyCorMatrix)]

	#boxplot(diseaseCor, healthyCor)
	indices <- which(upper.tri(diseaseCorMatrix) == TRUE, arr.ind=T)
	names1 <- geneList[indices[,1]]
	names2 <- geneList[indices[,2]]
	pair=paste0(names1,"/", names2)

	melted=as.data.frame(rbind(cbind(pair, diseaseCor, rep("Disease", length(diseaseCor))),
						cbind(pair, healthyCor, rep("Healthy", length(healthyCor)))))
	colnames(melted)=c("Pair", "Correlation", "Status")
	melted$Correlation=as.numeric(melted$Correlation)

	dc=diffCoexpression(cbind(t(diseaseMatrix), t(healthyMatrix)),
		c(rep(1, nrow(diseaseMatrix)),rep(2, nrow(healthyMatrix))), plot=F, FDR.threshold=alphaLevel)
	z_scores=dc[[1]]
	z_scores[z_scores > z.score.limit]=z.score.limit
	z_scores[z_scores < -z.score.limit]=-z.score.limit
	z.table=reshape2::melt(z_scores)
	z.table$p.value=as.list(dc[[2]])
	z.table$p.adj=as.list(dc[[3]])
	z.table$signif=FALSE
	if(p.adjust==T) {
		z.table$signif[z.table$p.adj<alphaLevel]=TRUE
	} else {
		z.table$signif[z.table$p.value<alphaLevel]=TRUE
	}

	meltedDs=reshape2::melt(diseaseCorMatrix)
	meltedDs$signif=z.table$signif
	meltedWt=reshape2::melt(healthyCorMatrix)
	meltedWt$signif=z.table$signif

	diseasePlot <- ggplot(meltedDs, aes(x=Var1, y=Var2, fill=value)) +
					geom_tile(aes(color = signif), size=0.1) +
					geom_vline(xintercept=c(13.5, 24.5))+
					geom_hline(yintercept=c(13.5, 24.5))+
					scale_fill_gradient2(name=NULL, low = "red", mid="white", high = "gold",
							limits=c(-1, 1), na.value="grey",
							guide = guide_colorbar(frame.colour = "black", frame.linewidth = 1,
							ticks.linewidth=1, ticks.colour = "black")) +
					scale_color_manual(values=c("white", "black"), guide = "none")+
					theme_classic()+
					theme(axis.text.x = element_text(angle = 90, vjust=(0.25), hjust=1),
							axis.title=element_blank(), axis.ticks=element_blank(),
							legend.key.width = unit(4, 'mm'), legend.key.height = unit(20, 'mm')) +
					{if(!label) theme(axis.text.x=element_blank(), axis.text.y=element_blank())}+
					coord_fixed()
	healthyPlot <- ggplot(meltedWt, aes(x=Var1, y=Var2, fill=value)) +
					geom_tile(aes(color = signif), size=0.1) +
					geom_vline(xintercept=c(13.5,24.5))+
					geom_hline(yintercept=c(13.5,24.5))+
					scale_fill_gradient2(name=NULL, low = "red", mid="white", high = "gold",
							limits=c(-1, 1), na.value="grey",
							guide = guide_colorbar(frame.colour = "black", frame.linewidth = 1,
							ticks.linewidth=1, ticks.colour = "black")) +
					scale_color_manual(values=c("white", "black"), guide = "none")+
					theme_classic()+
					theme(axis.text.x = element_text(angle = 90, vjust=(0.25), hjust=1),
							axis.title=element_blank(), axis.ticks=element_blank(),
							legend.key.width = unit(4, 'mm'), legend.key.height = unit(20, 'mm')) +
					{if(!label) theme(axis.text.x=element_blank(), axis.text.y=element_blank())}+
					coord_fixed(clip="off")
	dcPlot <- ggplot(z.table, aes(x=Var1, y=Var2, fill=value)) +
					geom_tile(aes(color = signif), size=0.1) +
					geom_vline(xintercept=c(13.5, 24.5), color="white")+
					geom_hline(yintercept=c(13.5, 24.5), color="white")+
					scale_fill_gradient2(name=NULL, low = "cyan", mid="black", high = "magenta",
							limits=c(-z.score.limit, z.score.limit), na.value="grey",
							guide = guide_colorbar(frame.colour = "black", frame.linewidth = 1,
							ticks.linewidth=1, ticks.colour = "black")) +
					scale_color_manual(values=c("black", "white"), guide = "none")+
					theme_classic()+
					theme(axis.text.x = element_text(angle = 90, vjust=(0.25), hjust=1),
							axis.title=element_blank(), axis.ticks=element_blank(),
							legend.key.width = unit(4, 'mm'), legend.key.height = unit(20, 'mm')) +
					{if(!label) theme(axis.text.x=element_blank(), axis.text.y=element_blank())}+
					coord_fixed(clip="off")
	boxPlot <- correlationComparisonBoxplot(diseaseDatExpr, healthyDatExpr, geneList)


	#ggarrange(diseasePlot, healthyPlot, dcPlot, boxPlot, ncol=4, align="hv", widths=c(3, 3, 3, 1), axis = "bt")
	plot_grid(diseasePlot, healthyPlot, dcPlot, boxPlot, ncol=4, align = "hv",
		axis = "tb", rel_widths = c(2, 2, 2, 1))
}

#' Differential co-expresison analysis 
#'
#' Performs a differential co-expression ananlysis given an expression data.frame 
#' and a conditions vector
#'
#' @param datExpr a data.frame containing expression values
#' @param conditions a vector containing conditions for the samples
#' @param geneList vector of genes, will use all genes if NULL (default)
#' @param plot plot a network?
#' @param method either "pearson" or "spearman" 
#' @param removeFreeNodes remove free nodes from network?
#' @param labelSize label size
#' @param labelDist distance from labels to nodes
#' @param shape shape of nodes
#' @param degreeForSize should node size correspond to degree?
#' @param label label nodes?
#' @param onlyPositive only draw positive correlations?
#' @param z.threshold z-score threshold
#' @param FDR.threshold FDR threshold
#' @param nodeSize size of node
#' 
#' @author Dario Tommasini
#'
#' @import dcanr
#' @export
diffCoexpression <- function(datExpr, conditions, geneList=NULL, plot=F, method="pearson", removeFreeNodes=TRUE, labelSize=0.5, labelDist=0,
						shape="circle", degreeForSize=F, label=F, onlyPositive=F, z.threshold=NULL, FDR.threshold=0.05, nodeSize=3){

	if(!is.null(geneList)) datExpr=datExpr[rownames(datExpr) %in% geneList,]
	cor1= cor(t(datExpr[, conditions==1]), method=method)
	cor2= cor(t(datExpr[, conditions==2]), method=method)
	z_scores <- dcScore(datExpr, conditions, dc.method = 'zscore', cor.method = method)
	raw_p <- dcTest(z_scores, datExpr, conditions)
	adj_p <- dcAdjust(raw_p, f = p.adjust, method = 'fdr')

	geneList=rownames(datExpr)
	indices <- which(upper.tri(z_scores) == TRUE, arr.ind=T)
	names1 <- geneList[indices[,1]]
	names2 <- geneList[indices[,2]]
	summaryDf=	data.frame(gene1=names1, gene2=names2, cor1= cor1[upper.tri(cor1)], cor2=cor2[upper.tri(cor2)],
					z.score=z_scores[upper.tri(z_scores)], p.value=raw_p[upper.tri(raw_p)], p.adj=adj_p[upper.tri(adj_p)])

	adj_mat= z_scores
	adj_mat[diag(adj_mat)]=0
	if(is.null(z.threshold)) z.threshold=min(abs(na.omit(adj_mat[adj_p<FDR.threshold])))
	graph=graph_from_adjacency_matrix(adj_mat, mode="undirected", weighted=T)
	graph=simplify(graph)
	graph <- delete.edges(graph, which(abs(E(graph)$weight) < z.threshold))
	if(onlyPositive) graph <- delete.edges(graph, which(E(graph)$weight < 0))
	if(removeFreeNodes) graph <- delete.vertices(simplify(graph), degree(graph)==0)
	ealpha=rescale(abs(E(graph)$weight), to=c(0.1, 0.5))
	ecol=lapply(ealpha, function(x) rgb(0, 0, 1, x))
	ecol[which(E(graph)$weight < 0)]=rgb(1, 0, 1, ealpha[which(E(graph)$weight < 0)])
	if(degreeForSize) {
		V(graph)$size = log(degree(graph, mode = "out") + nodeSize) * 2
	} else {
		V(graph)$size=nodeSize
	}

	if(!label) V(graph)$name=NA
	#ecol=lapply(E(graph)$weight, function(x) rgb(1-x, 1-x, 1-x))
	#ewidth=rescale(E(graph)$weight, from=c(0,1), to=c(0, 1.5))

	if(plot){
		plot(graph, layout=layout_with_fr, edge.width=1.5, vertex.shape=shape,
		vertex.color="white", vertex.label.dist=labelDist, vertex.label.color="black",
		#vertex.size=(strwidth(V(graph)$name)+strwidth("oo"))*100,
		vertex.label.cex=labelSize, edge.color=unlist(ecol))
	} else {
		list(z_scores, raw_p, adj_p, summaryDf %>% arrange(p.adj))
	}
}
