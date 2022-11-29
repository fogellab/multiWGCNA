#' computeOverlapsFromWGCNA
#'
#' Computes overlap between the modules of two objects of class WGCNA
#'
#' @param dataset1 an object of class WGCNA to compare with dataset2
#' @param dataset2 an object of class WGCNA to compare with dataset1
#' @param convertSymbols1 convert symbols for first WGCNA
#' @param convertSymbols2 convert symbols for second WGCNA
#' @author Dario Tommasini
#'
#' @import biomaRt
#' @export
computeOverlapsFromWGCNA <- function(dataset1, dataset2, convertSymbols1=F, convertSymbols2=F) {
	datExpr1= dataset1@datExpr
	datExpr2= dataset2@datExpr
	treatDat <- datExpr1
	controlDat <- datExpr2

	#computation of overlaps
	mod1=list()
	mod2=list()
	mod1_size=list()
	mod2_size=list()
	genes=list()
	pval=list()
	element=1

	if(convertSymbols1 | convertSymbols2){
		human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
		mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
	}

	for(treatment in sort(unique(treatDat$dynamicLabels))){
	        for(control in sort(unique(controlDat$dynamicLabels))){
	                mod1[[element]]=treatment
	                mod2[[element]]=control
	                if(convertSymbols1){
	                	mod1Genes=toupper(treatDat$X[treatDat$dynamicLabels==treatment])
	                	mod1Genes=toupper(unique(getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = mod1Genes,
	                		mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)))
	                } else {
	                	mod1Genes=toupper(treatDat$X[treatDat$dynamicLabels==treatment])
	                }
	                if(convertSymbols2){
	                		mod2Genes=toupper(controlDat$X[controlDat$dynamicLabels==control])
	                		mod2Genes=toupper(unique(getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = mod2Genes,
	                			mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)))
	                } else {
	                	mod2Genes=toupper(controlDat$X[controlDat$dynamicLabels==control])
	                }
	                mod1_size[[element]]=length(mod1Genes)
			mod2_size[[element]]=length(mod2Genes)
			genes[[element]]=length(intersect(mod1Genes, mod2Genes))
	                pval[[element]]=phyper(genes[[element]]-1,
	                		mod1_size[[element]],
	                		23000-mod1_size[[element]],
	                		mod2_size[[element]],
	                		lower.tail=FALSE, log.p=FALSE)
	                if(pval[[element]]==0) pval[[element]]=.Machine$double.xmin
	                element=element+1
	        }
	}
	data.frame(mod1=unlist(mod1), mod2=unlist(mod2), mod1.size=unlist(mod1_size),
		   	mod2.size=unlist(mod2_size), overlap=unlist(genes),
			p.value=unlist(pval), p.adj=unlist(p.adjust(pval, method='fdr')))
}

#' Plots an expression profile for a module
#'
#' A plotting function that returns a heatmap and barplot for a module
#'
#' @param overlapDf a data.frame resulting from a call to computeOverlapsFromWGCNA
#' @param dataset1 an object of class WGCNAobject to compare with dataset2
#' @param dataset2 an object of class WGCNAobject to compare with dataset1
#'
#' @author Dario Tommasini
#'
#' @import ggplot2
#' @import ggalluvial
#' @import stringr
#' @import ggrepel
#' @import ggpubr
#' @export
moduleComparisonPlot <- function(overlapDf, dataset1, dataset2) {
	data=overlapDf
	name1=str_split_fixed(data$mod1,"_",2)[,1][[1]]
	name2=str_split_fixed(data$mod2,"_",2)[,1][[1]]
	data$mod1=str_split_fixed(data$mod1,"_",2)[,2]
	data$mod2=paste0(str_split_fixed(data$mod2,"_",2)[,2]," ") #add spaces to avoid ambiguous plot labels
	categories=unique(c(dataset1@trait$trait, dataset2@trait$trait))
	palette=colors(length(categories), random=F)
	colors=palette[match(c(rev(dataset1@trait$trait), rev(dataset2@trait$trait)), categories)]
	totalGenes=sum(data$overlap)

 	#draw heatmap
 	heatmap <- moduleToModuleHeatmap(overlapDf)

	#draw flow plot
	flowPlot <- ggplot(data, aes(y = overlap, axis1 = mod1, axis3 = mod2)) +
		geom_flow(aes(fill = -log10(p.adj), alpha = -log10(p.adj)), width = .2, curve_type = "linear") +
		ylab("Genes") +
		scale_fill_gradient(name="Overlap", low="cyan", high="magenta")+
		theme(axis.ticks.x=element_blank(), legend.position = "none", axis.text.x=element_text(size=12),
			axis.title.y=element_text(size=15),
		panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
 		geom_stratum(width = .2, fill = colors) +
 		#scale_colour_identity('Trait cor', labels = categories, breaks=palette, guide = 'legend')+
		annotate("text", x=2.2, y=(1:length(categories))*(totalGenes/3)/length(categories),
			label=categories, vjust=0, hjust=0, color=palette)+
		geom_text(stat = "stratum", aes(label = gsub("^0+", "", after_stat(stratum))), size = 3, min.y=200) +
		coord_cartesian(xlim = c(0.9, 2.5), clip = "off") +
		#ggfittext::geom_fit_text(stat = "stratum", aes(label = after_stat(stratum)), min.size = 1) +
 		scale_x_discrete(expand=c(0,0), limits = c("Mod1", "Mod2"), labels=c(name1,name2))

	print(ggarrange(flowPlot, heatmap, nrow=1, ncol=2, widths=c(1,2.5)))
}


#for study of module member allocation over a continuous trait ie time
continuousFlowPlot <- function(WGCNAlist){

	allGenes=unique(ulist(lapply(WGCNAlist, function(WGCNAobject) WGCNAobject@datExpr$X)))
	df=data.frame(Genes=allGenes)
	for(WGCNA in WGCNAlist){
		df=cbind(df, WGCNA@datExpr$dynamicLabels[match(allGenes, WGCNA@datExpr$X)])
	}
	names(df)=c("Genes",names(WGCNAlist))
	head(df)

	combinations=apply(df[,2:5], 1, function(x) x)
	geneOverlap=list()
	for(element in 1:ncol(combinations)){
		geneOverlap[[element]]=nrow(df[apply(df[,2:5], 1, function(x) all(x==combinations[,element])),])
	}
	data=as.data.frame(cbind(geneOverlap, t(combinations)))
	colnames(data)=c("geneOverlap", paste0("mod", 1:length(WGCNAlist)))

	uniqueData=unique(data)

	mod1=(sort(unique(WGCNAlist[[1]]@datExpr$dynamicLabels)))
	mod2=(sort(unique(WGCNAlist[[2]]@datExpr$dynamicLabels)))
	mod3=(sort(unique(WGCNAlist[[3]]@datExpr$dynamicLabels)))
	mod4=(sort(unique(WGCNAlist[[4]]@datExpr$dynamicLabels)))

	n1=length(mod1)
	n2=length(mod2)
	n3=length(mod3)
	n4=length(mod4)

	finalData=data.frame(mod1=unlist(lapply(mod1, function(x) rep(x, n2*n3*n4))),
		mod2=rep(unlist(lapply(mod2, function(x) rep(x, n3*n4))), n1),
		mod3=rep(unlist(lapply(mod3, function(x) rep(x, n4))), n1*n2),
		mod4=rep(mod4, n1*n2*n3), stringsAsFactors=FALSE)

	finalData$identifier=apply(finalData, 1, function(x) paste0(x, sep="_"))


	ggplot(uniqueSortedData, aes(y = geneOverlap, axis1 = mod1, axis2 = mod2, axis3 = mod3, axis4=mod4)) +
		geom_flow(aes(fill = geneOverlap), width = .4, curve_type = "linear") + ylab("Genes") +
		scale_fill_manual(values=c("white","red")) +
		theme(axis.ticks.x=element_blank(),legend.position = "none", axis.text.x=element_text(size=12),
			axis.title.y=element_text(size=15), panel.background = element_blank(),
			panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
		geom_stratum(width = .4) +
		geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3, min.y=250) +
		scale_x_discrete(expand=c(0,0), limits = c("Mod2", "Mod3", "Mod1", "Mod4"), labels=names(WGCNAlist))

	if(length(WGCNAlist)==4){
		mod1=list()
		mod2=list()
		mod3=list()
		mod4=list()
		genes=list()
		element=1
		for(one in sort(unique(WGCNAlist[[1]]@datExpr$dynamicLabels))){
	        for(two in sort(unique(WGCNAlist[[2]]@datExpr$dynamicLabels))){
                for(three in sort(unique(WGCNAlist[[3]]@datExpr$dynamicLabels))){
               		for(four in sort(unique(WGCNAlist[[4]]@datExpr$dynamicLabels))){
                		mod1[[element]]=one
                        mod2[[element]]=two
                        mod3[[element]]=three
                        mod4[[element]]=four
                        #genes[[element]]=nrow(df[df[,2]==one & df[,3]==two & df[,4]==three & df[,5]==four,])
                        element=element+1
                	}
                }
			}
		}

	data=data.frame(genes=unlist(genes), mod1=unlist(mod1), mod2=unlist(mod2), mod3=unlist(mod3), mod4=unlist(mod4))
	data$ofInterest=FALSE
	data$ofInterest[as.character(data$mod3)==as.character(args[7])]=TRUE



	} else {
		error("only WGCNAlist length of four implemented at the moment")
	}


}

#' Module to module heatmap
#'
#' Returns a heatmap where color corresponds to FDR-adjusted overlap (hypergeometric test) and the label corresponds to the number of overlapping genes
#'
#' @param comparisonDf the output of computeOverlapFromWGCNA
#' @param dataset1 WGCNAobject for dataset 1
#' @param dataset2 WGCNAobject for dataset 2
#' @param trait1 subset to modules correlated to this trait for dataset 1
#' @param trait2 subset to modules correlated to this trait for dataset 2
#' @param list1 subset to this list of modules for dataset 1
#' @param list2 subset to this list of modules for dataset 2
#' @param filterByTrait only plot for modules that correlate with some trait?
#' @param alphaLevel the alpha level of significance for module-trait correlation, defaults to 0.05
#'
#' @export
moduleToModuleHeatmap <- function(comparisonDf, dataset1=NULL, dataset2=NULL, trait1=NULL, trait2=NULL, list1=NULL, list2=NULL, filterByTrait=FALSE, alphaLevel=0.05){

	#filter by trait if desired
	if(!is.null(dataset1)){
		if(!is.null(list1)){
			mod1ToKeep=list1
			comparisonDf=comparisonDf[comparisonDf$mod1 %in% mod1ToKeep,]
		} else if(!is.null(trait1)){
			association1=paste0("p.value.", trait1)
			index=which(colnames(dataset1@trait)==association1)
			mod1ToKeep=gsub("ME", "", dataset1@trait$Module[dataset1@trait[,index]<alphaLevel])
			comparisonDf=comparisonDf[comparisonDf$mod1 %in% mod1ToKeep,]
		} else if(filterByTrait){
			mod1ToKeep=gsub("ME", "", dataset1@trait$Module[dataset1@trait$p.value.interest_trait_code<alphaLevel])
			comparisonDf=comparisonDf[comparisonDf$mod1 %in% mod1ToKeep,]
		}
	}
	if(!is.null(dataset2)){
		if(!is.null(list2)){
			mod2ToKeep=list2
			comparisonDf=comparisonDf[comparisonDf$mod2 %in% mod2ToKeep,]
		} else if(!is.null(trait2)){
			association2=paste0("p.value.", trait2)
			index=which(colnames(dataset2@trait)==association2)
			mod2ToKeep=gsub("ME", "", dataset2@trait$Module[dataset2@trait[,index]<alphaLevel])
			comparisonDf=comparisonDf[comparisonDf$mod2 %in% mod2ToKeep,]
		} else if(filterByTrait){
			mod2ToKeep=gsub("ME", "", dataset2@trait$Module[dataset2@trait$p.value.interest_trait_code<alphaLevel])
			comparisonDf=comparisonDf[comparisonDf$mod2 %in% mod2ToKeep,]
		}
	}

	columns=str_split_fixed(comparisonDf$mod1,"_", 2)[,1][[1]]
	rows=str_split_fixed(comparisonDf$mod2,"_",2)[,1][[1]]
	comparisonDf$mod1=gsub("_"," ", gsub("^0+", "", str_split_fixed(comparisonDf$mod1,"_",2)[,2]))
	comparisonDf$mod2=gsub("_"," ", gsub("^0+", "", str_split_fixed(comparisonDf$mod2,"_",2)[,2]))

	ggplot(comparisonDf, aes(x = factor(mod1, levels=(unique(mod1))),
 				y = factor(mod2, levels=rev(unique(mod2))),
 				fill = (-log10(p.adj)), label = overlap)) +
				geom_tile(color = "black") +
				scale_fill_gradient(name="-log10(FDR)", low = "white",
					high = "red", na.value="red", #limits=c(0, 50),
					guide = guide_colorbar(frame.colour = "black", frame.linewidth = 1,
							ticks.linewidth=1, ticks.colour = "black")) +
				geom_text(aes(label = overlap), color = "black", size=2) +
				labs(x=columns, y=rows) +
				theme_classic()+
				theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=(0.5)),
					axis.text.y = element_text(hjust=1),
					panel.background=element_blank(), axis.line=element_blank(),
					axis.ticks=element_blank(), legend.key.size=unit(4, "mm")) +
				#geom_fit_text(reflow = TRUE)+
				coord_fixed()

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
#' @import ggplot2
#' @import dplyr
#' @export
bidirectionalBestMatches <- function(overlapDf, plot=TRUE){
  comparison=overlapDf
	name1=str_split_fixed(comparison$overlap$mod1,"_",2)[,1][[1]]
	name2=str_split_fixed(comparison$overlap$mod2,"_",2)[,1][[1]]
	comparison$overlap$mod1=str_split_fixed(comparison$overlap$mod1,"_",2)[,2]
	comparison$overlap$mod2=str_split_fixed(comparison$overlap$mod2,"_",2)[,2]
	columns=(unique(comparison$overlap$mod1))
	rows=(unique(comparison$overlap$mod2))
	padj_matrix=matrix(unlist(comparison$overlap$p.adj), ncol=length(columns),nrow=length(rows))
	colnames(padj_matrix)=columns
	rownames(padj_matrix)=rows
	bestMod1=list()
	bestMod2=list()
	p.adj=list()
	count=1
	for(column in 1:length(columns)){
		for(row in 1:length(rows)){
			if(which.min(padj_matrix[,column])==row & which.min(padj_matrix[row,])==column){
				bestMod1[[count]]=columns[[column]]
				bestMod2[[count]]=rows[[row]]
				p.adj[[count]]=min(padj_matrix[,column])
				count=count+1
			}
		}
	}
	bestMatches=data.frame(mod1=unlist(bestMod1), mod2=unlist(bestMod2), p.adj=unlist(p.adj))
	bestMatchesSorted=bestMatches %>% arrange(p.adj)
	subset_padj_matrix=padj_matrix[bestMatchesSorted$mod2,bestMatchesSorted$mod1]
	subsetComparisonDf=comparison$overlap[comparison$overlap$mod1 %in% bestMatches$mod1 &
											comparison$overlap$mod2 %in% bestMatches$mod2,]
	if(plot) {
		print(ggplot(subsetComparisonDf, aes(x = factor(mod1, levels=colnames(subset_padj_matrix)),
 				y = factor(mod2, levels=rev(rownames(subset_padj_matrix))),
 				fill = (-log10(p.adj)))) +
				geom_tile(color = "black") +
				scale_fill_gradient(low = "white", high = "red", na.value="red",
				                    guide = guide_colorbar(frame.colour = "black",
				                    frame.linewidth = 1, ticks.linewidth=1,
				                    ticks.colour = "black")) +
				geom_text(aes(label = overlap), color = "black") +
				labs(x=name1, y=name2) +
				theme(axis.text.x = element_text(angle = 90, vjust=(0.5)), panel.background=element_blank())+
				coord_fixed())
	}
	colnames(bestMatchesSorted)=c(name1, name2, "p.adj")
	return(bestMatchesSorted)
}

#' Overlap comparisons
#'
#' Compares modules between two objects of type WGCNAobjects
#' within a WGCNAobject list given the indices
#'
#' @param comparisonList a list passed by the iterate function
#' @param WGCNAlist list of objects of type WGCNAobject
#' @param first index of first WGCNAobject
#' @param second index of second WGCNAobject
#' @param element element position in the comparison list (passed by iterate function)
#' @param plot generate plots?
#'
#' @author Dario Tommasini
#'
#' @export
overlapComparisons <- function(comparisonList, WGCNAlist, first, second, element, plot=TRUE, write=FALSE){
	comparisonList[[element]]=list()
	comparisonList[[element]]=append(comparisonList[[element]],
					list(computeOverlapsFromWGCNA(WGCNAlist[[first]], WGCNAlist[[second]])))
	names(comparisonList[[element]])=append(names(comparisonList[[element]]), c("overlap"))
	names(comparisonList)[[element]]=paste0(names(WGCNAlist)[[first]], "_vs_", names(WGCNAlist)[[second]])
	if(write) write.csv(comparisonList[[element]]$overlap, paste0(names(WGCNAlist)[[first]], "_vs_", names(WGCNAlist)[[second]], ".csv"), row.names=F)
	cat("\n#### comparing ", names(WGCNAlist)[[first]], " and ", names(WGCNAlist)[[second]], "####\n")
	if(plot){
		moduleComparisonPlot(comparisonList[[element]]$overlap,
						WGCNAlist[[first]],
						WGCNAlist[[second]])
	}
	comparisonList[[element]]=append(comparisonList[[element]],
								list(bidirectionalBestMatches(comparisonList[[element]],
									plot=plot)))
	names(comparisonList[[element]])[[2]]=("bestMatches")
	comparisonList
}
