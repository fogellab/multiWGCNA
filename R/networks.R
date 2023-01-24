
#plot(g3, layout=t(cbind(c(2,2),c(1,1),c(3,1),c(0,3),c(1.33,3),c(2.66,3),c(4,3))), vertex.shape="rectangle", edge.width=rescale(c(100, 20, 30, 40, 50, 60), to=c(1,10)), vertex.size=(strwidth(V(g3)$name)+strwidth("oo")) * 100, vertex.size2=15, vertex.color="white", vertex.label.family="Arial", vertex.label.color="black")

flowNetwork <- function(comparisonList, continuousTraits, overlapCutoff=NULL, padjCutoff=NULL){
	
	series=lapply(1:3, function(element) paste0(continuousTraits[[element]],"_vs_", continuousTraits[[element+1]]))
	cleanComparisonList=comparisonList[names(comparisonList) %in% series]
	overlapList=lapply(cleanComparisonList, function(x) x$overlap)
	overlapList=do.call(rbind, overlapList)
	if(missing(overlapCutoff) | missing(padjCutoff)){
		filteredOverlapList=overlapList
	} else {
		filteredOverlapList=overlapList[overlapList$overlap>overlapCutoff & overlapList$p.adj<padjCutoff, ]
	}
	graph=graph_from_data_frame(d=filteredOverlapList, directed = FALSE)
	
	#sankey-like layout
	allModules=V(graph)$name
	moduleSizes=lapply(1:length(continuousTraits), function(element) length(allModules[startsWith(allModules, continuousTraits[[element]]) ])) 
	maxSize=max(unlist(moduleSizes))
	maxSize=1000
	layout_as_flow=rbind(cbind(0, seq(maxSize, 0, length.out=moduleSizes[[1]])), cbind(1, seq(maxSize, 0, length.out=moduleSizes[[2]])), 
						cbind(2, seq(maxSize, 0, length.out=moduleSizes[[3]])), cbind(3, seq(maxSize, 0, length.out=moduleSizes[[4]])))
	
	#node and edge attributes
	vcol=str_split_fixed(V(graph)$name, "_", 2)[,1]
	conditions=unique(vcol)
	V(graph)$name=str_split_fixed(V(graph)$name, "_", 2)[,2]
	palette=brewer.pal(n = length(conditions), name = "Dark2")
	for(condition in 1:length(conditions)){
		vcol[vcol==conditions[[condition]] ]=palette[[condition]]
	}
	V(graph)$color=vcol
	
	#edge attributes
	E(graph)$width=rescale(E(graph)$overlap, to=c(0,5))
	ealpha=rescale(-log10(E(graph)$p.adj), to=c(0,1))
	#ecol= rep("black", length(E(graph)))
	if(!missing(moduleOfInterest)) ecol[filteredOverlapList$mod1==moduleOfInterest | filteredOverlapList$mod2==moduleOfInterest]="red"
	ecol=lapply(ealpha, function(x) rgb(0, 0, 0, x))
	E(graph)$color=rgb(0,0,0,0.5)
	
	#plot graph
	plot(graph, vertex.label.color="black", vertex.label.size=3, vertex.size=5, vertex.label.cex=0.6, 
		vertex.shape="crectangle", vertex.size2=2, edge.color=unlist(ecol), layout=layout_as_flow)
	legend("right", legend = conditions, pch=21,
       col=palette, pt.bg=palette, pt.cex=1, cex=.8, bty="n", ncol=1)
	
	#second call to plot over with edges of interest
	E(graph)[E(graph)$color=="grey"]$color  <- NA
	plot(graph, vertex.label.color="black", vertex.size=5,
		vertex.label.cex=.5, layout=layout_in_circle, add=TRUE)
}

#' constructNetworks: Construct all the weighted gene correlation networks
#'
#' A high level function that returns all networks
#' possible for a given experimental design
#'
#' @param WGCNAlist list of WGCNA objects
#' @param comparisonList the list of overlap comparisons ie from iterate(myNetworks, overlapComparisons, ...) 
#' @param moduleOfInterest module of interest, ie "combined_001"
#' @param design the sampleTable design matrix
#' @param overlapCutoff cutoff to remove module correspondences with less than this number of genes
#' @param padjCutoff cutoff to remove module correspondences above this significance value
#' @param removeOutliers remove outlier modules? 
#' @param alpha alpha level of significance
#' @param layout layout of network to be passed to plot function of igraph object, defaults to multiWGCNA custom layout
#' @param hjust horizontal justification of labels
#' @param vjust vertical justification of labels
#' @param width width of labels
#' @param colors colors to use for modules, should be the same length as the number of WGCNA objects in the WGCNAlist. Defaults to random colors for each condition. 
#'
#' @author Dario Tommasini
#'
#' @import igraph
#' @import stringr 
#' @import scales
#' @export
drawMultiWGCNAnetwork <- function(WGCNAlist, comparisonList, moduleOfInterest, design = sampleTable, 
                                  overlapCutoff = 0, padjCutoff = 1, removeOutliers = T, alpha = 1e-50, 
                                  layout = NULL, hjust = 0.4, vjust = 0.3, width = 0.5, colors = NULL){

  	#extract the overlaps objects into a list
	overlapList=lapply(comparisonList, function(x) x$overlap)
	overlapList=do.call(rbind, overlapList)
	#filteredOverlapList=overlapList[overlapList$overlap>overlapCutoff & overlapList$p.adj<padjCutoff, ]
	filteredOverlapList=overlapList

	#remove outliers if necessary
	if(removeOutliers) {
		for(WGCNA in WGCNAlist){
			filteredOverlapList= filteredOverlapList[!filteredOverlapList$mod1 %in% WGCNA@outlierModules,]
			filteredOverlapList= filteredOverlapList[!filteredOverlapList$mod2 %in% WGCNA@outlierModules,]
		}
	}
	
	admittedModules=unique(c(filteredOverlapList$mod1, filteredOverlapList$mod2))
	#print(admittedModules)
	
	#generate the multiWGCNA layout
	if(is.null(layout)){
		myCoords=list()
		for(level in 1:3){
			WGCNAs=getLevel(level)
			from=0-width*length(WGCNAs)/2
			to=0+width*length(WGCNAs)/2
			x.coordinates=seq(from, to, length.out=length(WGCNAs))
			if(level==1) x.coordinates=0
			for(nWGCNA in 1:length(WGCNAs)){
				nModules=length(admittedModules[startsWith(admittedModules, WGCNAs[[nWGCNA]]) ])
				myCoords=append(myCoords, list(cbind(runif(nModules, x.coordinates[[nWGCNA]], x.coordinates[[nWGCNA]]+hjust), 3-level+runif(nModules, -vjust, vjust))))
			}
		}
		layout=do.call(rbind, myCoords)
	}
	
	#make the igraph object
	graph=graph_from_data_frame(d=filteredOverlapList, directed = FALSE)
	
	#node and edge attributes
	vcol=str_split_fixed(V(graph)$name, "_", 2)[,1]
	conditions=unique(vcol)
	#V(graph)$name=str_split_fixed(V(graph)$name, "_", 2)[,2]
	
	# Colors of modules by condition
	if(is.null(colors)) palette = colors(length(conditions), random = T)
	if(!is.null(colors)) palette = colors
	# if(random.colors) palette = colors(length(conditions), random = F)
	
	for(condition in 1:length(conditions)){
		vcol[vcol==conditions[[condition]] ]=palette[[condition]]
	}
	V(graph)$color=vcol
	
	#edge attributes
	E(graph)$weight=-log10(E(graph)$p.adj)
	#E(graph)$width=rescale(E(graph)$overlap, to=c(0, 5))
	E(graph)$width=rescale(E(graph)$weight, from=c(0, 320), to=c(0,5))
	#ecol= rep("grey", length(E(graph)))
	#ecol[filteredOverlapList$mod1==moduleOfInterest | filteredOverlapList$mod2==moduleOfInterest]="red"
	ealpha=rescale(-log10(E(graph)$p.adj), from=c(0, 320), to=c(0,1))
	ecol=lapply(ealpha, function(x) rgb(1, 0, 0, x))
	E(graph)$color=unlist(ecol)
	
	#plot graph
	#plot(graph, vertex.label.color="black", vertex.size=5,
	#	vertex.label.cex=.5, layout=layout_in_circle)
	#legend("right", legend = conditions, pch=21,
    #   col=palette, pt.bg=palette, pt.cex=1, cex=.8, bty="n", ncol=1)
	
	#second call to plot over with edges of interest
	#E(graph)[!(filteredOverlapList$mod1==moduleOfInterest | filteredOverlapList$mod2==moduleOfInterest)]$color  <- NA
	modulesOfInterest=unique(c(filteredOverlapList$mod2[filteredOverlapList$mod1==moduleOfInterest & filteredOverlapList$p.adj<alpha], 
		filteredOverlapList$mod1[filteredOverlapList$mod2==moduleOfInterest & filteredOverlapList$p.adj<alpha]))
	# print(modulesOfInterest)
	#E(graph)[!(filteredOverlapList$mod1 %in% modulesOfInterest | filteredOverlapList$mod2 %in% modulesOfInterest)]$color  <- NA
	
	#delete edges
	#graph <- delete.edges(graph, which(E(graph)$p.adj>padjCutoff | E(graph)$overlap<overlapCutoff))
	graph <- delete.edges(graph, which(!(filteredOverlapList$mod1 %in% modulesOfInterest | filteredOverlapList$mod2 %in% modulesOfInterest)))
	# print(which(!(E(graph)$mod1 %in% modulesOfInterest | E(graph)$mod2 %in% modulesOfInterest)))
	#neighborhood <- neighborhood(graph, 1, nodes=modulesOfInterest, mindist=1)
	#res.tokeep <- lapply(res, function(x) which(E(graph)[x]$weight>100))
	#res.todelete <- lapply(res, function(x) which(E(graph)[x]$weight<=100))
	#ntwrk <- delete.vertices(graph, unique(unlist(res.todelete)))

	plot = plot(graph, vertex.label.color="black", vertex.size=3, vertex.label=NA, 
	            vertex.label.cex=.5, layout=layout) 
	legend("right", legend = conditions, pch=21, col=palette, pt.bg=palette, 
	       pt.cex=1, cex=.8, bty="n", ncol=1)
  return(plot)

	#filteredOverlapList[filteredOverlapList$mod1 %in% modulesOfInterest | filteredOverlapList$mod2 %in% modulesOfInterest,]
	
	#ego(graph, order=1, nodes="combined_001", mindist=300)
}
	# old_graph=graph
	# new_graph=graph_from_data_frame(d=as_data_frame(old_graph, what="edges") %>%
    	# arrange(desc(ordering_attribute)),
    	# vertices=as_data_frame(old_graph, what="vertices"))

	#V(graph)$size will equal module size
	#E(graph)$color will equal -log10 p.adj value
	
	#plot(graph, edge.arrow.size=.5, vertex.color="gold", vertex.size=3, 
    #vertex.frame.color="gray", vertex.label.color="black", 
    #vertex.label.cex=.5, vertex.label.dist=2, edge.curved=0.5,layout=layout_with_lgl)


estimateTOM <- function(datExpr, geneList, softPower=12){
	filteredDatExpr=datExpr[datExpr$X %in% geneList,]
	cleanDatExpr=cleanDatExpr(filteredDatExpr)
	adjacency = adjacency(cleanDatExpr, power=softPower, type="signed");
	TOM = TOMsimilarity(adjacency);
	dimnames(TOM)=dimnames(adjacency)
	TOM[TOM==1]=0
	TOM
}

#draw a basic network
drawNetwork <- function(matrix, threshold=0, nodeList=NULL, edgeList=NULL, layout=layout_with_fr, removeFreeNodes=TRUE){
	graph=graph_from_adjacency_matrix(matrix, mode="undirected", weighted=TRUE)
	graph <- delete.edges(graph, which(E(graph)$weight < threshold))
	if(removeFreeNodes) graph <- delete.vertices(simplify(graph), degree(graph)==0)
	ealpha=rescale(E(graph)$weight, from=c(0,1), to=c(0,0.4))
	ecol=lapply(ealpha, function(x) rgb(0.2, 0.2, 0.2, x))
	#ecol=lapply(E(graph)$weight, function(x) rgb(1-x, 1-x, 1-x))
	#ewidth=rescale(E(graph)$weight, from=c(0,1), to=c(0, 1.5))
	
	if(length(nodeList)>0){
		graph <- delete.vertices(graph, which(!V(graph)$name %in% nodeList))
		layout=layout[match(V(graph)$name, nodeList),]
	}
	
	if(length(edgeList)>0){
		graph <- delete.edges(graph, which(!apply(as_edgelist(graph), 1, function(x) 
			any(unlist(lapply(1:nrow(edgeList), function(y) {
				x[[1]]==edgeList[y,1] & x[[2]]==edgeList[y,2]
			}))))))
	}
		
	plot(graph, layout=layout, vertex.size=4, edge.width=1.5, 
		vertex.color="white", vertex.label=NA, edge.color=unlist(ecol))
	list(V(graph)$name, layout.fruchterman.reingold(graph), as_edgelist(graph))
}


