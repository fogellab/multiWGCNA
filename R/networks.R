#' Draw multiWGCNA network
#'
#' Draw a network where nodes are modules and edges represent significant gene overlap. 
#' Modules are sorted by levels 1, 2, and 3.  
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
#' @importFrom igraph graph_from_data_frame delete.edges V V<- E E<- plot.igraph
#' @import stringr 
#' @importFrom scales rescale
#' @export
#' 
#' @examples
#' library(ExperimentHub)
#' eh = ExperimentHub()
#' eh_query = query(eh, c("multiWGCNAdata"))
#' astrocyte_se = eh_query[["EH8223"]]
#' sampleTable = colData(astrocyte_se)
#' astrocyte_networks = eh_query[["EH8222"]]
#' results = list()
#' results$overlaps = iterate(astrocyte_networks, overlapComparisons, plot=FALSE)
#' drawMultiWGCNAnetwork(astrocyte_networks, 
#'   results$overlaps, 
#'   "combined_013", 
#'   sampleTable)
#'   
drawMultiWGCNAnetwork <- function(WGCNAlist, comparisonList, moduleOfInterest, design, 
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
			filteredOverlapList=filteredOverlapList[!filteredOverlapList$mod1 %in% WGCNA@outlierModules,]
			filteredOverlapList=filteredOverlapList[!filteredOverlapList$mod2 %in% WGCNA@outlierModules,]
		}
	}
	
	admittedModules=unique(c(filteredOverlapList$mod1, filteredOverlapList$mod2))
	print(admittedModules)
	
	#generate the multiWGCNA layout
	if(is.null(layout)){
		myCoords=list()
		for(level in 1:3){
			WGCNAs=getLevel(level, design)
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
	list(V(graph)$name, igraph::layout.fruchterman.reingold(graph), igraph::as_edgelist(graph))
}


