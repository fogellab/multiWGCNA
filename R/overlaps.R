#' TOMFlowPlot
#'
#' Plots a sankey flow diagram showing the movement of genes from one WGCNA
#' to another WGCNA. Uses the ggalluvial framework. 
#'
#' @param WGCNAlist list of WGCNA objects
#' @param networks list of network names of length 2
#' @param toms a list of TOM distance objects of length 2
#' @param genes_to_label genes to label across two networks
#' @param alpha alpha of flows
#' @param color color of flows
#' @param width width of the strata
#' 
#' @return a ggplot object
#' 
#' @author Dario Tommasini
#'
#' @import flashClust
#' @import ggalluvial
#' @import WGCNA
#' @import stringr
#' @export
TOMFlowPlot = function(WGCNAlist, networks, toms, genes_to_label, alpha = 0.1, color = 'black', width = 0.05){
  
  stopifnot(length(networks) == 2)
  stopifnot(length(toms) == 2)
  stopifnot(length(genes_to_label) > 0)
  
  # Read in files (equivalent to object@datExpr slot of WGCNA object)
  datasets = WGCNAlist[networks]
  
  # Make common list of genes
  allCommonGenes=datasets[[1]]@datExpr$X[order(datasets[[1]]@datExpr$X)]
  
  # First extract the order from hierarchical clustering of TOM dissimilarity matrix
  orderList=list()
  positionList=list()
  for(element in seq_along(networks)){
    
    # Get TOM
    TOM=toms[[element]]
    dissTOM=1-TOM
    
    # Cluster TOM
    message("Clustering TOM tree ", element, "...")
    geneTree = flashClust(as.dist(dissTOM), method="average")
    
    # Compile data
    geneOrder=data.frame(Gene=datasets[[element]]@datExpr$X[geneTree$order], 
                         Color=datasets[[element]]@datExpr$dynamicColors[geneTree$order])
    geneOrder$order=1:length(allCommonGenes)
    sortedOrder=geneOrder[order(geneOrder$Gene),]
    orderList[[element]]=sortedOrder
  }
  
  # Combine results
  df=do.call(cbind, orderList)
  colnames(df)=c("one.gene","one.color", "one.order", "two.gene", "two.color", "two.order")
  df$Count=1
  head(df)
  df$one.gene=paste0("1_", df$one.gene)
  df$one.gene=factor(df$one.gene, levels=df$one.gene[order(df$one.order)])
  df$two.gene=paste0("2_", df$two.gene)
  df$two.gene=factor(df$two.gene, levels=df$two.gene[order(df$two.order)])
  colors=c(rev(df$one.color[order(df$one.order)]), rev(df$two.color[order(df$two.order)]))
  
  # Choose module to label
  df$module="other"
  df$module[allCommonGenes %in% genes_to_label]="labeled"
  df$module = factor(df$module, levels = c("other", "labeled"))
  
  # Set unlabeled genes to NAs to avoid labeling them, speeds up plotting as well
  temp1 = subset(df, module != "labeled")
  temp1$two.gene = NA
  temp2 = subset(df, module != "labeled")
  temp2$one.gene = NA
  
  new_df = rbind(subset(df, module == "labeled"), 
                 temp1, temp2)
  
  # Flow plot
  plt = ggplot(new_df, aes(y = Count, axis1 = one.gene, axis2 = two.gene)) +
      geom_flow(aes(fill = module), width=width, curve_type = "cubic", alpha = alpha, fill = color) +
      geom_stratum(width = width, fill = colors, size=0, alpha = 1) +
      ylab("Genes") + 
      theme(axis.ticks.x = element_blank(), panel.background = element_blank(),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
      geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 0) +
      scale_y_continuous(limits = c(0, nrow(df)))+
      scale_x_discrete(expand=c(0,0), limits = networks, labels = networks)
  
  return(plt)
}

#' BuildTOMFlowDF
#'
#' Preprocess for plotting a sankey flow diagram showing the movement of genes from one WGCNA
#' to another WGCNA. Uses the flashClust framework. 
#'
#' @param WGCNAlist list of WGCNA objects
#' @param networks list of network names of length 2
#' @param toms a list of TOM distance objects of length 2
#' @param genes_to_label genes to label across two networks
#' 
#' @return a data.frame
#' 
#' @author Dario Tommasini, Xinye Li
#'
#' @import flashClust
#' @import WGCNA
#' @import stringr
#' @export
BuildTOMFlowDF <- function(WGCNAlist, networks, toms, genes_to_label) {
  
  stopifnot(length(networks) == length(toms))
  stopifnot(length(genes_to_label) > 0)
  
  datasets <- WGCNAlist[networks]
  allCommonGenes <- datasets[[1]]@datExpr$X[order(datasets[[1]]@datExpr$X)]
  
  # Calculate the ranking information for each network
  orderList <- list()
  colorList <- list()
  
  for(i in seq_along(networks)){
    TOM <- toms[[i]]
    dissTOM <- 1 - TOM
    
    message("Clustering TOM tree ", networks[i], "...")
    geneTree <- flashClust::flashClust(as.dist(dissTOM), method="average")
    
    geneOrder <- data.frame(
      Gene = datasets[[i]]@datExpr$X[geneTree$order],
      Color = datasets[[i]]@datExpr$dynamicColors[geneTree$order],
      Order = seq_along(allCommonGenes)
    )
    sortedOrder <- geneOrder[order(geneOrder$Gene), ]
    orderList[[i]] <- sortedOrder
    colorList[[i]] <- sortedOrder$Color
  }
  
  # Merge all orderList into one dataframe
  combined_df <- data.frame(Gene = orderList[[1]]$Gene)
  
  for(i in seq_along(orderList)){
    combined_df[[paste0("G", i)]] <- paste0(i, "_", orderList[[i]]$Gene)
    combined_df[[paste0("Color", i)]] <- colorList[[i]]
    combined_df[[paste0("Order", i)]] <- orderList[[i]]$Order
  }
  
  combined_df$Count <- 1
  
  # Construct the information of the module tag
  combined_df$module <- "other"
  combined_df$module[combined_df$Gene %in% genes_to_label] <- "labeled"
  combined_df$module <- factor(combined_df$module, levels = c("other", "labeled"))
  
  # Construct new_df (long format)
  long_df_list <- list()
  for(i in seq_along(networks)){
    temp_df <- combined_df[, c("Gene", paste0("G", i), paste0("Order", i), paste0("Color", i))]
    colnames(temp_df) <- c("Gene", "Label", "Order", "Color")
    temp_df$Network <- paste0("Network", i)
    temp_df$module <- combined_df$module
    long_df_list[[i]] <- temp_df
  }
  
  new_df <- do.call(rbind, long_df_list)
  
  # Change to wide format
  wide_df <- tidyr::pivot_wider(
    new_df,
    id_cols = Gene,
    names_from = Network,
    values_from = c(Label, Order, Color)
  )
  
  # Restore module information
  wide_df$module <- combined_df$module[match(wide_df$Gene, combined_df$Gene)]
  wide_df$Count <- 1
  wide_df <- as.data.frame(wide_df)
  
  return(wide_df)
}

#' PlotMultiNodesTOMflow
#'
#' Plots a sankey flow diagram showing the movement of genes from one WGCNA
#' to multi-WGCNA networks. Uses the ggalluvial framework. 
#'
#' @param TOMDF created by BuildTOMFlowDF
#' @param alpha alpha of flows
#' @param width width of the strata
#' @param color color of flows
#' 
#' @return a ggplot object
#' 
#' @author Dario Tommasini, Xinye Li
#'
#' @import flashClust
#' @import ggalluvial
#' @import WGCNA
#' @import stringr
#' @export
PlotMultiNodesTOMflow <- function(TOMDF, alpha = 0.1, width = 0.05, color = "black") {
  
  # Split data
  networks <- grep("^Label_Network", colnames(TOMDF), value = TRUE)
  orders <- grep("^Order_Network", colnames(TOMDF), value = TRUE)
  colors <- grep("^Color_Network", colnames(TOMDF), value = TRUE)
  
  # Equal the number of networks
  network_count <- length(networks)
  if(network_count < 2) {
    stop("At least 2 networks are needed for drawing a Sankey plot.")
  }
  if(network_count > 5) {
    stop("Currently only 2-5 networks are supported.")
  }
  
  # Build basic df
  df_columns <- c(networks, orders, colors, "module", "Count")
  df <- TOMDF[, df_columns]
  
  # Sortting
  for (i in seq_along(networks)) {
    df[[networks[i]]] <- factor(df[[networks[i]]], 
                                levels = df[[networks[i]]][order(df[[orders[i]]])])
  }
  
  # Adjust connections
  new_dfs <- list()
  
  # Keep connections between labeled genes
  new_dfs[[1]] <- subset(df, module == "labeled")
  
  # Rm connections between others
  for (i in seq_along(networks)) {
    temp <- subset(df, module != "labeled")
    # Set all other nodes' genes as NA, only keep current node for only drawing the current node
    for (j in seq_along(networks)) {
      if (j != i) {
        temp[[networks[j]]] <- NA
      }
    }
    new_dfs[[length(new_dfs) + 1]] <- temp
  }
  
  # Combine the df
  new_df <- do.call(rbind, new_dfs)
  
  # Set colors
  all_node_colors <- list()
  for (i in seq_along(networks)) {
    sorted_df <- df[order(df[[orders[i]]]), ]
    all_node_colors[[i]] <- rev(sorted_df[[colors[i]]])
  }
  
  # combine colors
  node_colors <- unlist(all_node_colors)
  
  # Build plot
  plt <- ggplot(new_df)
  
  if (network_count == 2) {
    plt <- plt + aes(y = Count, axis1 = !!sym(networks[1]), axis2 = !!sym(networks[2]))
  } else if (network_count == 3) {
    plt <- plt + aes(y = Count, axis1 = !!sym(networks[1]), axis2 = !!sym(networks[2]), 
                     axis3 = !!sym(networks[3]))
  } else if (network_count == 4) {
    plt <- plt + aes(y = Count, axis1 = !!sym(networks[1]), axis2 = !!sym(networks[2]), 
                     axis3 = !!sym(networks[3]), axis4 = !!sym(networks[4]))
  } else if (network_count >= 5) {
    plt <- plt + aes(y = Count, axis1 = !!sym(networks[1]), axis2 = !!sym(networks[2]), 
                     axis3 = !!sym(networks[3]), axis4 = !!sym(networks[4]), 
                     axis5 = !!sym(networks[5]))
  }
  
  # Nodes & flows
  plt <- plt +
    geom_flow(aes(fill = module), width = width, curve_type = "cubic", alpha = alpha, fill = color) +
    geom_stratum(width = width, fill = node_colors, size = 0, alpha = 1) +
    ylab("Genes") +
    theme(
      axis.ticks.x = element_blank(),
      panel.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 0) +
    scale_y_continuous(limits = c(0, nrow(df))) +
    scale_x_discrete(
      expand = c(0, 0), 
      limits = paste0("Network ", 1:network_count),
      labels = paste0("Network ", 1:network_count)
    )
  
  return(plt)
}

#' computeOverlapsFromWGCNA
#'
#' Computes overlap between the modules of two objects of class WGCNA
#'
#' @param dataset1 an object of class WGCNA to compare with dataset2
#' @param dataset2 an object of class WGCNA to compare with dataset1
#' 
#' @return Returns a data.frame showing the overlap results for modules from 
#' dataset1 with dataset2
#' 
#' @author Dario Tommasini
#'
#' @import stringr
#' @export
#' 
#' @examples
#' library(ExperimentHub)
#' eh = ExperimentHub()
#' eh_query = query(eh, c("multiWGCNAdata"))
#' astrocyte_networks = eh_query[["EH8222"]]
#' computeOverlapsFromWGCNA(astrocyte_networks$EAE, astrocyte_networks$WT)
#' 
computeOverlapsFromWGCNA <- function(dataset1, dataset2) {
  
  # Check input
  stopifnot(inherits(dataset1, "WGCNA") & inherits(dataset2, "WGCNA"))
  
	datExpr1 <- dataset1@datExpr
	datExpr2 <- dataset2@datExpr
	treatDat <- datExpr1
	controlDat <- datExpr2
	
	# computation of overlaps
	all.genes = sort(unique(c(treatDat$X, controlDat$X)))
	sorted.modules <- data.frame(mod1 = treatDat$dynamicLabels[match(all.genes, treatDat$X)],
	                             mod2 = controlDat$dynamicLabels[match(all.genes, controlDat$X)])

	overlap.count <- table(sorted.modules$mod1, sorted.modules$mod2)
  
	pval = overlap.count
	for(column in seq_along(colnames(overlap.count))){
	  for(row in seq_along(rownames(overlap.count))){
	    pval[row,column]=phyper(overlap.count[row,column]-1,
	                            sum(overlap.count[row,]),
	                            sum(rowSums(overlap.count))-sum(overlap.count[row,]),
	                            sum(overlap.count[,column]),
	                            lower.tail=FALSE,
	                            log.p=FALSE)
	    if(pval[row,column]==0) pval[row,column]=.Machine$double.xmin
	  }
	}

	# Adjust for multiple comparisons
	pval.unlist = unlist(as.list(pval))
	pval.adj = p.adjust(pval.unlist, method='fdr')

	# Get module sizes
	mod1.size = table(sorted.modules$mod1)
	mod2.size = table(sorted.modules$mod2)

	# Summarize results
	output.df = reshape2::melt(overlap.count) 
	colnames(output.df) = c("mod1", "mod2", "overlap")
	output.df$mod1 = as.character(output.df$mod1)
	output.df$mod2 = as.character(output.df$mod2)
	output.df$p.value = pval.unlist
	output.df$p.adj = pval.adj
	output.df = output.df %>% arrange(mod1, mod2)
	output.df$mod1.size = mod1.size[match(output.df$mod1, names(mod1.size))]
	output.df$mod2.size = mod2.size[match(output.df$mod2, names(mod2.size))]

	# Return columns in proper order
	return(output.df[,c(1:2,6:7,3:5)] %>% arrange(mod1, mod2))
}

#' Module comparison plot
#'
#' A plotting function that returns a heatmap and barplot for a module
#'
#' @param overlapDf a data.frame resulting from a call to computeOverlapsFromWGCNA
#' @param dataset1 an object of class WGCNA to compare with dataset2
#' @param dataset2 an object of class WGCNA to compare with dataset1
#'
#' @return Returns a ggplot object (flowplot and heatmap) showing the module correspondence 
#' between two objects of class WGCNA
#'
#' @author Dario Tommasini
#'
#' @import ggplot2
#' @import ggalluvial
#' @import stringr
#' @importFrom cowplot plot_grid
#' 
#' @export
#' 
#' @examples
#' library(ExperimentHub)
#' eh = ExperimentHub()
#' eh_query = query(eh, c("multiWGCNAdata"))
#' astrocyte_networks = eh_query[["EH8222"]]
#' overlapDf = computeOverlapsFromWGCNA(astrocyte_networks$EAE, astrocyte_networks$WT)
#' moduleComparisonPlot(overlapDf, astrocyte_networks$EAE, astrocyte_networks$WT)
#' 
moduleComparisonPlot <- function(overlapDf, dataset1, dataset2) {
  
  # Check input
  stopifnot(inherits(dataset1, "WGCNA") & inherits(dataset2, "WGCNA"))
  
	data=overlapDf
	name1=str_split_fixed(data$mod1,"_",2)[,1][[1]]
	name2=str_split_fixed(data$mod2,"_",2)[,1][[1]]
	data$mod1=str_split_fixed(data$mod1,"_",2)[,2]
	data$mod2=paste0(str_split_fixed(data$mod2,"_",2)[,2]," ") # add spaces to avoid ambiguous plot labels
	categories=unique(c(dataset1@trait$trait, dataset2@trait$trait))
	palette=colors(length(categories), random=FALSE)
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
		annotate("text", x=2.2, y=(1:length(categories))*(totalGenes/3)/length(categories),
			label=categories, vjust=0, hjust=0, color=palette)+
		geom_text(stat = "stratum", aes(label = gsub("^0+", "", after_stat(stratum))), size = 3, min.y=50) +
		coord_cartesian(xlim = c(0.9, 2.5), clip = "off") +
 		scale_x_discrete(expand=c(0,0), limits = c("Mod1", "Mod2"), labels=c(name1,name2))

	plot = plot_grid(flowPlot, NULL, heatmap, rel_widths = c(1, 0.2, 2.5), nrow = 1)
	
	return(plot)
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
#' @param comparisonDf the data.frame output of computeOverlapFromWGCNA
#' @param dataset1 optional; WGCNA object for dataset 1
#' @param dataset2 optional; WGCNA object for dataset 2
#' @param trait1 optional; subset to modules correlated to this trait for dataset 1
#' @param trait2 optional; subset to modules correlated to this trait for dataset 2
#' @param list1 subset to this list of modules for dataset 1
#' @param list2 subset to this list of modules for dataset 2
#' @param filterByTrait only plot for modules that correlate with some trait?
#' @param alphaLevel the alpha level of significance for module-trait correlation, defaults to 0.05
#'
#' @return A ggplot object
#'
#' @import stringr
#' @import ggplot2
#' @export
#' 
#' @examples
#' library(ExperimentHub)
#' eh = ExperimentHub()
#' eh_query = query(eh, c("multiWGCNAdata"))
#' astrocyte_networks = eh_query[["EH8222"]]
#' overlapDf = computeOverlapsFromWGCNA(astrocyte_networks$EAE, astrocyte_networks$WT)
#' moduleToModuleHeatmap(overlapDf)
#' 
moduleToModuleHeatmap <- function(comparisonDf, dataset1=NULL, dataset2=NULL, trait1=NULL, trait2=NULL, list1=NULL, list2=NULL, filterByTrait=FALSE, alphaLevel=0.05){
  
  # Check input
  stopifnot(inherits(comparisonDf, "data.frame"))
  
	# filter by trait if desired
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

	plot = ggplot(comparisonDf, aes(x = factor(mod1, levels=(unique(mod1))),
	                         y = factor(mod2, levels=rev(unique(mod2))),
 				fill = (-log10(p.adj)), label = overlap)) +
				geom_tile(color = "black") +
				scale_fill_gradient(name="-log10(FDR)", low = "white",
					high = "red", na.value="red", #limits=c(0, 50),
					guide = guide_colorbar(frame.colour = "black", frame.linewidth = 0.5,
							ticks.linewidth = 0.5, ticks.colour = "black")) +
				geom_text(aes(label = overlap), color = "black", size=2) +
				labs(x=columns, y=rows) +
				theme_classic()+
	      ggtitle("All matches") +
				theme(axis.text.x = element_text(angle = 0, hjust=(0.5)), #angle = 90, hjust=1, vjust=(0.5)),
					axis.text.y = element_text(hjust=1), plot.title = element_text(hjust = 0.5),
					panel.background=element_blank(), axis.line=element_blank(),
					axis.ticks=element_blank(), legend.key.size=unit(4, "mm")) +
				coord_fixed()
	
	return(plot)
}

#' Best matching modules
#'
#' Find all the modules from dataset1 that have a best match to a module in dataset2
#' if that module in dataset2 is also a best match to the module in dataset1
#'
#' @param comparisonList a list with an elemnt "overlap", which is a data.frame resulting from a call to computeOverlapsFromWGCNA
#' @param plot whether to generate a heatmap; default is TRUE
#'
#' @return A ggplot object
#' 
#' @author Dario Tommasini
#'
#' @import ggplot2
#' @import stringr
#' @import dplyr
#' @export
#' 
#' @examples
#' library(ExperimentHub)
#' eh = ExperimentHub()
#' eh_query = query(eh, c("multiWGCNAdata"))
#' astrocyte_networks = eh_query[["EH8222"]]
#' comparisonList = list()
#' comparisonList$overlaps = computeOverlapsFromWGCNA(astrocyte_networks$EAE, astrocyte_networks$WT)
#' bidirectionalBestMatches(comparisonList)
#' 
bidirectionalBestMatches <- function(comparisonList, plot=TRUE){
  
  # Check input
  stopifnot(inherits(comparisonList, "list"))
  
  comparison=comparisonList
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
	if(plot){
		plt = ggplot(subsetComparisonDf, 
		             aes(x = factor(mod1, levels=bestMatchesSorted$mod1),
		                y = factor(mod2, levels=rev(bestMatchesSorted$mod2)),
 				fill = (-log10(p.adj)))) +
				geom_tile(color = "black") +
			  scale_x_discrete(labels = function(x) gsub("^0+", "", x))+
				scale_y_discrete(labels = function(x) gsub("^0+", "", x))+
				scale_fill_gradient(name="-log10(FDR)", low = "white",
				                      high = "red", na.value="red", #limits=c(0, 50),
				                      guide = guide_colorbar(frame.colour = "black", frame.linewidth = 0.5,
				                                             ticks.linewidth = 0.5, ticks.colour = "black")) +
				geom_text(aes(label = overlap), color = "black") +
				labs(x=name1, y=name2) +
			  ggtitle("Best matches") +
  		  theme(axis.text.x = element_text(angle = 0, hjust=(0.5)), #angle = 90, hjust=1, vjust=(0.5)),
  		        axis.text.y = element_text(hjust=1), plot.title = element_text(hjust = 0.5),
  		        panel.background=element_blank(), axis.line=element_blank(),
  		        axis.ticks=element_blank(), legend.key.size=unit(4, "mm"))
				coord_fixed()
		print(plt)
	}
	colnames(bestMatchesSorted)=c(name1, name2, "p.adj")
	return(bestMatchesSorted)
}

#' Overlap comparisons
#'
#' Compares modules between two objects of type WGCNAobjects
#' within a WGCNAobject list given the indices. Recommended to be used in 
#' conjunction with the iterate function. 
#'
#' @param comparisonList a list passed by the iterate function
#' @param WGCNAlist list of objects of class WGCNA
#' @param first index of first WGCNA object
#' @param second index of second WGCNA object
#' @param element element position in the comparison list (passed by iterate function)
#' @param plot generate plots?
#' @param write write results to file?
#' 
#' @return A list, in which the first element is a data.frame showing the 
#' overlap results and the second element is a data.frame showing the best 
#' matching modules between the two WGCNA objects. 
#'
#' @author Dario Tommasini
#'
#' @export
#' 
#' @examples
#' library(ExperimentHub)
#' eh = ExperimentHub()
#' eh_query = query(eh, c("multiWGCNAdata"))
#' astrocyte_networks = eh_query[["EH8222"]]
#' results = list()
#' results$overlaps = iterate(astrocyte_networks, overlapComparisons, plot=FALSE)
#' 
overlapComparisons <- function(comparisonList, WGCNAlist, first, second, element, plot=TRUE, write=FALSE){
  
  # Check input
  stopifnot(inherits(comparisonList, "list"))
  stopifnot(inherits(WGCNAlist[[1]], "WGCNA"))
  
	comparisonList[[element]]=list()
	comparisonList[[element]]=append(comparisonList[[element]],
					list(computeOverlapsFromWGCNA(WGCNAlist[[first]], WGCNAlist[[second]])))
	names(comparisonList[[element]])=append(names(comparisonList[[element]]), c("overlap"))
	names(comparisonList)[[element]]=paste0(names(WGCNAlist)[[first]], "_vs_", names(WGCNAlist)[[second]])
	if(write) write.csv(comparisonList[[element]]$overlap, paste0(names(WGCNAlist)[[first]], "_vs_", names(WGCNAlist)[[second]], ".csv"), row.names=FALSE)
	message("\n#### comparing ", names(WGCNAlist)[[first]], " and ", names(WGCNAlist)[[second]], " ####\n")
	if(plot){
		print(
		  moduleComparisonPlot(comparisonList[[element]]$overlap,
						WGCNAlist[[first]],
						WGCNAlist[[second]])
		  )
	}
	comparisonList[[element]]=append(comparisonList[[element]],
								list(bidirectionalBestMatches(comparisonList[[element]],
									plot=plot)))
	names(comparisonList[[element]])[[2]]=("bestMatches")
	
	return(comparisonList)
}
