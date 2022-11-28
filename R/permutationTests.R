#' Disease preservation permutation test
#'
#' Tests whether the disease
#'
#' @param referenceDatExpr data.frame containing expression levels
#' @param design a data.frame with the columns corresponding to sample traits of referenceDatExpr
#' @param nPermutations number of permutations to perform
#' @param refColumn column of design table with reference traits
#' @param testColumn column of design table with test traits
#' @param nCores number of cores to use, default is 1
#' @param saveRData save objects to an RData file?
#' 
#' @author Dario Tommasini
#'
#' @import WGCNA
#' @import stringr
#' @export
diseasePreservationPtest <- function(referenceDatExpr, design, nPermutations=20, refColumn=3, testColumn=2, nCores=1, saveRData=T){

  registerDoParallel(cores=nCores)
  enableWGCNAThreads(nThreads=nCores)

  datExpr=referenceDatExpr[,match(design$Sample, colnames(referenceDatExpr))]

  wtDatExpr=list()
  dsDatExpr=list()
  WGCNAobjects=list()
  preservationData=list()
  filteredPreservationData=list()

  for(permutation in 1:nPermutations){
    
    #assign phenotype labels randomly
    WT.indices=list()
    conditions=unique(design[, refColumn])
    for(condition in conditions){
      conditionalDesign=design[design[, refColumn]==condition,]
      nSamples=nrow(conditionalDesign)
      nWT=nrow(conditionalDesign[conditionalDesign[, testColumn]=="WT",])
      sampleIndices=which(design[, refColumn]==condition)
      WT.indices=append(WT.indices, sampleIndices[sample(1:nSamples, nWT, replace=F)])
    }
    WT.indices=unlist(WT.indices)
    randomHealthy=datExpr[, WT.indices]
    randomHealthy=data.frame(X=referenceDatExpr$X, randomHealthy)
    randomDisease=datExpr[, !(1:ncol(datExpr) %in% WT.indices)]
    randomDisease=data.frame(X=referenceDatExpr$X, randomDisease)

    #perform WGCNA
    dir.create("WGCNAs")
    setwd("WGCNAs")
    wtDatExpr[[permutation]] <- runWGCNA_minimal(randomHealthy, paste0("WtP", permutation), softPowerCalculation=F)
    dsDatExpr[[permutation]] <- runWGCNA_minimal(randomDisease, paste0("DsP", permutation), softPowerCalculation=F)
    setwd("..")
    gc()

    #perform preservation
    dir.create("preservation")
    setwd("preservation")
    preservationData[[permutation]]=list()
    preservationData[[permutation]][[1]] <- getPreservation(wtDatExpr[[permutation]], dsDatExpr[[permutation]], write=T)
    preservationData[[permutation]][[2]] <- getPreservation(dsDatExpr[[permutation]], wtDatExpr[[permutation]], write=T)
    setwd("..")
    save.image("presPermutation.RData")
    WGCNAobjects[[permutation]]=findOutlierModules(findModuleEigengenes(dsDatExpr[[permutation]]))
    filteredPreservationData[[permutation]]=preservationData[[permutation]][[2]][!rownames(preservationData[[permutation]][[2]]) %in% WGCNAobjects[[permutation]]@outlierModules,]
    dir.create("filteredPreservation")
    write.csv(filteredPreservationData[[permutation]], paste0("filteredPreservation/filt", permutation, ".csv"), row.names=F)
    
    save.image("presPermutation.RData")
  }

  z.summary.dist=list()
  module.size=list()
  for(permutation in 1:nPermutations){
    myPerm=filteredPreservationData[[permutation]]
    z.summary.dist[[permutation]]=myPerm$Zsummary[which.min(abs(myPerm$moduleSize-297))]
    module.size[[permutation]]=myPerm$moduleSize[which.min(abs(myPerm$moduleSize-297))]
  }

  if(saveRData) save.image("presPermutation.RData")

  return(data.frame(Z.summary=z.summary.dist, module.size=module.size))
}

minSamplesForPreservation <- function(){

	filterPreservationList <- function(preservationList){
	        preservationList=lapply(preservationList, function(x) x[!rownames(x) %in% c("gold", "WtRand50_00"),])
	        #preservationList=lapply(preservationList, function(x) x[!rownames(x) %in% paste0("WtRand50_", outliers),])
	        preservationList
	}

	minimizePermutation <- function(permutationTable, NSamples){
		modules=rownames(permutationTable[[1]])
		collapsedList=do.call(rbind, permutationTable)
		collapsedList$Module=rep(modules, length(permutationTable))
		meltedDf=melt(collapsedList)
		meltedDf= meltedDf[meltedDf[,2]=="Zsummary.pres",]
		meltedDf[,2]=paste0("N", NSamples)
		meltedDf
	}

	minimizeN24=minimizePermutation(filterPreservationList(permutationN24[[1]]), 24)
	minimizeN15=minimizePermutation(filterPreservationList(permutationN15[[1]]), 15)
	minimizeN12=minimizePermutation(filterPreservationList(permutationN12[[1]]), 12)
	minimizeN10=minimizePermutation(filterPreservationList(permutationN10[[1]]), 10)
	minimizeN8=minimizePermutation(filterPreservationList(permutationN8[[1]]), 8)
	minimizeN5=minimizePermutation(filterPreservationList(permutationN5[[1]]), 5)

	mergedData=rbind(minimizeN24, minimizeN15, minimizeN12, minimizeN10, minimizeN8, minimizeN5)
	colnames(mergedData)=c("Module", "NSamples", "Zsummary")
	mergedData=mergedData[!mergedData$Module %in% c("WtRand50_00", "gold"),]

	pdf("minSamplesForPresBoxplot.pdf", width=10, height=5)
	ggplot(data=mergedData, aes(x=Module, y=Zsummary, color=factor(NSamples, levels=c("N24", "N15", "N12", "N10", "N8", "N5")))) +
	    				labs(title="", y="Preservation Z-summary", x="Module") +
	    				geom_boxplot() + #ylim(min(mergedData$Zsummary)-1, max(mergedData$Zsummary)+1)+
	    				scale_color_manual(values=c("red", "gold", "green", "cyan", "blue", "magenta")) +
	    				guides(color=guide_legend(title="Sample size")) +
	    				annotate("rect", xmin=0, xmax=20.5, ymin = 2, ymax = 10, alpha = .1, fill = "orange") +
	    				annotate("rect", xmin=0, xmax=20.5, ymin = min(mergedData$Zsummary)-1, ymax = 2, alpha = .1, fill = "red") +
	    				theme(panel.background = element_blank(), axis.line = element_line(colour = "black"),
	    					panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
	    					axis.text.x=element_text(angle=45, hjust=1), plot.title=element_text(hjust=0.5),
	    					legend.key = element_rect(fill = "white"))
	dev.off()
}
