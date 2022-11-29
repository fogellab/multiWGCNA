#!/usr/bin/env Rscript

convertMouseGeneList <- function(x){
	human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
	mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
	genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
	humanx <- unique(genesV2[, 2])
	print(head(humanx))
	return(humanx)
}

convertHumanGeneList <- function(x){
	human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
	mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
	genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
	humanx <- unique(genesV2[, 2])
	print(head(humanx))
	return(humanx)
}

#' Run goseq gene ontology
#'
#' Performs a goseq gene ontology analysis of a module
#'
#' @author Dario Tommasini
#'
#' @import goseq
#' @export
#'
run_goseq <- function(WGCNAobject, module, organism="mm10", id="geneSymbol"){
  datExpr=WGCNAobject@datExpr
  genes=as.numeric(datExpr$dynamicLabels==module)
  names(genes)=datExpr$X
  pwf=nullp(genes, organism, id, plot.fit=FALSE)
  GO.wall=goseq(pwf, organism, id)
  enriched.GO=GO.wall[p.adjust(GO.wall$over_represented_pvalue, method="BH")<.05,]
  print(head(enriched.GO))
  return(enriched.GO)
}

