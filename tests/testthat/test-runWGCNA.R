library(multiWGCNA)
library(ExperimentHub)
eh = ExperimentHub()

# Note: this requires the SummarizedExperiment package to be installed
eh_query = query(eh, c("multiWGCNAdata"))
autism_se = eh_query[["EH8219"]]

# Collect the metadata in the sampleTable
sampleTable = colData(autism_se)

# Randomly sample 2000 genes from the expression matrix
set.seed(1)
autism_se = autism_se[sample(rownames(autism_se), 1500),]

# Check the data
assays(autism_se)[[1]][1:5, 1:5]
sampleTable

# Define our conditions for trait 1 (disease) and 2 (brain region)
conditions1 = unique(sampleTable[,2])
conditions2 = unique(sampleTable[,3])

# Construct the combined networks and all the sub-networks (autism only, controls only, FC only, and TC only)
# The input data can be either a SummarizedExperiment object or a data.frame with genes as rows and samples as columns

test_that("Can run on WGCNA on SummarizedExperiment", {
  autism_networks = constructNetworks(autism_se, sampleTable, conditions1[[1]], conditions2[[1]],
                                      networkType = "unsigned", power = 10,
                                      minModuleSize = 40, maxBlockSize = 25000,
                                      reassignThreshold = 0, minKMEtoStay = 0.7,
                                      mergeCutHeight = 0.10, numericLabels = TRUE,
                                      pamRespectsDendro = FALSE, verbose=3)
  expect_s4_class(autism_networks[[1]], "WGCNA")
})

test_that("Can run on WGCNA on data.frame", {
  autism_networks = constructNetworks(assays(autism_se)[[1]], sampleTable, conditions1[[1]], conditions2[[1]],
                                      networkType = "unsigned", power = 10,
                                      minModuleSize = 40, maxBlockSize = 25000,
                                      reassignThreshold = 0, minKMEtoStay = 0.7,
                                      mergeCutHeight = 0.10, numericLabels = TRUE,
                                      pamRespectsDendro = FALSE, verbose=3)
  expect_s4_class(autism_networks[[1]], "WGCNA")
})
