# multiWGCNA: an R package for deep mining gene co-expression networks in multi-trait expression data

# Installation 
The multiWGCNA R package can be installed from GitHub like this: 
```
if(!require(devtools)) install.packages("devtools")
devtools::install_github("fogellab/multiWGCNA")
```

If you would like to build the vignettes, you must specify build_vignettes = TRUE. 
```
if(!require(devtools)) install.packages("devtools")
devtools::install_github("fogellab/multiWGCNA", build_vignettes = TRUE)
```


# Dependencies

multiWGCNA has been tried and tested successfully using the following dependency versions. multiWGCNA may not work for older or newer versions of these dependencies. 

stringr_1.4.0
readr_2.0.2
WGCNA_1.70-3
vegan_2.5-7
dplyr_1.0.7
reshape2_1.4.4
data.table_1.14.6
patchwork_1.1.1
doParallel_1.0.17
scales_1.2.1
igraph_1.2.7
flashClust_1.01-2
filesstrings_3.2.2
ggplot2_3.4.0
biomaRt_2.46.3
dcanr_1.6.0
goseq_1.42.0
readr_2.0.2
vegan_2.5-7 
GO.db_3.12.1
cowplot_1.1.1
ggalluvial_0.12.3

# Quick Start

Please refer to the vignette, generalWorkflow.Rmd, for a detailed example of how to use multiWGCNA.

For a tutorial using the astrocyte Ribotag data discussed in the manuscript (Tommasini and Fogel. BMC Bioinformatics. 2023), please refer to astrocyte_map2.Rmd.  
