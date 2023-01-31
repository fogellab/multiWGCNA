# multiWGCNA: an R package for deep mining gene co-expression networks in multi-trait expression data

# Installation 
The multiWGCNA R package can be installed from GitHub like this: 
```
if(!require(devtools)) install.packages("devtools")
devtools::install_github("fogellab/multiWGCNA")
```

# Dependencies

Imports: 
    stringr,
    readr,
    WGCNA,
    ggrepel,
    dplyr,
    reshape2,
    data.table,
    patchwork,
    doParallel,
    scales,
    igraph,
    flashClust,
    filesstrings,
    ggplot2,
    biomaRt,
    dcanr,
    goseq,
    limma,
    vegan,
    GO.db,
    cowplot
Depends: 
    ggalluvial

# Quick Start

Please refer to the vignette, generalWorkflow.Rmd, for a detailed example of how to use multiWGCNA.

For a tutorial using the astrocyte Ribotag data discussed in the manuscript (Tommasini and Fogel. BMC Bioinformatics. 2023), please refer to astrocyte_map2.Rmd.  
