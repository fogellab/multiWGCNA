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

1. stringr_1.4.0
2. readr_2.0.2
3. WGCNA_1.70-3
4. vegan_2.5-7
5. dplyr_1.0.7
6. reshape2_1.4.4
7. data.table_1.14.6
8. patchwork_1.1.1
9. doParallel_1.0.17
10. scales_1.2.1
11. igraph_1.2.7
12. flashClust_1.01-2
13. filesstrings_3.2.2
14. ggplot2_3.4.0
15. biomaRt_2.46.3
16. dcanr_1.6.0
17. goseq_1.42.0
18. readr_2.0.2
19. vegan_2.5-7 
20. GO.db_3.12.1
21. cowplot_1.1.1
22. ggalluvial_0.12.3

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

![My Image](images/drawmultiWGCNA.png)

# Citation
