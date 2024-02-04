# multiWGCNA: an R package for deep mining gene co-expression networks in multi-trait expression data

The multiWGCNA R package builds on the existing weighted gene co-expression 
network analysis (WGCNA) package by extending workflows to expression data with 
two dimensions. multiWGCNA is especially useful for the study of 
disease-associated modules across time or space. For more information, please 
see the multiWGCNA paper available at https://doi.org/10.1186/s12859-023-05233-z. 

# Installation 

The multiWGCNA R package can be installed from Bioconductor like this: 
```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("multiWGCNA")
```

The development version of multiWGCNA can be installed from GitHub like this:
```
if (!require("devtools", quietly = TRUE))
    install.packages("devtools")
    
devtools::install_github("fogellab/multiWGCNA")
```

# Vignettes

We recommend running through both of the vignettes before applying multiWGCNA to your own data:

* The autism_full_workflow.Rmd vignette provides a quick example of how to use multiWGCNA.
* The astrocyte_map2.Rmd vignette provides a more in-depth tutorial discussing the preservation analyses from the manuscript and functions for visualization. 

# Citation

To cite multiWGCNA in publications, please use:

  Tommasini, D, Fogel, BL (2023). multiWGCNA: an R package for deep mining gene
  co-expression networks in multi-trait expression data. BMC Bioinformatics, 24,
  1:115.

For LaTeX users, a BibTeX entry is available here: 

```
  @Article{,
    title = {multiWGCNA: an R package for deep mining gene co-expression 
    networks in multi-trait expression data},
    author = {Dario Tommasini and Brent L. Fogel},
    journal = {BMC Bioinformatics},
    year = {2023},
    volume = {24},
    number = {1},
    pages = {115},
  }
```
