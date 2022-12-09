# multiWGCNA
an R package for deep mining gene co-expression networks in multi-trait expression data

```
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("fogellab/multiWGCNA")
```

Alternatively, multiWGCNA can be installed from source by downloading this repository as a zip file, extracting the contents, and performing the following command:

```
R CMD INSTALL --preclean --no-multiarch --with-keep.source multiWGCNA
```

Please refer to the vignette, generalWorkflow.Rmd, for a detailed example of how to use multiWGCNA!
