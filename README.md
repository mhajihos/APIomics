# APIomics
## Overview
APIomics is a bioinformatics analysis pipeline designed to process and analyze high-throughput expression data. This tool enables users to perform preprocessing, differential expression analysis, gene set enrichment, regulatory network analysis, pubmed and clinical trial search.

## Install Bioconductor dependencies
```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("limma", "edgeR", "clusterProfiler", "org.Hs.eg.db", "ComplexHeatmap","impute", "preprocessCore","DOSE"))
```

## Install the Package
```
install.packages("devtools") 
devtools::install_github("mhajihos/APIomics",force = TRUE)
```

## Load the package
```
library(APIomics)
```

## Launch the app
```
APIomics()
```


![Flowchart](https://github.com/mhajihos/APIomics/blob/master/www/FLowchart.tif)

## Support
If you need help, please refer to the user guide or contact support. Email:morteza.hajihosseini@appliedpharma.ca


