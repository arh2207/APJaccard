---
title: "APJaccard Tutorial"
authors: 
 - Andrew Howe
 - Kelly Jones
---
**Authors:** Andrew Howe, Kelly Jones
**Contacts:**

* Andrew Howe: arh2207@columbia.edu
* Kelly Jones: kaj2165@columbia.edu

### Overview

A package with functions that facilitate analysis of AP Clustering solutions
using the Jaccard similarity index

### Installation

Here's how you can install the APJaccard package:

```
install.packages(c("cluster", "ggplot2", "devtools", "Seurat", 
                   "pheatmap", "BiocManager", "RColorBrewer", 
		   "stringr", "apcluster", "igraph", "Matrix",
                   "tidyverse", "mcclust", "philentropy"))
BiocManager::install("viper")
BiocManager::install("biomaRt")
devtools::install_github("JEFworks/MUDAN")
devtools::install_github(repo = "califano-lab/PISCES", force = TRUE, build_vignettes = TRUE)
devtools::install_github('https://github.com/arh2207/APTestPipeline')
devtools::install_github('https://github.com/jchiquet/aricode')
```
=======
### References

1. 	Frey, J. and Dueck, D. (2007) *Clustering by Passing Messages Between Data Points*. Science, 135, 972-6.
2. 	Kozakura, Y. et al. (2017) *Comparison of Methods for Single-Cell Transcrip-tome Analysis Information Processing*. Society of Japan Technical Report, BIO-51, 1-6.

#### Acknowledgements

Philippe Chlenski - TA of Computational Genomics at Columbia University 
Dr. Itsik Pe'er - Professor of Computational Genomics at Columbia University
