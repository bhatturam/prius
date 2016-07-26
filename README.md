# Prius
This is an R package for analysis of affected pathways in differential gene expression experiments
For more details refer the preprint From Differentiated Genes to Affected Pathways which can be found at http://biorxiv.org/content/early/2016/02/07/038901.

This package can be installed in R using [devtools](https://cran.r-project.org/web/packages/devtools/index.html) using the function [install_github](http://www.inside-r.org/packages/cran/devtools/docs/install_github) by running the following commands

```R
  library(devtools)
  install_github("bhatturam/prius")
```
For an example of individual and cohort analysis (Datasets used in the paper), please look at demo.R. Refer the manual for documentation
of the R functions in the package.

## NOTE: 
The package TCGAbiolinks, that is used in demo.R to fetch data from TCGA has stopped working. This is because the TCGA data has been moved. (See post here - https://support.bioconductor.org/p/83999/). 
