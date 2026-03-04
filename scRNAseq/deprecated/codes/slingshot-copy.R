if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("clusterExperiment", version = "3.8")
BiocManager::install("scone", version = "3.8")
BiocManager::install("zinbwave", version = "3.8")
BiocManager::install("kstreet13/slingshot")


install.packages("doParallel")
install.packages("gam")
install.packages("remotes")
install.packages("copula")
install.packages("gsl")
# Bioconductor
library(BiocParallel)
library(clusterExperiment)
library(scone)
library(zinbwave)

# GitHub
library(slingshot)

# CRAN
library(doParallel)
library(gam)
library(RColorBrewer)

set.seed(20)
register(SerialParam())












