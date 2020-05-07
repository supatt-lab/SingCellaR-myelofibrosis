# SingCellaR-myelofibrosis.
This is the shiny application software that can visualise single-cell RNA-seq data. We aim to visualise scRNA-seq data derived from myelofibrosis patients and healthy donors (Psaila and Wang et al., Mol. Cell, 2020). 

# Installation Requirement.
1)shiny
2) SingleCellExperiment
3) Matrix
4) matrixStats
5) ggplot2
6) gridExtra
7) igraph
8) ComplexHeatmap
9) circlize
10) DT

# Important data.

The data can be downloaded from this URL : http://sara.molbiol.ox.ac.uk/public/supatt/shared_dat/Psaila_et_al_MolCell_2020.tar.gz
After downloading, please extract the data and move the "Psaila_et_al_MolCell_2020" into the main package folder!.

# How to run the App.
In R, set up your current working directory and then :

> library(shiny)
> runApp()

# Contact

supat.thongjuea@imm.ox.ac.uk


