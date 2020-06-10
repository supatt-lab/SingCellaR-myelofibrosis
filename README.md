# SingCellaR-myelofibrosis.
This is the shiny application software that can visualise single-cell RNA-seq data. We aim to visualise scRNA-seq data derived from myelofibrosis patients and healthy donors (Psaila and Wang et al., Mol. Cell, 2020). 

# Installation requirement.
1) shiny <p>
2) SingleCellExperiment <p>
3) Matrix <p>
4) matrixStats <p>
5) ggplot2 <p>
6) gridExtra <p>
7) igraph <p>
8) ComplexHeatmap <p>
9) circlize <p>
10) DT

# Important data.

The data can be downloaded from this URL : http://sara.molbiol.ox.ac.uk/public/supatt/shared_data/projects.tar.gz <p>

Please extract the data and move the folder "projects" into the main package folder!.

# How to run the App.
In R, set up your current working directory to the main package folder and then :

> library(shiny) <p>
> runApp()

You would have to click "Load" button to load the R object into the software. The process of loading might take a bit of time depending on the capacity of the computer. 


# Contact

supat.thongjuea@imm.ox.ac.uk


