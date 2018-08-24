#Please change these according to your paths.
.libPaths(c(.libPaths(),"/n/home13/yasinkaymaz/biotools/Rlibs/"))
workingDirectory='/n/home13/yasinkaymaz/LabSpace/testdata/DulacLab/GSEA/MDT_mn'
source("~/codes/Brainformatics/scripts/functions.R")

### Run
library(Matrix)
library(Seurat)
load(paste(workingDirectory,"mn.RObj",sep=""))
#load(paste(workingDirectory,"ACC.Neurons.RObj",sep=""))

setwd(workingDirectory)


mmDatasets=c("Calvin_manual_genesets.gmt")

for (ds in mmDatasets){
  
  for (exp.Seu.obj in c(mn)){
    
    RunGSEAforClusters(SeuratObj = exp.Seu.obj, Cond1 = "Dominant", Cond2 = "Subordinate", GeneSet = ds, outputDir = './')
    RunGSEAforClusters(SeuratObj = exp.Seu.obj, Cond1 = "Dominant", Cond2 = "Control", GeneSet = ds, outputDir = './')
    RunGSEAforClusters(SeuratObj = exp.Seu.obj, Cond1 = "Subordinate", Cond2 = "Control", GeneSet = ds, outputDir = './')
    
  }
  
}