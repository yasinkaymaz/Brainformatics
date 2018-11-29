# Brainformatics
This Repository is the collection of developing tools and programs to be potentially used for analyzing single cell NGS data in several Brain Projects at Harvard University.

```
git clone https://github.com/yasinkaymaz/Brainformatics.git
```

### Running GSEA on scRNAseq clusters.
This function takes an Seurat object which stores normalized gene expression for every cell in the dataset. Creates required input files, .gct and .cls, out of expression matrix. Then, runs GSEA on between the given conditions for each cluster against each geneset provided.

```
RunGSEAforClusters(SeuratObj = Seurat.obj, Cond1 = "Condition-1", Cond2 = "Condition-2", GeneSet = "geneset_in_gmt_format.gmt", outputDir = 'outputDir/')
```

#### Genesets database:
I took the datasets created for mouse from http://ge-lab.org/gskb/

#### An example Rscript to run the function:
Please don't forget to modify the directories before running the script below. Also store the Seurat object file in the working directory.

```{r}
#Please change these according to your paths.

workingDirectory='~/'
source("~/Brainformatics/scripts/functions.R")

### Run
library(Matrix)
library(Seurat)
load(paste(workingDirectory,"Neurons.RObj",sep=""))

setwd(workingDirectory)

#Available mouse datasets are:
mmDatasets=c('MousePath_Co-expression_gmt.gmt',
             'MousePath_GO_gmt.gmt',
             'MousePath_Location_gmt.gmt',
             'MousePath_Metabolic_gmt.gmt',
             'MousePath_Other_gmt.gmt',
             'MousePath_Pathway_gmt.gmt',
             'MousePath_TF_gmt.gmt',
             'MousePath_miRNA_gmt.gmt')

# Use all available mouse genesets one by one:
for (ds in mmDatasets){
    RunGSEAforClusters(SeuratObj = Neurons, Cond1 = "Dominant", Cond2 = "Subordinate", GeneSet = ds, outputDir = './')
    RunGSEAforClusters(SeuratObj = Neurons, Cond1 = "Dominant", Cond2 = "Control", GeneSet = ds, outputDir = './')
    RunGSEAforClusters(SeuratObj = Neurons, Cond1 = "Subordinate", Cond2 = "Control", GeneSet = ds, outputDir = './')
}

```

### Running the script on [Odyssey](https://www.rc.fas.harvard.edu/resources/running-jobs/)

First of all, put the above code in an R script: "scRNAseq-to-GSEA.R". Then, create a bash script with specific resource requests which will submit a job to the scheduler. "job_script.sh":

```{bash}
#!/bin/bash
#SBATCH -n 1
#SBATCH --mem 96000
#SBATCH -p general
#SBATCH -e mrf_err.%j.txt
#SBATCH -o mrf_out.%j.txt
#SBATCH -t 5-24:00

module load R/3.4.2-fasrc01
module load gcc/7.1.0-fasrc01

Rscript scRNAseq-to-GSEA.R
```
Then, submit the script to run the job:

```{bash}
sbatch job_script.sh
```
