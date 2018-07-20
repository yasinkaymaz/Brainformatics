# Brainformatics
This Repository is the collection of developing tools and programs to be potentially used for analyzing single cell NGS data in several Brain Projects at Harvard University.


### Running GSEA on scRNAseq clusters.
This function takes an Seurat object which stores normalized gene expression for every cell in the dataset. Creates required input files, .gct and .cls, out of expression matrix. Then, runs GSEA on between the given conditions for each cluster against each geneset provided.

```
RunGSEAforClusters(SeuratObj = Seurat.obj, Cond1 = "Condition-1", Cond2 = "Condition-2", GeneSet = "geneset_in_gmt_format.gmt", outputDir = 'outputDir/')
```

#### Genesets database:
I took the datasets created for mouse from http://ge-lab.org/gskb/

```
mmDatasets.avail=c('MousePath_Co-expression_gmt.gmt',
                  'MousePath_GO_gmt.gmt',
                  'MousePath_Location_gmt.gmt',
                  'MousePath_Metabolic_gmt.gmt',
                  'MousePath_Other_gmt.gmt',
                  'MousePath_Pathway_gmt.gmt',
                  'MousePath_TF_gmt.gmt',
                  'MousePath_miRNA_gmt.gmt')
```
