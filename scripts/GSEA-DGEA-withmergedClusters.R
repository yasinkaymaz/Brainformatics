#This script is used to generate GSEA and DGEA results after merging some of the clusters.
library(Matrix)
library(Seurat)
library(dplyr)
library(tidyverse)
source(here::here("scripts/functions.R"))# you can get the functions from https://github.com/yasinkaymaz/Brainformatics.git

load("mn.RObj")

head(mn@meta.data)
PlotClusterTree(mn)
node.scores <- AssessNodes(mn)
print(node.scores)
mn.merged <- mn
#manually pick the nodes to merge: Decided with Adam and Eric.
nodes.to.merge=c(44,46)
for (n in nodes.to.merge){
  mn.merged <- MergeNode(mn.merged, n, rebuild.tree =T)
}

PlotClusterTree(mn.merged)

newclusterids <- levels(mn.merged@ident)

newids <- data.frame( newids=paste(mn.merged@meta.data$States, mn.merged@ident, sep = '-'),  row.names = names(mn.merged@ident) )
mn.merged@meta.data$newids <- as.vector(newids[,1])

mn.merged = SetAllIdent(mn.merged, id = "newids")

save(mn.merged, file="mn.merged.RObj")#After merging Nodes: 44 46


mmDatasets=c("Calvin_manual_genesets.gmt")
i=1
exp.Seu.obj=mn.merged

MDTclusters=c(1, 13, 15)
for (ds in mmDatasets){
  for (exp.Seu.obj in c(mn.merged)){
    i=i+1
    RunGSEAforClusters(SeuratObj = exp.Seu.obj, Cond1 = "Dominant", Cond2 = "Subordinate", Clusteridlist = MDTclusters, GeneSet = ds, outputDir = './')
    RunGSEAforClusters(SeuratObj = exp.Seu.obj, Cond1 = "Dominant", Cond2 = "Control", Clusteridlist = MDTclusters, GeneSet = ds, outputDir = './')
    RunGSEAforClusters(SeuratObj = exp.Seu.obj, Cond1 = "Subordinate", Cond2 = "Control", Clusteridlist = MDTclusters, GeneSet = ds, outputDir = './')

    Sig.Enrichment.ESTable <- SummarizeGSEAoutputs(GSEAoutputDir = "./")

    DEGs.dc <- RunDGEA(SeuratObj = exp.Seu.obj, Cond1 = "Dominant", Cond2 = "Control", Clusteridlist = MDTclusters)
    DEGs.sd <- RunDGEA(SeuratObj = exp.Seu.obj, Cond1 = "Dominant", Cond2 = "Subordinate", Clusteridlist = MDTclusters)
    DEGs.sc <- RunDGEA(SeuratObj = exp.Seu.obj, Cond1 = "Subordinate", Cond2 = "Control", Clusteridlist = MDTclusters)

    #Combine results:
    DEGs <- dplyr::bind_rows(DEGs.dc, DEGs.sd, DEGs.sc, .id='id')
    DEGs$fdr <- p.adjust(DEGs$p_val, method = 'fdr')
    DEGs$id <- factor(DEGs$id)
    levels(DEGs$id) <- list('Dominant_Control'='1', 'Dominant_Subordinate'='2', 'Subordinate_Control'='3')
    DEGs <- arrange(DEGs,fdr)
    save(DEGs,file = paste("DEGs",i,"Rdata",sep="."))
    merged <- inner_join(Sig.Enrichment.ESTable, DEGs, by = c("GENE" = "gene", "Comparison" = "id", "cluster" = "cluster" )) %>% as.data.frame()
    write_tsv(merged,file.path(paste("GSEA-DGE-mergedtables.significantEnrichments",i,"txt",sep=".")))
  }
}
