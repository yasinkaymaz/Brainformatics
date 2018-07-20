
srcdir <- dirname(sys.frame(1)$ofile)
print(srcdir)

RunGSEAforClusters <- function(SeuratObj, Cond1, Cond2, GeneSet='MousePath_GO_gmt.gmt', outputDir=getwd(), ...){
  
  SeuratObj = SetAllIdent(SeuratObj, id = 'cluster.states')
  clusters = sort(unique(as.numeric(SeuratObj@meta.data$tree.ident)))
  pb <- txtProgressBar(min = 0, max = length(clusters), style = 3)
  
  for (i in clusters) {
    #Create GSEA input files
    filePrefix= paste(Cond1,Cond2,"Cluster",i,"ExpMatrix",sep="_")
    CreateGSEAinput(SeuratObj = SeuratObj, Cond1 = Cond1, Cond2 = Cond2, clusterid = i, filePrefix = filePrefix, ...)
    #Run GSEA for the cluster i comparing cells from Cond1 and Cond2    
    RunGSEA(InputPrefix = filePrefix, GeneSet=GeneSet, ...)
    
    setTxtProgressBar(pb,i)
    print(c('GSEA with Cluster', i,' is completed.'))
  }
}


CreateGSEAinput <- function(SeuratObj, Cond1, Cond2, outputDir=getwd(), clusterid, filePrefix, ...){
  
  ident.1 = paste(Cond1,"-",clusterid, sep = '')
  ident.2 = paste(Cond2,'-',clusterid, sep = '')
  #Subset expression data into cells from two conditions:
  sub.data.id1 <- as.data.frame(as.matrix(x = SeuratObj@data[, WhichCells(object = SeuratObj, ident = ident.1)]))
  sub.data.id2 <- as.data.frame(as.matrix(x = SeuratObj@data[, WhichCells(object = SeuratObj, ident = ident.2)]))
  #Store cell numbers in each condition here:
  c1 <- dim(sub.data.id1)[2]
  c2 <- dim(sub.data.id2)[2]
  #merge the datasets from both conditions into a df
  tpm <- cbind( sub.data.id1,sub.data.id2)
  new.df <- cbind(Description="",tpm)
  #get the matrix dimensions
  dimensions <- dim(tpm)
  
  #Create a GCT file: for details: https://software.broadinstitute.org/cancer/software/genepattern/file-formats-guide#GCT
  header1="#1.2"
  header2=paste(dimensions[1],dimensions[2],sep="\t")
  write.table(header1, file=paste(filePrefix,".gct",sep=""), sep="\t", quote=FALSE,col.names=FALSE,row.names=FALSE   )
  write.table(header2, file=paste(filePrefix,".gct",sep=""), sep="\t", quote=FALSE,col.names=FALSE,row.names=FALSE , append=TRUE  )
  write.table(new.df, file=paste(filePrefix,".gct",sep=""), sep="\t", quote=FALSE, append=TRUE,col.names=NA   )
  
  #create a CLS file: for details: https://software.broadinstitute.org/cancer/software/genepattern/file-formats-guide#CLS
  conditions = c(ident.1,ident.2)
  header=paste(dimensions[2], "2", "1",sep=" ")
  line2=paste("#",conditions[1],conditions[2], sep=" ")
  line3=paste( rep(c(conditions[1],conditions[2]), c(c1,c2)),sep = " " )
  write.table(header, file=paste(filePrefix,".cls",sep=""), sep=" ",quote=FALSE,col.names=FALSE,row.names=FALSE)
  write.table(line2,file=paste(filePrefix,".cls",sep=""), sep=" ", quote=FALSE,col.names=FALSE,row.names=FALSE , append=TRUE)
  linex=line3[1]
  for (i in 2:length(line3)){
    linex <- paste(linex,line3[i],sep =" ")}
  write.table(linex,file=paste(filePrefix,".cls",sep=""), sep=" ", quote=FALSE,col.names=FALSE,row.names=FALSE , append=TRUE)
}


RunGSEA <- function(InputPrefix, GeneSet, outputDir=getwd(), ...){
  GSEA.program.location <- paste(srcdir,"/GSEA.1.1.R",sep="")
  source(GSEA.program.location, verbose=F, max.deparse.length=9999)
  library(gskb)
  data(mm_miRNA)
  data(mm_GO)
  data(mm_metabolic)
  data(mm_pathway)
  
  doc.STRING= paste(InputPrefix,substr(GeneSet,1,nchar(GeneSet)-4), sep="_")
  print(doc.STRING)
  print(InputPrefix)
  
  GSEA(   # Input/Output Files :-------------------------------------------
          input.ds =  paste(outputDir,"/",InputPrefix,".gct",sep = ""),           # Input gene expression Affy dataset file in RES or GCT format
          input.cls = paste(outputDir,"/",InputPrefix,".cls",sep = ""),           # Input class vector (phenotype) file in CLS format
          gs.db =   paste(srcdir,"/../GeneSetDatabases/",GeneSet,sep=""),         # Gene set database in GMT format
          output.directory      = paste(outputDir,"/",sep = ""),        # Directory where to store output and results (default: "")
          #  Program parameters :-------------------------------------------------------------------------------------------------------------------------
          doc.string            = doc.STRING,   # Documentation string used as a prefix to name result files (default: "GSEA.analysis")
          non.interactive.run   = F,               # Run in interactive (i.e. R GUI) or batch (R command line) mode (default: F)
          reshuffling.type      = "sample.labels", # Type of permutation reshuffling: "sample.labels" or "gene.labels" (default: "sample.labels" 
          nperm                 = 1000,            # Number of random permutations (default: 1000)
          weighted.score.type   =  1,              # Enrichment correlation-based weighting: 0=no weight (KS), 1= weigthed, 2 = over-weigthed (default: 1)
          nom.p.val.threshold   = -1,              # Significance threshold for nominal p-vals for gene sets (default: -1, no thres)
          fwer.p.val.threshold  = -1,              # Significance threshold for FWER p-vals for gene sets (default: -1, no thres)
          fdr.q.val.threshold   = 0.25,            # Significance threshold for FDR q-vals for gene sets (default: 0.25)
          topgs                 = 20,              # Besides those passing test, number of top scoring gene sets used for detailed reports (default: 10)
          adjust.FDR.q.val      = F,               # Adjust the FDR q-vals (default: F)
          gs.size.threshold.min = 15,              # Minimum size (in genes) for database gene sets to be considered (default: 25)
          gs.size.threshold.max = 500,             # Maximum size (in genes) for database gene sets to be considered (default: 500)
          reverse.sign          = F,               # Reverse direction of gene list (pos. enrichment becomes negative, etc.) (default: F)
          preproc.type          = 0,               # Preproc.normalization: 0=none, 1=col(z-score)., 2=col(rank) and row(z-score)., 3=col(rank). (def: 0)
          random.seed           = 3338,            # Random number generator seed. (default: 123456)
          perm.type             = 0,               # For experts only. Permutation type: 0 = unbalanced, 1 = balanced (default: 0)
          fraction              = 1.0,             # For experts only. Subsampling fraction. Set to 1.0 (no resampling) (default: 1.0)
          replace               = F,               # For experts only, Resampling mode (replacement or not replacement) (default: F)
          save.intermediate.results = F,           # For experts only, save intermediate results (e.g. matrix of random perm. scores) (default: F)
          OLD.GSEA              = F,               # Use original (old) version of GSEA (default: F)
          use.fast.enrichment.routine = T          # Use faster routine to compute enrichment for random permutations (default: T)
  )
  
}

