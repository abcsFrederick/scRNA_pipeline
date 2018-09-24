#!/usr/bin/env Rscript
# Version: 1.0 beta release
# Date: September 17, 2018
# The script was developed by CCR-SF Bioinformatics Group at FNLCR
# Contact: CCRSF_IFX@nih.gov if you have any questions
#Process seurat 10x pipeline output for characterization of cluster results


library(clusterProfiler)
library(dplyr)

#These libraries loaded conditionally (see below)
#library(org.Mm.eg.db)
#library(org.Hs.eg.db)

args<-commandArgs(TRUE);
show(args);
workdir<-args[1];
genome<-args[2];
seuratPath<-args[3];
seuratResolutionFile<-args[4];

if (genome == "mm10") {
  library(org.Mm.eg.db)
  orgdb = "org.Mm.eg.db"
  organism = "mmu"
}else if(genome == "hg38"){
  library(org.Hs.eg.db)
  orgdb = "org.Hs.eg.db"
  organism = "hsa"
}


res <- as.double(read.table(seuratResolutionFile))
degs <- readRDS(paste(seuratPath, "/markers_res", res, ".rds", sep=""))

setwd(workdir)

maxGenes <- 100
categories <- 50
clusters <- unique(degs$cluster)

#up <- degs[degs["avg_logFC"] > 0,]
#down <- degs[degs["avg_logFC"] < 0,]
temp <- degs %>% group_by(cluster) %>% top_n(-maxGenes, p_val)

result <- vector("list")
for(i in clusters){
  result[i] <- list(temp$gene[grep(i, temp$cluster)])
}

ids = list()
for(i in clusters){
  test <- bitr(result[[i]], fromType = "SYMBOL", toType = "ENTREZID", OrgDb = orgdb)
  ids[paste('X', i, sep="")] <- list(test$ENTREZID)
}

try({
  ck <- compareCluster(geneCluster = ids, fun = "enrichKEGG", organism=organism)#, OrgDb = "org.Mm.eg.db")
  ht = max(7, min(ceiling(sum(ck@compareClusterResult$Count > 0)/4), ceiling(categories*length(clusters)/4)))
  pdf("enrichKegg.pdf", height=ht)
  print(dotplot(ck, showCategory=categories))
  dev.off()
})
try({
  ck <- compareCluster(geneCluster = ids, fun = "enrichGO", OrgDb = orgdb)
  ht = max(7, min(ceiling(sum(ck@compareClusterResult$Count > 0)/4), ceiling(categories*length(clusters)/4)))
  wt = max(ceiling(max(nchar(ck@compareClusterResult$Description))/8), 7)
  pdf("enrichGO.pdf", height=ht, width=wt)
  print(dotplot(ck, showCategory=categories))
  dev.off()
})
try({
  ck <- compareCluster(geneCluster = ids, fun = "groupGO", ont="BP", OrgDb = orgdb)
  #ht = max(7, min(ceiling(sum(ck@compareClusterResult$Count > 0)/4), ceiling(categories*length(clusters)/4)))
  wt = max(ceiling(max(nchar(as.character(ck@compareClusterResult$Description)))/8), 7)
  pdf("groupGO.pdf", width=wt)
  print(dotplot(ck, showCategory=categories))
  dev.off()
 })
