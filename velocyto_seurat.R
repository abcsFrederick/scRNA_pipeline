#!/usr/bin/env Rscript
# Version: 1.0 beta release
# Date: September 17, 2018
# The script was developed by CCR-SF Bioinformatics Group at FNLCR
# Contact: CCRSF_IFX@nih.gov if you have any questions
#Process single cell sample and loom results to generate velocyto images

library(velocyto.R)
library(pagoda2)
library(igraph)
library(Seurat)

args<-commandArgs(TRUE);
show(args);
workdir<-args[1];
loomFile<-args[2];
seuratRdsFile<-args[3];

ldat <- read.loom.matrices(loomFile)
seur <- readRDS(seuratRdsFile)

setwd(workdir)

emat <- ldat$spliced
rownames(emat) <- make.unique(rownames(emat))
nmat <- ldat$unspliced
rownames(nmat) <- make.unique(rownames(nmat))

colnames(emat) <- gsub("x", "", gsub(paste(strsplit(colnames(emat)[1], ":")[[1]][1], ":", sep=""), "", colnames(emat)))
colnames(nmat) <- gsub("x", "", gsub(paste(strsplit(colnames(nmat)[1], ":")[[1]][1], ":", sep=""), "", colnames(nmat)))

emat <- emat[,rownames(seur@meta.data)]; nmat <- nmat[,rownames(seur@meta.data)];
#res <- grep("res", colnames(seur@meta.data), value=TRUE)[1]
#seur <- SetIdent(seur, ident.use=seur@meta.data[[res]])
seur <- SetIdent(seur, ident.use=seur@meta.data[["res.0.6"]])
cluster.label <- seur@ident
cell.colors <- pagoda2:::fac2col(cluster.label)
emb <- seur@dr$tsne@cell.embeddings
cell.dist <- as.dist(1-armaCor(t(seur@dr$pca@cell.embeddings)))
emat <- filter.genes.by.cluster.expression(emat,cluster.label,min.max.cluster.average = 0.5)
nmat <- filter.genes.by.cluster.expression(nmat,cluster.label,min.max.cluster.average = 0.05)
length(intersect(rownames(emat),rownames(nmat)))
write(length(intersect(rownames(emat),rownames(nmat))), "velocyto_geneOverlap.txt")
fit.quantile <- 0.02
rvel.cd <- gene.relative.velocity.estimates(emat,nmat,deltaT=1,kCells=20,cell.dist=cell.dist,fit.quantile=fit.quantile)

pdf("velocyto.pdf")
show.velocity.on.embedding.cor(emb,rvel.cd,n=300,scale='sqrt',cell.colors=ac(cell.colors,alpha=0.5),cex=0.8,arrow.scale=5,show.grid.flow=TRUE,min.grid.cell.mass=0.5,grid.n=40,arrow.lwd=1,do.par=F,cell.border.alpha = 0.1)
dev.off()
