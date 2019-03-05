#!/usr/bin/env Rscript
# Version: 1.2 beta release
# Date: January 9, 2019
# The script was developed by CCR-SF Bioinformatics Group at FNLCR
# Contact: CCRSF_IFX@nih.gov if you have any questions
#Process 10x genomics output with Seurat package for single cell RNA-seq data
#usage: Final_Seurat.R /workingdirectory/ /path/to/10x/output/ sample_name genome

library(dplyr)
library(Matrix)
library(Seurat)
library(URD)
library(cluster)
library(scClustViz)

roundUp <- function(x,to=10)
{
  to*(x%/%to + as.logical(x%%to))
}

args<-commandArgs(TRUE);
show(args);
show(str(args))
workdir<-args[1];
path<-args[2];
name<-args[3];
genome<-args[4]

data <- Read10X(path)
setwd(workdir)


#make project name 10x_name of file (from arg[1])
seur <- CreateSeuratObject(raw.data = data, min.cells = 3, min.genes = 100, project = name)
seur <- NormalizeData(object=seur)

#make a csv with the barcodes and how many genes per cell
#genematrix <- as.matrix(data)
genematrix <- as.matrix(seur@data)
write.csv(genematrix, file="Gene_Expression_Matrix.csv", quote=FALSE, row.names=TRUE, col.names=TRUE)
tmp <- data.frame(colSums(genematrix !=0))
tmpnames <- row.names(tmp)
tmpall <- data.frame(cbind(tmpnames, tmp))
colnames(tmpall) <- c("Barcodes", "Genes")
write.csv(tmpall, file="Number_of_genes_per_barcode.csv", quote=FALSE, row.names=TRUE)

#make edited plot
x <- tmpall$Genes
mx <- mean(x)
med <- median(x)
d<- density(x)
h<-hist(x, breaks=20, col="grey",  labels = TRUE, yaxt='n', xlab = "Number of Genes", ylab = "Number of Cells", main = name)
#xmin <- min(h$breaks)
xmin <- 0
xmax <- max(roundUp(h$breaks, 500))
png("GenesPerBarcodeEditedPlot.png", height=7, width=7, units='in', res=300)
h<-hist(x, breaks=20, col="grey",  labels = TRUE, yaxt='n', xlab = "Number of Genes", ylab = "Number of Cells", main = name, xlim=c(xmin,xmax), xaxp=c(xmin, xmax, xmax/500))
abline(v=mx, lwd=2, col = "blue")
abline(v=med, lwd=2, col = "orange")
par(new = T)
#h<-hist(x,breaks=20, col="grey",freq = FALSE, yaxt='n', ann=FALSE)
h<-hist(x,breaks=20, col="grey",freq = FALSE, yaxt='n', ann=FALSE, xlim=c(xmin,xmax), xaxp=c(xmin, xmax, xmax/500))
lines(d, col = "red", lwd =2)
mtext(paste("Mean = ", round(mx,2), sep=""), side =3, adj = 0.24,  col= "blue", lwd =2)
#mtext(round(mx,2), side =3, adj = 0.35, col= "blue")
mtext(paste("Median = ", round(med,2), sep=""), side =3, adj = 0.65,  col= "orange")
#mtext(round(med,2), side =3, adj = 0.76, col= "orange")
par (new = F)
dev.off()

#make barcodes per genes
tmpl <- data.frame(rowSums(genematrix !=0))
tmplnames <- row.names(tmpl)
tmplall <- data.frame(cbind(tmplnames, tmpl))
colnames(tmplall) <- c("Genes", "Cells")
write.csv(tmplall, file="Number_of_cells_per_gene.csv", quote=FALSE, row.names=TRUE)

#remove 0's before plotting from tmplall

tmplallplot <- tmplall[tmplall$Cells > 0,]
x <- tmplallplot$Cells
y <- nrow(tmplallplot)
mx <- mean(x)
med <- median(x)
d<- density(x)
h<-hist(x, breaks=20, col="grey",  labels = TRUE, yaxt='n', xlab = "Number of Cells", ylab = "Number of Genes", main = name)
#xmin <- min(h$breaks)
xmin <- 0
xmax <- max(roundUp(h$breaks, 100))
png("CellsPerGenePlot.png", height=7, width=7, units='in', res=300)
h<-hist(x, breaks=20, col="grey",  labels = TRUE, yaxt='n', xlab = "Number of Cells", ylab = "Number of Genes", main = name, xlim=c(xmin, xmax), xaxp=c(xmin, xmax, xmax/100))
abline(v=mx, lwd=2, col = "blue")
abline(v=med, lwd=2, col = "orange")
par(new = T)
h<-hist(x,breaks=20, col="grey",freq = FALSE, yaxt='n', ann=FALSE, xlim=c(xmin, xmax), xaxp=c(xmin, xmax, xmax/100))
lines(d, col = "red", lwd =2)
mtext(paste("Mean = ", round(mx,2), sep=""), side =3, adj = 0.1,  col= "blue", lwd =2)
#mtext(round(mx,2), side =3, adj = 0.25, col= "blue")
mtext(paste("Median = ", round(med,2), sep=""), side =3, adj = 0.4,  col= "orange")
#mtext(round(med,2), side =3, adj = 0.5, col= "orange")
mtext(paste("Total Genes = ", y, sep=""), side =3, adj = 0.80, col = "purple")
#mtext(y, side =3, adj = 0.88, col= "purple")
par (new = F)
dev.off()

#find mitochondrial genes
if (genome == "mm10") {
	mito.genes <- grep("^mt-", x=rownames(x=seur@data), value = T)
	percent.mito <- colSums(as.matrix(seur@data[mito.genes, ]))/colSums(as.matrix(seur@data))
}
if (genome == "hg19") {
 	mito.genes <- grep("^MT-", x=rownames(x=seur@data), value = T)
	percent.mito <- colSums(as.matrix(seur@data[mito.genes, ]))/colSums(as.matrix(seur@data))
}
if (genome == "hg38") {
  mito.genes <- grep("^MT-", x=rownames(x=seur@data), value = T)
	percent.mito <- colSums(as.matrix(seur@data[mito.genes, ]))/colSums(as.matrix(seur@data))
}
seur <- AddMetaData(seur, percent.mito, "percent.mito")

#write plots to one PDF
pdf("VlnPlot.pdf")
par(mfrow = c(10,1))
VlnPlot(seur, c("nGene", "nUMI", "percent.mito"), nCol = 3)
dev.off()

#new plot
pdf("MitoPlot.pdf")
GenePlot(seur, "nUMI", "percent.mito")
dev.off()

#We filter out cells that have unique gene counts over 6000 or higher mitochondrial or ribosomal protein percentages
seur <- SubsetData(seur, subset.name = "nGene", accept.high = 6000)
seur <- SubsetData(seur, subset.name = "percent.mito", accept.high = 0.05)

pdf("VlnPlot-Filtered.pdf")
par(mfrow = c(10,1))
VlnPlot(seur, c("nGene", "nUMI", "percent.mito"), nCol = 3)
dev.off()

pdf("Gene_UMI_plot.pdf")
GenePlot(seur, "nUMI", "nGene")
dev.off()

write.csv(as.matrix(seur@raw.data[,colnames(seur@data)]), file="Filtered_Gene_Expression_Matrix.csv", quote=FALSE, row.names=TRUE, col.names=TRUE)
write.csv(as.matrix(seur@data), file="Filtered_Normalized_Gene_Expression_Matrix.csv", quote=FALSE, row.names=TRUE, col.names=TRUE)

saveRDS(seur, file = "seur_10x_preprocessed_object.rds")



###URD
#inputTags <- as.matrix(read.csv(expressionFile, row.names = 1))
mat1 <- as.matrix(seur@data)
test <- createURD(count.data = mat1, min.cells=3, min.counts=3)
test <- calcPCA(test, mp.factor = 2)
write.table(test@pca.sig,"URD.txt")
pdf("URD.pdf")
pcSDPlot(test)
dev.off()

numPCs <- sum(test@pca.sig)

##Continue Seurat analysis

seur@misc[["udr_nPCs"]] <- sum(test@pca.sig)

seur <- FindVariableGenes(object = seur, x.low.cutoff = 0, y.cutoff = 0.8)
#regress out mitochondrial (unwanted source of variation)
#seur <- RegressOut(seur, latent.vars = c("nUMI", "percent.mito"))
seur <- ScaleData(object = seur, vars.to.regress = c("percent.mito", "nUMI"), genes.use = seur@var.genes, model.use = "negbinom")


#plot variable genes
#you have to run this before runnning PCA function

#cannot run this function without above MeanVarPlot function
#Perform linear dimensional reduction
seur <- RunPCA(seur, pc.genes = seur@var.genes, pcs.compute=max(20,numPCs), do.print = TRUE, pcs.print = 5, genes.print = 5)
seur <- ProjectPCA(seur)
#PrintPCA(seur, pcs.print = 1:5, genes.print = 5, use.full = TRUE)
pdf("VizPCAPlot.pdf")
for (i in c(1:10)) {
  j = i + 1
  VizPCA(seur, i:j)
}
dev.off()

#plot PCA
pdf("AllPCAPlot.pdf")
for (i in c(1:10)) {
  PCAPlot(seur, i, i+1)
}
dev.off()

#PC heatmap
pdf("PC_HeatmapPlot.pdf")
for (i in c(1:10)) {
  PCHeatmap(seur, pc.use = i, cells.use = 100, do.balanced = TRUE)
}
dev.off()

#make PC elbow plot
pdf("PC_ElbowPlot.pdf")
PCElbowPlot(seur)
dev.off()

#determine statistically significant PCs
#seur <- JackStraw(seur, num.replicate = 100, do.print = FALSE)
#seur <- JackStraw(seur, num.replicate = 100, display.progress = FALSE)
#pdf("JackStrawPlot.pdf")
#JackStrawPlot(seur, PCs = 1:18)
#JackStrawPlot(seur, PCs = 1:(nPC+5))
#dev.off()

#run second steps
resolutions <- c(0.1, 0.3, 0.6, 0.8)

#seur <- RunTSNE(seur, dims.use=1:16, do.fast=T)
seur <- RunTSNE(seur, dims.use=1:numPCs, do.fast=T)
write.csv(seur@dr$tsne@cell.embeddings, file = "tSNECoordinates.csv")
seur <- RunUMAP(seur, dims.use=1:numPCs)
write.csv(seur@dr$umap@cell.embeddings, file = "UMAPCoordinates.csv")

runRes <- c()
tsnePlots <- list()
umapPlots <- list()
for (res in resolutions) {
  seur <- FindClusters(seur, dims.use=1:numPCs, resolution = res, print.output = 0, save.SNN = T)
  tsne <- TSNEPlot(seur) + ggtitle(paste(numPCs,"PCs_res", res, sep="")) +
    theme(plot.title = element_text(hjust = 0.5))
  tsnePlots[[as.character(res)]] <- tsne
  png(paste("TSNEPlotwith",numPCs,"PCs_", res, ".png", sep=""), height=7, width=7, units='in', res=300)
  print(tsne)
  dev.off()
  umap <- DimPlot(seur, reduction.use="umap") + ggtitle(paste(numPCs,"PCs_res", res, sep="")) +
    theme(plot.title = element_text(hjust = 0.5))
  umapPlots[[as.character(res)]] <- umap
  png(paste("UMAPPlotwith",numPCs,"PCs_", res, ".png", sep=""), height=7, width=7, units='in', res=300)
  print(umap)
  dev.off()

  try({
    seur.markers <- FindAllMarkers(object = seur, thresh.use = 0.25, only.pos=TRUE)
    write.csv(seur.markers %>% group_by(cluster) %>% top_n(-100,
                                                           p_val), paste("top100markers_pc", numPCs, "_res", res, ".csv", sep = ""))
    saveRDS(seur.markers, paste("markers_res", res, ".rds", sep = ""))
    runRes <- append(runRes, res)})
}

#save object
saveRDS(seur, file = "seur_10x_cluster_object.rds")

pdf("TSNEPlots.pdf")
for (res in tsnePlots){
  print(res)
}
dev.off()

pdf("UMAPPlots.pdf")
for (res in umapPlots){
  print(res)
}
dev.off()


##Create Silhoutte Plots
for (res in runRes){
  coord <- seur@dr$pca@cell.embeddings[,1:numPCs]
  seur <- SetIdent(seur, ident.use=seur@meta.data[[paste("res.", res, sep="")]])
  clusters <- seur@ident
  d <- dist(coord, method="euclidean")
  sil<-silhouette(as.numeric(clusters), dist=d)
  #silPlot <- recordPlot()
  pdf(paste0("SilhouettePlot_res",res,".pdf"))#, height=7, width=7, units='in', res=300)
  plot(sil, col=as.factor(clusters[order(clusters, decreasing=FALSE)]), main=paste("Silhouette plot of Seurat clustering - resolution ", res, sep=""), lty=2)
  abline(v=mean(sil[,3]), col="red4", lty=2)
  dev.off()
}

##Remove resolutions that failed marker generation
print(runRes)
for (res in setdiff(resolutions, runRes)){
  seur@meta.data[paste("res.", res, sep="")] <- NULL
}

write(min(runRes), "minRes.txt")


##Generate scClustViz object
tryCatch({
  data_for_scClustViz <- readFromSeurat(seur)
  DE_for_scClustViz <- clusterWiseDEtest(data_for_scClustViz,exponent=exp(1))
  save(data_for_scClustViz,DE_for_scClustViz,file="for_scClustViz.RData")
},
error = function(cond) {
  print(cond)
})
