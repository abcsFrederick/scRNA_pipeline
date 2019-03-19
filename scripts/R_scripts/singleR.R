#!/usr/bin/env Rscript
# Version: 1.1 beta release
# Date: February 4, 2019
# The script was developed by CCR-SF Bioinformatics Group at FNLCR
# Contact: CCRSF_IFX@nih.gov if you have any questions
# Process seurat 10x pipeline output for identifying cell types

library(reshape2)
library(Seurat)
library(SingleR)

sessionInfo()

args<-commandArgs(TRUE);
show(args);
workdir<-args[1];
rdsPath<-args[2];
genome<-args[3]
markerList<-args[4]

if (genome == "mm10") { 
  species = "Mouse"
}
if (genome == "hg38") {
  species = "Human"
}
geneList <- read.csv(markerList, header=FALSE, row.names=1, stringsAsFactors = FALSE)


#######LOAD SEURAT OBJECT
tenx<-readRDS(rdsPath)


####SET WORKING DIRECTORY
setwd(workdir)


####PLOT GENERATION
plot_generation <- function(singler, type){
  #### TSNE – CLUSTERS FROM SEURAT (RES 0.6)
  
  png(paste0("mainType_", type, "_tsne.png"), height = 8, width =8, units='in', res=300)
  out = SingleR.PlotTsne(singler$singler[[1]]$SingleR.single,
                         singler$meta.data$xy,do.label = T,
                         do.letters = T,
                         labels=singler$seurat@meta.data[singler$singler[[1]]$SingleR.single$cell.names,]$res.0.6,
                         dot.size = 1.3,label.size = 5,alpha=0.8)
  print(out$p)
  dev.off()
  
  #### TSNE WITH ANNOTATIONS
  png(paste0("mainType_", type, "_clusters.png"), height = 8, width =8, units='in', res=300)
  out1 = SingleR.PlotTsne(singler$singler[[2]]$SingleR.single,
                          singler$meta.data$xy,do.label=FALSE,
                          do.letters =T,labels=singler$singler[[2]]$SingleR.single$labels,
                          dot.size = 1.6, font.size = 8)
  print(out1$p)
  dev.off()
  
  ### HEATMAP
  if(type == "true"){
  try({
    png(paste0("mainType_", type, "_heatmap_celltypes.png"), height=7, width=7, units='in', res=300)
    print(SingleR.DrawHeatmap(singler$singler[[2]]$SingleR.single,top.n=25,order.by.clusters = TRUE,
                              clusters = singler$singler[[2]]$SingleR.single$labels))
    dev.off()
    
    png(paste0("mainType_", type, "_heatmap_clusters.png"), height=7, width=7, units='in', res=300)
    print(SingleR.DrawHeatmap(singler$singler[[2]]$SingleR.single,top.n=25, order.by.cluster=TRUE, 
                              clusters=singler$seurat@meta.data[singler$singler[[1]]$SingleR.single$cell.names,]$res.0.6))
    dev.off()
  })}else if (type == "false") {
  	try({
    png(paste0("mainType_", type, "_heatmap_celltypes.png"), height=7, width=7, units='in', res=300)
    print(SingleR.DrawHeatmap(singler$singler[[2]]$SingleR.single,top.n=Inf,order.by.clusters = TRUE,
                              clusters = singler$singler[[2]]$SingleR.single$labels))
    dev.off()
    png(paste0("mainType_", type, "_heatmap_clusters.png"), height=7, width=7, units='in', res=300)
    print(SingleR.DrawHeatmap(singler$singler[[2]]$SingleR.single,top.n=Inf, order.by.cluster=TRUE, 
                              clusters=singler$seurat@meta.data[singler$singler[[1]]$SingleR.single$cell.names,]$res.0.6))
    dev.off()
  })}
  graphics.off()
  
  #### FEATURE PLOT
  try({
    png(paste0("mainType_", type, "_featureplot.png"), height=7, width=7, units='in', res=300)
    print(SingleR.PlotFeature(singler$singler[[2]]$SingleR.single,singler$seurat,
                              plot.feature = -log10(singler$singler[[2]]$SingleR.single.main$pval)))
    dev.off()
  })
  graphics.off()
  
  #### DATA FOR BOX PLOT####
  try({
    df = data.frame(nGene=singler$seurat@meta.data[singler$singler[[1]]$SingleR.single$cell.names,]$nGene,
                  pval=-log10(singler$singler[[2]]$SingleR.single.main$pval),
                  Identity=singler$singler[[2]]$SingleR.single$labels)
  
    png("cell_type_pval_boxplot.png",height = 8, width =12, units='in', res=300)
    print(ggplot(df,aes(x=Identity,y=pval,color=Identity))+geom_boxplot()+
      ylab('-log10(p-value)')+theme(axis.text.x = element_text(angle = 45, hjust = 1) ) +
      geom_hline(yintercept=1.3, linetype="dashed", color = "chartreuse4") + geom_hline(yintercept=1, linetype="dashed", color = "red"))
    dev.off()
  })
  graphics.off()
}

##### SINGLR OBJECT
singler = CreateSinglerObject(tenx@data, annot = NULL, 'project.name', min.genes = 100,
                              technology = "10X", species = species,  normalize.gene.length = F, variable.genes = "de",
                              fine.tune = T, do.signatures = T, clusters = NULL, do.main.types = T,
                              reduce.file.size = T, numCores =16)
singler$seurat = SubsetData(tenx, cells.use=singler$singler[[1]]$SingleR.single$cell.names) # (optional)
singler$meta.data$xy = tenx@dr$tsne@cell.embeddings[singler$singler[[1]]$SingleR.single$cell.names,] # the tSNE coordinates
singler$meta.data$clusters = tenx@meta.data[singler$singler[[1]]$SingleR.single$cell.names,]$res.0.6 # the Seurat clusters (if 'clusters' not provided)
saveRDS(singler, 'singleR_types_true.rds')


tenx@meta.data[singler$singler[[1]]$SingleR.single$cell.names,paste0('singler_', singler$singler[[2]]$about$RefData)] <- 'N/A'
tenx@meta.data[singler$singler[[1]]$SingleR.single$cell.names,paste0('singler_', singler$singler[[2]]$about$RefData)] <- singler$singler[[2]]$SingleR.single$labels[singler$singler[[1]]$SingleR.single$cell.names,]
tenx@meta.data[singler$singler[[1]]$SingleR.single$cell.names,paste0('singler_', singler$singler[[1]]$about$RefData)] <- 'N/A'
tenx@meta.data[singler$singler[[1]]$SingleR.single$cell.names,paste0('singler_', singler$singler[[1]]$about$RefData)] <- singler$singler[[1]]$SingleR.single$labels[singler$singler[[1]]$SingleR.single$cell.names,]
saveRDS(tenx, 'seur_10x_cluster_singler.rds')


plot_generation(singler, "true")
temp <- table(singler$singler[[1]]$SingleR.single$labels[singler$singler[[1]]$SingleR.single$cell.names,],singler$seurat@meta.data[singler$singler[[1]]$SingleR.single$cell.names,'orig.ident'])
names(temp)[1] <- "Cell Type"
write.csv(temp[order(temp, decreasing=TRUE), , drop = FALSE], file=paste0('annotations_', singler$singler[[1]]$about$RefData, '.csv'))
temp <- table(singler$singler[[1]]$SingleR.single.main$labels[singler$singler[[1]]$SingleR.single$cell.names,],singler$seurat@meta.data[singler$singler[[1]]$SingleR.single$cell.names,'orig.ident'])
names(temp)[1] <- "Cell Type"
write.csv(temp[order(temp, decreasing=TRUE), , drop = FALSE], file=paste0('annotations_', singler$singler[[1]]$about$RefData, '_general.csv'))
temp <- table(singler$singler[[2]]$SingleR.single$labels[singler$singler[[1]]$SingleR.single$cell.names,],singler$seurat@meta.data[singler$singler[[1]]$SingleR.single$cell.names,'orig.ident'])
names(temp)[1] <- "Cell Type"
write.csv(temp[order(temp, decreasing=TRUE), , drop = FALSE], file=paste0('annotations_', singler$singler[[2]]$about$RefData, '.csv'))
temp <- table(singler$singler[[2]]$SingleR.single.main$labels[singler$singler[[1]]$SingleR.single$cell.names,],singler$seurat@meta.data[singler$singler[[1]]$SingleR.single$cell.names,'orig.ident'])
names(temp)[1] <- "Cell Type"
write.csv(temp[order(temp, decreasing=TRUE), , drop = FALSE], file=paste0('annotations_', singler$singler[[2]]$about$RefData, '_general.csv'))

#####################
### Plot Genes on clusters
######################
dir.create("gene_list_plots")
for (i in 1:dim(geneList)[1]) {
  genes.use <- as.character(geneList[row.names(geneList)[i],])
  if (sum(genes.use %in% row.names(singler$seurat@data)) > 1) {
    df = data.frame(x=singler$seurat@dr$tsne@cell.embeddings[,1],
                    y=singler$seurat@dr$tsne@cell.embeddings[,2],
                    t(as.matrix(singler$seurat@data[genes.use[genes.use %in% row.names(singler$seurat@data)],])))
    df = melt(df,id.vars = c('x','y'))
    png(paste0("gene_list_plots/",row.names(geneList)[i],".png"),height = 2, width = (length(genes.use[genes.use %in% row.names(singler$seurat@data)])+1), units='in', res=300)
    print(ggplot(df,aes(x=x,y=y,color=value)) +
            geom_point(size=0.5)+scale_color_gradient(low="gray", high="blue") +
            facet_wrap(~variable,nrow=1) +theme_classic()+xlab('')+ylab('')+
            theme(strip.background = element_blank()))
    dev.off()
  }else if(sum(genes.use %in% row.names(singler$seurat@data)) == 1){
  	png(paste0("gene_list_plots/",row.names(geneList)[i],".png"),height = 2, width = 3, units='in', res=300)
  	print(FeaturePlot(singler$seurat, genes.use[genes.use %in% row.names(singler$seurat@data)], cols.use=c("gray", "blue"), no.legend=FALSE))
  	dev.off()
  }
}


