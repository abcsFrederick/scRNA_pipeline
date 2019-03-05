#!/usr/bin/env Rscript
# Version: 1.1 beta release
# Date: September 17, 2018
# The script was developed by CCR-SF Bioinformatics Group at FNLCR
# Contact: CCRSF_IFX@nih.gov if you have any questions
#Process 10x genomics output with scrnaPipeline for single cell data
#this describes what needs to be done to 10x or other single cell matrix type data in order to run cyclone (where they use ensembl id rather than symbol)

library(scater)
library(dplyr)

args<-commandArgs(TRUE);
show(args);
show(str(args))
path_to_matrix<-args[1];
workdir<-args[2];
genome<-args[3];

####sample ct35-2 located in 10x/kelly/ct35_2

data <- data.frame(read.csv(path_to_matrix, header=TRUE, row.names=1))
setwd(workdir)
#use as input original data frame.
library(MASS)
library(Matrix)
library(RColorBrewer)
library(scales)
library(viridis)
library(scran) # from Bioconductor
species <- "human"
library(org.Hs.eg.db) # from Bioconductor (for human data)
library(org.Mm.eg.db) # from Bioconductor (for mouse data)

if (genome == "mm10") {
        species <- "mouse"
} else if (genome == "hg19") {
  	species  <- "human"
} else if (genome == "hg38") {
	species <- "human"
} else { }

eb_raw <- data

temp_r <- rownames(eb_raw)[which(duplicated(toupper(rownames(eb_raw))))]
temp_r <- lapply(temp_r,function(X) grep(paste0("^",X,"$"),rownames(eb_raw),ignore.case=T))
temp_r <- which(rownames(eb_raw) %in%
                  names(sapply(temp_r,function(X) which.min(apply(eb_raw[X,],1,function(Y) sum(Y>0))))))
if (length(temp_r) > 0) { eb_raw <- eb_raw[-temp_r,] }

eb_raw <- eb_raw[Matrix::rowSums(eb_raw) > 0,]
eb_raw <- eb_raw[,Matrix::colSums(eb_raw) > 0]

cS <- data.frame(libSize=Matrix::colSums(eb_raw),
                 geneDetect=apply(eb_raw,2,function(X) sum(X>0)))
p_hi <- 1e-3 #p-value for filtering doublets
p_lo <- 1e-2 #p-value for filtering poor libraries
fitLibSize <- fitdistr(cS$libSize,"negative binomial")
c_hiL <- qnbinom(p_hi,size=fitLibSize$estimate["size"],
                 mu=fitLibSize$estimate["mu"],lower.tail=F)
c_loL <- qnbinom(p_lo,size=fitLibSize$estimate["size"],
                 mu=fitLibSize$estimate["mu"],lower.tail=T)
fitGeneDetect <- fitdistr(cS$geneDetect,"negative binomial")
c_hiG <- qnbinom(p_hi,size=fitGeneDetect$estimate["size"],
                 mu=fitGeneDetect$estimate["mu"],lower.tail=F)
c_loG <- qnbinom(p_lo,size=fitGeneDetect$estimate["size"],
                 mu=fitGeneDetect$estimate["mu"],lower.tail=T)
                                                                                                                        
temp_doublets <- (cS$libSize > c_hiL) | (cS$geneDetect > c_hiG) #doublets IDed based on high library size or genes detected
temp_crapLibs <- (cS$libSize < c_loL) | (cS$geneDetect < c_loG) #poor libraries IDed based on low library size or genes detected

eb_rawF <- eb_raw[,!(temp_doublets | temp_crapLibs)]
temp_postFgenes <- Matrix::rowSums(eb_rawF) > 0

#eb1p <- newSCESet(countData=eb_rawF)
eb1p <- SingleCellExperiment(assays=list(counts = as.matrix(eb_rawF)), rowData=list(feature_symbol=rownames(eb_rawF)), colData=list(barcode=colnames(eb_rawF)))

#drop mito 
if (species == "human") {
  mitoGenePrefix <- "^MT-"
} else if (species == "mouse") {
  mitoGenePrefix <- "^mt-"
} else { } #Mito gene IDs and cell cycle prediction might need your attention.

eb1p <- calculateQCMetrics(eb1p,feature_controls=list(Mt=grepl(mitoGenePrefix,rownames(eb1p))))
forText_mitoPct <- round(median(eb1p$pct_counts_feature_control))
drop_mitoMads <- 4
drop_mito <- isOutlier(eb1p$pct_counts_feature_control,nmads=drop_mitoMads,type="higher")

eb1 <- eb1p[,!drop_mito]
eb1 <- eb1[rowSums(counts(eb1)) > 0,]

forText_dropMito <- sum(drop_mito)
forText_mitoMads <- drop_mitoMads
forText_genesLost <- nrow(eb1p) - nrow(eb1)
forText_genesRemain <- nrow(eb1)

#run cyclone and get cell states
if (species == "human") {
	anno <- select(org.Hs.eg.db, keys=rownames(eb1), keytype="SYMBOL", column="ENSEMBL")
	cycScores <- cyclone(eb1,gene.names=anno$ENSEMBL[match(rownames(eb1), anno$SYMBOL)],
                       pairs=readRDS(system.file("exdata", "human_cycle_markers.rds", package="scran")))
} else if (species == "mouse") {
	anno <- select(org.Mm.eg.db, keys=rownames(eb1), keytype="SYMBOL", column="ENSEMBL")
	cycScores <- cyclone(eb1,gene.names=anno$ENSEMBL[match(rownames(eb1), anno$SYMBOL)],
                       pairs=readRDS(system.file("exdata", "mouse_cycle_markers.rds", package="scran")))
} else { } 

cycScores$phases <- as.factor(cycScores$phases)

png("CellCycle.png",width=16, height=8, units='in', res=300)

### Visualize cell cycle phase per cell
cycDlibSize <- tapply(eb1$total_counts,cycScores$phases,function(X) density(X))
#cycDgeneDetect <- tapply(eb1$total_features,cycScores$phases,function(X) density(X))
cycDgeneDetect <- tapply(eb1$total_features_by_counts,cycScores$phases,function(X) density(X))

layout(matrix(c(2,1,0,3,5,4,0,6),2),
       widths=c(3.6,.6,3.6,.6),heights=c(.6,3.6))
par(mar=c(3,3,0,0),mgp=2:0)
plot(cycScores$score$G1,cycScores$score$G2M,pch=21,cex=1.2,
     col=viridis(3,.5)[c(3,1,2)][cycScores$phases],
     bg=viridis(3,0.2)[c(3,1,2)][cycScores$phases],
     xlab="G1 score", ylab="G2/M score")
par(mar=c(0,3,0,0))
hist(cycScores$score$G1,breaks=100,col="grey",main=NULL,xaxt="n")
par(mar=c(3,0,0,0))
barplot(hist(cycScores$score$G2M,breaks=100,plot=F)$counts,
        horiz=T,space=0,col="grey",main=NULL,xlab="Frequency")

par(mar=c(3,3,0,0),mgp=2:0)
#plot(total_features~total_counts,data=colData(eb1),pch=21,cex=1.2,
plot(total_features_by_counts~total_counts,data=colData(eb1),pch=21,cex=1.2,
     col=viridis(3,.5)[c(3,1,2)][cycScores$phases],
     bg=viridis(3,0.2)[c(3,1,2)][cycScores$phases],
     xlab="Library Size",ylab="Genes detected",
     main=NULL)
legend("bottomright",bty="n",pch=21,legend=levels(cycScores$phases),xpd=NA,
       col=viridis(3)[c(3,1,2)],pt.bg=viridis(3,0.5)[c(3,1,2)])
par(mar=c(0,3,0,0))
plot(x=NULL,y=NULL,ylab="Density",xaxt="n",
     xlim=range(eb1$total_counts),
     ylim=c(min(sapply(cycDlibSize,function(X) min(X$y))),
            max(sapply(cycDlibSize,function(X) max(X$y)))))
for (x in 1:length(cycDlibSize)) {
  lines(cycDlibSize[[x]],col=viridis(3)[c(3,1,2)][x],lwd=3)
}
par(mar=c(3,0,0,0))
plot(x=NULL,y=NULL,xlab="Density",yaxt="n",
     xlim=c(min(sapply(cycDgeneDetect,function(X) min(X$y))),
            max(sapply(cycDgeneDetect,function(X) max(X$y)))),
     #ylim=range(eb1$total_features))
     ylim=range(eb1$total_features_by_counts))
for (x in 1:length(cycDgeneDetect)) {
  lines(x=cycDgeneDetect[[x]]$y,y=cycDgeneDetect[[x]]$x,col=viridis(3)[c(3,1,2)][x],lwd=3)
}

### Data frame of cell cycle distribution
CellCyclePhase <- sapply(table(cycScores$phases)/length(cycScores$phases),percent)
as.data.frame(rbind(CellCyclePhase))

## Filter out low abundance genes
geneStatsR <- with(rowData(eb1),data.frame(DR=n_cells_by_counts/ncol(eb1),
                                         MDTC=total_counts/n_cells_by_counts,
                                         MTC=total_counts/ncol(eb1),
                                         sumTC=total_counts))
geneStatsR$cellMax <- apply(counts(eb1),1,max)

lowCellNum <- 3 # filter out genes detected in less than this many cells
DRcut <- lowCellNum/ncol(eb1)
drop_lowCell <- geneStatsR$DR < DRcut
eb1F1 <- eb1[!drop_lowCell,]

dev.off()

#make other plots
png("FilteringPlot.png", width=8, height=8, units='in', res=300)
layout(matrix(c(2,1,0,3),2),c(7,1.4),c(1.4,7))
par(mar=c(3,3,0,0),mgp=2:0)
plot(geneDetect~libSize,data=cS,
     pch=21,col=alpha("black",0.2),bg=alpha("black",0.1),cex=1.2,
     xlim=range(cS$libSize),ylim=range(cS$geneDetect),
     xlab="Library Size",ylab="Genes Detected")
points(geneDetect~libSize,cex=1.2,pch=4,col="red",data=cS[temp_doublets | temp_crapLibs,])
legend("topleft",bty="n",inset=c(-.02,0),legend=c(paste("Total genes:",nrow(eb_raw)),
                                                  paste("Post-filter:",sum(temp_postFgenes))))
legend("bottomright",bty="n",pch=4,col="red",
       legend=c(paste(sep="\n","Poor-quality libraries",
                      paste0("(p<",p_lo,"): ",sum(temp_crapLibs))),
                paste(sep="\n","Predicted doublets",
                      paste0("(p<",p_hi,"): ",sum(temp_doublets)))))

par(mar=c(0,3,1,0))
tempD <- density(rnbinom(10000,size=fitLibSize$estimate["size"],mu=fitLibSize$estimate["mu"]))
hist(cS$libSize,breaks=100,freq=F,col="grey",main=NULL,xaxt="n",ylab="Density")
lines(tempD,lwd=2,col=alpha("red",0.5))
abline(v=c_hiL,lty=2,lwd=2,col="darkred")
abline(v=c_loL,lty=2,lwd=2,col="darkred")

par(mar=c(3,0,0,1))
tempD <- density(rnbinom(10000,size=fitGeneDetect$estimate["size"],mu=fitGeneDetect$estimate["mu"]))
tempH <- hist(cS$geneDetect,breaks=100,plot=F)
tempB <- barplot(tempH$density,horiz=T,space=0,col="grey",main=NULL,xlab="Density")
tempSF <- (max(tempB) - min(tempB)) / (max(tempH$mids) - min(tempH$mids))
lines(y=tempD$x * tempSF + (min(tempB) - min(tempH$mids) * tempSF),
      x=tempD$y,lwd=2,col=alpha("red",0.5))
abline(h=c_hiG * tempSF + (min(tempB) - min(tempH$mids) * tempSF),lty=2,lwd=2,col="darkred")
abline(h=c_loG * tempSF + (min(tempB) - min(tempH$mids) * tempSF),lty=2,lwd=2,col="darkred")
dev.off()


#mito filtering plot
png("MitoFiltering.png", width=16, height=8, units='in', res=300)
layout(matrix(c(1,2,4,3,0,5),2),c(4.3,4.2,0.5),c(0.5,4))
par(mar=c(0,3,3,1),mgp=2:0)
plot.new()
title(main="Filtering cells based on\nmitochondrial gene proportion")
par(mar=c(3,3,0,1))
#plot(total_features~total_counts,data=colData(eb1p)[!drop_mito,],
plot(total_features_by_counts~total_counts,data=colData(eb1p)[!drop_mito,],
     pch=21,cex=1.2,col=alpha("black",0.2),bg=alpha("black",0.1),
     xlab="Library Size",ylab="Genes detected",
     #xlim=range(eb1p$total_counts),ylim=range(eb1p$total_features))
     xlim=range(eb1p$total_counts),ylim=range(eb1p$total_features_by_counts))
#points(total_features~total_counts,data=colData(eb1p)[drop_mito,],
points(total_features_by_counts~total_counts,data=colData(eb1p)[drop_mito,],
       pch=21,cex=1.2,col=alpha("red",0.5),bg=alpha("red",0.3))
legend("topleft",bty="n",pch=21,col=alpha("red",0.5),pt.bg=alpha("red",0.3),
       legend=paste(drop_mitoMads,"MADs above median"))
par(mar=c(3,3,0,0))
plot(pct_counts_feature_control~total_counts,data=colData(eb1p),
     pch=21,cex=1.2,col=alpha("black",0.2),bg=alpha("black",0.1),
     xlab="Library Size",ylab="Mitochondrial Transcript Percent")
with(colData(eb1p),abline(h=median(pct_counts_feature_control)+drop_mitoMads*mad(pct_counts_feature_control),
                        lwd=2,lty=2,col=alpha("red",0.5)))
legend("topright",lty=2,lwd=2,col=alpha("red",0.5),bty="n",
       legend=paste(drop_mitoMads,"MADs above median"))
par(mar=c(0,3,0,0))
hist(eb1p$total_counts,breaks=100,col="grey",main=NULL,xaxt="n")
par(mar=c(3,0,0,0))
barplot(hist(eb1p$pct_counts_feature_control,breaks=100,plot=F)$counts,
        horiz=T,space=0,col="grey",main=NULL,xlab="Frequency")
dev.off()

##expression plot
png("ExpressionPlot.png", width=16, height=8, units='in', res=300)
iH <- 101-cut(log10(geneStatsR[order(geneStatsR$cellMax,decreasing=F),"cellMax"]),breaks=100,labels=F)
layout(matrix(c(2,1,4,0,3,0,6,5,8,0,7,0),3),c(3.6,.6,3.6,.6),c(.6,3.6,.4))
par(mar=c(3,3,0,0),mgp=2:0)
plot(log10(MTC)~log10(DR),data=geneStatsR[order(geneStatsR$cellMax,decreasing=F),],
     xlim=log10(range(geneStatsR$DR)),ylim=log10(range(geneStatsR$MTC)),
     pch=21,col=viridis(100,0.5)[iH],bg=viridis(100,0.3)[iH],
     xlab=expression(Log[10]~"Proportion of cells detecting gene"),
     ylab=expression(Log[10]~"Mean transcript count (MTC)"))
points(log10(MTC)~log10(DR),data=geneStatsR[drop_lowCell,],
       pch=4,col=alpha("red",0.5),cex=1.2)
legend("topleft",bty="n",pch=c(4,NA),col=c("red",NA),
       legend=c(paste("Genes in <",lowCellNum,"cells"),
                paste(sum(!drop_lowCell),"genes remain")))
par(mar=c(0,3,.1,0))
hist(log10(geneStatsR$DR),breaks=100,col="grey",main=NULL,xaxt="n")
abline(v=log10(DRcut-1/(2*ncol(eb1))),lty=2,lwd=2,col=alpha("red",0.5))
par(mar=c(3,0,0,.1))
barplot(hist(log10(geneStatsR$MTC),breaks=100,plot=F)$counts,
        horiz=T,space=0,col="grey",main=NULL,xlab="Frequency")
par(mar=c(0.1,3,.5,0))
barplot(rep(1,100),col=viridis(100,begin=1,end=0),space=0,border=NA,axes=F,ylim=c(-3,1))
text(c(1,50,100),rep(-1,3),labels=c(bquote(10^.(log10(min(geneStatsR$cellMax)))),
                                    expression(Log[10]~bold(max)~transcript~count),
                                    bquote(10^.(round(log10(max(geneStatsR$cellMax)),1)))))

par(mar=c(3,3,0,0),mgp=2:0)
plot(log10(MDTC)~log10(DR),data=geneStatsR[order(geneStatsR$cellMax,decreasing=F),],
     xlim=log10(range(geneStatsR$DR)),ylim=log10(range(geneStatsR$MDTC)),
     pch=21,col=viridis(100,0.5)[iH],bg=viridis(100,0.3)[iH],
     xlab=expression(Log[10]~"Proportion of cells detecting gene"),
     ylab=expression(Log[10]~"Mean transcript count of detected genes (MDTC)"))
points(log10(MDTC)~log10(DR),data=geneStatsR[drop_lowCell,],
       pch=4,col=alpha("red",0.5),cex=1.2)
par(mar=c(0,3,.1,0))
hist(log10(geneStatsR$DR),breaks=100,col="grey",main=NULL,xaxt="n")
abline(v=log10(DRcut-1/(2*ncol(eb1))),lty=2,lwd=2,col=alpha("red",0.5))
par(mar=c(3,0,0,.1))
barplot(hist(log10(geneStatsR$MDTC),breaks=100,plot=F)$counts,
        horiz=T,space=0,col="grey",main=NULL,xlab="Frequency")
par(mar=c(0.1,3,.5,0))
barplot(rep(1,100),col=viridis(100,begin=1,end=0),space=0,border=NA,axes=F,ylim=c(-3,1))
text(c(1,50,100),rep(-1,3),labels=c(bquote(10^.(log10(min(geneStatsR$cellMax)))),
                                    expression(Log[10]~bold(max)~transcript~count),
                                    bquote(10^.(round(log10(max(geneStatsR$cellMax)),1)))))
dev.off()

