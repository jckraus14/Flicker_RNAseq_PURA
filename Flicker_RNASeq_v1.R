# ---- Preliminaries ---- 
#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# 
# 
# install.packages("readxl")
# install.packages("backports")
# library("backports")
# install.packages("devtools")
# library("devtools")
# install.packages("xlsx")
# library("xlsx")
# install.packages("heatmap3")
# install.packages("compare")
# install.packages("gridGraphics")
# BiocManager::install("GSVA")
# install.packages("ggfortify") #for PCA
# install_github("kevinblighe/EnhancedVolcano") 
# install.packages("scales")
# library("scales")
# install.packages("ggplot2")
# install.packages("Rmisc") 
# install.packages("pillar") 


options(java.parameters = "-Xmx8g")
library("fastcluster")
library("gplots")
library("dendextend")
library(reshape)
library(matrixStats)
library("readxl")
library("heatmap3")
library("compare")
library("gridGraphics")
library("ggplot2")
library("Rmisc") 
library("ggpubr")
#library("GSVA")
library("RColorBrewer")
library("Hmisc")
#library("xlsx")
library("ggfortify")#for PCA
#library("EnhancedVolcano")
library('dplyr')
library("gplots")
source("summarySE.R")
library("DESeq2")
library(EnhancedVolcano)

# ---- Set heatmap3 color bar parameters ---- 
breakBarColors=c(-200,seq(-1.5, 1.5, 0.01),200) #Outside numbers clip outliers. This is for zscoring.
barColors = colorpanel(length(breakBarColors)-1, "snow4", "white", "mediumorchid1")

breakBarColorsCor=c(-200,seq(-1, 1.5, 0.01),200) #Outside numbers clip outliers. This is for zscoring.
barColorsCor = colorpanel(length(breakBarColorsCor)-1, "black", "white", "orange")


# ---- Read and Plot RNAseq Data ----
dataIn=read.csv("RNAseq data 21days_raw counts.csv", header = TRUE)

yList=dataIn[, 2:(ncol(dataIn)-1)]

y=matrix(unlist(yList), ncol=ncol(yList), byrow=F)
rownames(y)=toupper(dataIn$geneName)
colnames(y)=colnames(dataIn)[2:(ncol(dataIn)-1)]

uniqueGenes=unique(rownames(y)) 

#Chose just one of the repeate genes
indU=match(uniqueGenes,rownames(y)) #which of the gene pairs is selected for repeat genes?

y=y[indU,]

#Remove lowly expressed genes
indKeepGenes=which(rowSums(y>=4)>=5)
y=y[indKeepGenes,]




pheno=c(rep("A",3),rep("B",3),rep("C",3))





colorBar=matrix(,nrow=length(pheno), ncol=1)

colorBar[pheno=="A"]="black"
colorBar[pheno=="B"]="slateblue4"
colorBar[pheno=="C"]="orangered4"


# breakBarColorsLPS=c(-1,seq(1, 4, 0.1))*10 #Outside numbers clip outliers. This is for z-scoring.
# barColorsAgeVec = c("black",colorpanel(length(breakBarColorsLPS)-2, "white", "orangered"))
# indMatchColor=match(round(dataPhenoIn$X50ug.ml.LPS[-37]*10),round(breakBarColorsLPS))
# barColorsLPS=barColorsAgeVec[indMatchColor]
# 
# breakBarColorsLPS=c(-1,seq(1, 4, 0.1))*10 #Outside numbers clip outliers. This is for z-scoring.
# barColorsAgeVec = c("black",colorpanel(length(breakBarColorsLPS)-2, "white", "orangered"))
# indMatchColor=match(round(dataPhenoIn$X100nM.PMA[-37]*10),round(breakBarColorsLPS))
# barColorsPMA=barColorsAgeVec[indMatchColor]
# 

colorBarCombined=c(colorBar)
#colnames(colorBarCombined)=c("Group")


#pdfPath=paste("./FiguresOut/nocluster.pdf",sep='')
#pdf(pdfPath, width=6,height=8,pointsize = 14)
heatmap3(y, ColSideColors =colorBarCombined, ColSideWidth=1, 
         col=barColors, breaks=breakBarColors, legendfun=function()showLegend(legend=c(NA),col=c(NA),cex=1.5), 
         Rowv=NA, Colv=NA,  scale="row",
         labRow = NA, labCol = NA) 
#dev.off()

yZ=t(apply(y,1,scale))
rownames(yZ)=rownames(y)

cort=cor(t(yZ), method="pearson")

hr = hclust(as.dist(1-cort), method="complete")

cutThreshCTD=1.95 #should be below max(hr$height)
myclCTD = cutree(hr, h=cutThreshCTD); 
#mycolhcCTD = rainbow(length(unique(myclCTD)), start=0.1, end=0.9);
colorBrewerInterp=colorRampPalette(brewer.pal(8,"Accent"))

mycolhcCTD =  sample(colorBrewerInterp(length(unique(myclCTD))))
mycolhcCTD = mycolhcCTD[as.vector(myclCTD)]

t1=as.dendrogram(hr)
rowOrder=order.dendrogram(t1);



#Prepare clustered data for writing to table
yMatClustCTD=y[rowOrder,];
#geneLabelsClustCTD=geneLabelsCTD[colOrder];
clusterNumber=myclCTD[rowOrder];

#analytesCluster=data.frame(clusterNumber, yMatClustCTD)
#write.csv(analytesCluster, file="analytesCluster.csv", row.names=TRUE)


pdfPath=paste("./FiguresOut/Heatmap.pdf",sep='')
pdf(pdfPath, width=6,height=8,pointsize = 14)
heatmap3(yZ, ColSideColors =c(colorBarCombined), ColSideWidth=1, 
         col=barColors, breaks=breakBarColors, legendfun=function()showLegend(legend=c(NA),col=c(NA),cex=1.5), 
         Rowv=as.dendrogram(hr), Colv=NA,  scale="row",
         labRow = NA, labCol = NA,
         RowSideColors = mycolhcCTD)
dev.off()

--------
#Normalize RNAseq data
#https://lashlock.github.io/compbio/R_presentation.html
#scroll down to airway scaled counts paart
coldata=data.frame(sample=colnames(y), condition=pheno, row.names="sample")
coldata$condition=pheno
coldata$condition = factor(coldata$condition)
aby=y[which(!is.na(rownames(y))),] #have to delete row 14446 first col 327
dds = DESeqDataSetFromMatrix(countData = round(aby), colData = coldata, 
                             design = ~ condition)
dds = estimateSizeFactors(dds)
sizeFactors(dds) #gives the scale factor for each sample
yNorm = counts(dds, normalized=TRUE)
#dds$condition <- relevel(dds$condition, ref = "A")
dds = DESeq(dds)

#Plot the normalized data


#pdfPath=paste("./FiguresOut/nocluster.pdf",sep='')
#pdf(pdfPath, width=6,height=8,pointsize = 14)
heatmap3(yNorm, ColSideColors =colorBarCombined, ColSideWidth=1, 
         col=barColors, breaks=breakBarColors, legendfun=function()showLegend(legend=c(NA),col=c(NA),cex=1.5), 
         Rowv=NA, Colv=NA,  scale="row",
         labRow = NA, labCol = NA) 
#dev.off()


#Cluster and re-plot normalized data
yNormZ=t(apply(yNorm,1,scale))
rownames(yNormZ)=rownames(yNorm)

cort=cor(t(yNormZ), method="pearson")

hr = hclust(as.dist(1-cort), method="complete")

cutThreshCTD=1.47 #should be below max(hr$height)
myclCTD = cutree(hr, h=cutThreshCTD); 
#mycolhcCTD = rainbow(length(unique(myclCTD)), start=0.1, end=0.9);
colorBrewerInterp=colorRampPalette(brewer.pal(8,"Accent"))

mycolhcCTD =  sample(colorBrewerInterp(length(unique(myclCTD))))
mycolhcCTD = mycolhcCTD[as.vector(myclCTD)]

t1=as.dendrogram(hr)
rowOrder=order.dendrogram(t1);



#Prepare clustered data for writing to table
yMatClustCTD=y[rowOrder,];
#geneLabelsClustCTD=geneLabelsCTD[colOrder];
clusterNumber=myclCTD[rowOrder];

analytesCluster=data.frame(clusterNumber, yMatClustCTD)
write.csv(analytesCluster, file="analytesCluster.csv", row.names=TRUE)


pdfPath=paste("./FiguresOut/HeatmapDESeq.pdf",sep='')
pdf(pdfPath, width=6,height=8,pointsize = 14)
heatmap3(yNormZ, ColSideColors =colorBarCombined, ColSideWidth=1, 
         col=barColors, breaks=breakBarColors, ColSideLabs = NA, RowSideLabs = NA,
         Rowv=as.dendrogram(hr), Colv=NA,  scale="row",
         labRow = NA, labCol = NA,
         RowSideColors = mycolhcCTD)
dev.off()




# ---- GSVA ----

#Choose the data to run GSVA on.
yDataforGSVA=y;



#Choose whether or not to compute p values
compPermute=1

# Astrocyte gene list
astrocyteFXN1 = read.csv("GeneSets_Astrocytes_FXN1.csv")
astrocyteFXN1_list = as.list(as.data.frame(astrocyteFXN1))
astrocyteFXN1_set = lapply(astrocyteFXN1_list, function(x){ x[!is.na(x) & x != ""]})
astrocyteFXN1_set = na.omit(astrocyteFXN1_set)
astro_geneSetEnrich=gsva(yNorm, astrocyteFXN1_set, min.sz=5, mx.diff=TRUE, kcdf="Poisson")

astrocyteFXN2 = read.csv("GeneSets_Astrocytes_FXN2.csv")
astrocyteFXN2_list = as.list(as.data.frame(astrocyteFXN2))
astrocyteFXN2_set = lapply(astrocyteFXN2_list, function(x){ x[!is.na(x) & x != ""]})
astrocyteFXN2_set = na.omit(astrocyteFXN2_set)
astro2_geneSetEnrich=gsva(yNorm, astrocyteFXN2_set, min.sz=5, mx.diff=TRUE, kcdf="Poisson")

neuronFXN1 = read.csv("GeneSets_Neurons_FXN1.csv")
neuronFXN1_list = as.list(as.data.frame(neuronFXN1))
neuronFXN1_set = lapply(neuronFXN1_list, function(x){ x[!is.na(x) & x != ""]})
neuronFXN1_set = na.omit(neuronFXN1_set)
neuron_gse = gsva(yNorm, neuronFXN1_set, min.sz = 5, mx.diff = TRUE, kcdf = "Poisson")

#____C2 Gene Sets####
geneSetsMat=read.csv("c2.all.v7.4.symbols.csv")
geneSetsMat=geneSetsMat[,-2]
geneSetsMat_T=t(geneSetsMat)
colnames(geneSetsMat_T)=geneSetsMat_T[1,]
geneSetsMat_T=geneSetsMat_T[-1,]
geneSets3=as.list(as.data.frame(geneSetsMat_T))
#geneSets=lapply(geneSets3, function(x) x[!is.na(x)])
geneSets=lapply(geneSets3, function(x){ x[!is.na(x) & x != ""]})
geneSets=na.omit(geneSets)
#SP=split(geneSetsMat, col(geneSetsMat))

#geneSets2=list(Archana =geneSetsMat[,1], Levi=geneSetsMat[,2] )

#geneSets2=lapply(geneSets2, function(x) x[x=="#N/A"])

geneSetEnrich=gsva(round(yNorm), astrocyteFXN2_set, min.sz>1, mx.diff=TRUE, kcdf="Poisson")

#write.csv(geneSetEnrich, file="Gene_Scores_C7.csv", col.names=TRUE, row.names=TRUE, append=FALSE)

#----gene analysis
AvsB = results(dds, contrast = c("condition", "A", "B"))
AvsB = AvsB[order(AvsB$padj),]
write.csv(subset(AvsB, padj < 0.05 & log2FoldChange > 0), file="AvsBup.csv", row.names=TRUE)
write.csv(subset(AvsB, padj < 0.05 & log2FoldChange < 0), file="AvsBdown.csv", row.names=TRUE)
CvsB = results(dds, contrast = c("condition", "C", "B"))
CvsB = CvsB[order(CvsB$padj),]
write.csv(subset(CvsB, padj < 0.05 & log2FoldChange > 0), file="CvsBup.csv", row.names=TRUE)
write.csv(subset(CvsB, padj < 0.05 & log2FoldChange < 0), file="CvsBdown.csv", row.names=TRUE)



#----GSVA heatmap----


geneSetEnrichSortZ=t(apply(astro_geneSetEnrich[,],1,scale)); #z-score the data

cort=cor(t(geneSetEnrichSortZ), method="pearson")
hrGSVA = hclust(as.dist(1-cort), method="complete")

pdfPath=paste("./FiguresOut/GSVA_Cluster.pdf",sep='')
pdf(pdfPath, width=6,height=8,pointsize = 14)
heatmap3(geneSetEnrichSortZ, ColSideColors =colorBarCombined, ColSideWidth=1, ColSideLabs=NA,
         col=barColors, breaks=breakBarColors, legendfun=function()showLegend(legend=c(NA),col=c(NA)), 
         Rowv=as.dendrogram(hrGSVA), Colv=NA,  scale="none",
         labCol = NA, cexRow=1.2) 
title("GSVA for MDSIG gene sets", adj = 0.5, line = 2)

dev.off()
