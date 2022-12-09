# ---- Introduction ----
# This script analyzes RNAseq count data by normalizing and z-scoring the data
# with the intention of rendering a heatmap. Normalized, z-scored data is
# used to identify differentially expressed genes (DEGs), which are
# then processed in a frequency-specific, directionality-specific manner.
#
# Jacob Kraus, Singer Lab
# Fall 2022 PURA

# ---- Initialize Libraries ----
options(java.parameters = "-Xmx8g")
library("heatmap3")
library("RColorBrewer")
library("EnhancedVolcano")
library("gplots")
library("DESeq2")

# ---- Set heatmap3 color bar parameters ---- 
breakBarColors=c(-200,seq(-1.5, 1.5, 0.01),200) # Color range for z-scoring
barColors = colorpanel(length(breakBarColors)-1, "snow4", "white", "mediumorchid1")

# ---- Read and plot RNAseq data ----
dataIn=read.csv("RNAseq data 21days_raw counts.csv", header = TRUE)

# ---- Filter and format RNAseq data ----
yList=dataIn[, 2:(ncol(dataIn)-1)]

y=matrix(unlist(yList), ncol=ncol(yList), byrow=F)
rownames(y)=toupper(dataIn$geneName)
colnames(y)=colnames(dataIn)[2:(ncol(dataIn)-1)]

uniqueGenes=unique(rownames(y)) 

repG=match(uniqueGenes,rownames(y))

y=y[repG,]

indKeepGenes=which(rowSums(y>=4)>=5)
y=y[indKeepGenes,]

# ---- Setup phenotype column color bar ----
pheno=c(rep("A",3),rep("B",3),rep("C",3))

colorBar=matrix(,nrow=length(pheno), ncol=1)

colorBar[pheno == "A"]= "black"
colorBar[pheno == "B"]= "slateblue4"
colorBar[pheno == "C"]= "orangered4"

colorBarCombined=c(colorBar)

# ---- Normalize RNAseq data using DESeq2 ----
coldata=data.frame(sample=colnames(y), condition=pheno, row.names="sample")
coldata$condition=pheno
coldata$condition = factor(coldata$condition)
remo=y[which(!is.na(rownames(y))),]
dds = DESeqDataSetFromMatrix(countData = round(remo), colData = coldata, 
                             design = ~ condition)
dds = estimateSizeFactors(dds)
yNorm = counts(dds, normalized=TRUE)
dds = DESeq(dds)

# ---- Z-score normalized data and create dendrogram ----
yNormZ=t(apply(yNorm,1,scale))
rownames(yNormZ)=rownames(yNorm)

cort=cor(t(yNormZ), method="pearson")

hr = hclust(as.dist(1-cort), method="complete")

cutThreshCTD=1.47
myclCTD = cutree(hr, h=cutThreshCTD); 
colorBrewerInterp=colorRampPalette(brewer.pal(8,"Accent"))

mycolhcCTD =  sample(colorBrewerInterp(length(unique(myclCTD))))
mycolhcCTD = mycolhcCTD[as.vector(myclCTD)]

t1=as.dendrogram(hr)
rowOrder=order.dendrogram(t1);

yMatClustCTD=y[rowOrder,];
clusterNumber=myclCTD[rowOrder];

analytesCluster=data.frame(clusterNumber, yMatClustCTD)
write.csv(analytesCluster, file="analytesCluster.csv", row.names=TRUE)

# ---- Render heatmap of normalized, z-scored data ----
pdfPath=paste("./FiguresOut/HeatmapDESeq.pdf",sep='')
pdf(pdfPath, width=6,height=8,pointsize = 14)
heatmap3(yNormZ, ColSideColors =colorBarCombined, ColSideWidth=1, 
         col=barColors, breaks=breakBarColors, ColSideLabs = NA,
         RowSideLabs = NA, Rowv=as.dendrogram(hr), Colv=NA,  scale="row",
         labRow = NA, labCol = NA, RowSideColors = mycolhcCTD)
dev.off()


# ---- Identify and analyze 20 Hz and 40 Hz DEGs ----
AvsB = results(dds, contrast = c("condition", "A", "B"))
AvsB = AvsB[order(AvsB$padj),]
write.csv(subset(AvsB, padj < 0.05 & log2FoldChange > 0), file="AvsBup.csv",
          row.names=TRUE)
write.csv(subset(AvsB, padj < 0.05 & log2FoldChange < 0), file="AvsBdown.csv",
          row.names=TRUE)
EnhancedVolcano(AvsB, lab = rownames(AvsB), x = 'log2FoldChange', y = 'padj',
                title = '20 Hz versus Control', 
                pCutoff = 0.05, FCcutoff = 0, pointSize = 3.0, labSize = 6.0)
CvsB = results(dds, contrast = c("condition", "C", "B"))
CvsB = CvsB[order(CvsB$padj),]
write.csv(subset(CvsB, padj < 0.05 & log2FoldChange > 0), file="CvsBup.csv",
          row.names=TRUE)
write.csv(subset(CvsB, padj < 0.05 & log2FoldChange < 0), file="CvsBdown.csv",
          row.names=TRUE)
EnhancedVolcano(CvsB, lab = rownames(CvsB), x = 'log2FoldChange', y = 'padj',
                title = '40 Hz versus Control', 
                pCutoff = 0.05, FCcutoff = 0, pointSize = 3.0, labSize = 6.0)