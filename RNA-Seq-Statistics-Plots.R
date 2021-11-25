# Plotting script by Simeon Hebrew and Zoe Alahouzou



#Loading all required packages

library(ggplot2)
library(scales)
library(lattice)
library(prettyR)
library(reshape)
library(plyr)
library(limma)
library(GGally)
library(RColorBrewer)
library(sva)
library(EDASeq)
library(edgeR)
library(directlabels)
library(ShortRead)
library(data.table)
library(splines)
library(VennDiagram)
library(statmod)
library(pheatmap)
library(ggrepel)
library(factoextra)
library(dplyr)
library(stringr)
library(gridExtra)

#*#*#*#*#*#*#* HISTOGRAMS #*#*#*#*#*#*
#*#*#*#**#*#*#*#*#*#*#*#*#*#*#**#*#*#*

#Set working directory
setwd("/Users/Mac/Desktop/Master-International/GSK3-Analysis/")

#Importing normalized counts
normalized.counts <- read.table("CrevsControl-Retina-normalized_counts.txt", header=TRUE, row.names=1)
head(normalized.counts)

#Plotting all histograms to confirm distribution profile 
tiff("Hist-General-normalized.counts.tiff",width = 1024, height = 1024, units = "px", pointsize = 32)
op <- par(mfrow=c(3, 2))  
lapply(seq(normalized.counts), function(x) 
  hist(x=normalized.counts[[x]], xlab=names(normalized.counts)[x],col = "green", main=paste("Histogram", names(normalized.counts)[x])))

dev.off()

#*#*#*#*#*#*#* PRINCIPAL COMPONENT ANALYSIS *#*#*#*#*#*#*##
#*#*#*#*#*#*#**#*#*#*#*#*#*#*#*#*#*#*#**#*#*#*#*#*#*#*#*#*#*#*
#*
#Importing required file
gene_length <-  read.table("../Master_Files_RNA-Seq/STAR-Ensembl-GRCm38.94-gene-length.txt", header=TRUE, row.names=1) 


#Importing raw counts, organizing column order and grouping data based on control vs mutant
raw_counts <- read.table("GSK3-alphaCre-Counts", header=TRUE, row.names=1)[,c(5:11)]
colnames(raw_counts)[2:7] <- c("Control1","Control2","Control3","Cre1","Cre2","Cre3")
group <- factor(c(rep("Control",3), rep("Cre",3)))
raw.DGEList <- DGEList(counts=raw_counts[,c(2:7)], group=group, genes = gene_length)

# Filtering Count-per-Million (CPM)
keep <- rowSums(cpm(raw.DGEList)>2) >=2
raw.DGEList <- raw.DGEList[keep, keep.lib.sizes=FALSE]
raw.DGEList$samples$lib.size <- colSums(raw.DGEList$counts)
raw.DGEList$samples

# Normalization of count data
TMM.DGEList <- calcNormFactors(raw.DGEList)
TMM.DGEList$samples

cpm.normalized <- round(cpm(TMM.DGEList, log=FALSE), digits = 1)

#Introducing transpose
counts.pca <- prcomp(t(log2(cpm.normalized + 0.01)), center=TRUE, scale=FALSE)

#Plotting PCA plot
tiff("Normalized-Counts-Cre-vs-Control_PCA2-3-v2.tiff",width = 960, height = 960, units = "px", pointsize = 16)
fviz_pca_ind(counts.pca, axes = c(2,3),
             col.ind = group, # Color by the quality of representation
             title = "Principal Component Analysis",
             subtitle = "Cre vs Control",
             palette = c("#1f78b4", "#6a3d9a"),
             legend.title = "Groups",
             repel = TRUE,    # Avoid text overlapping
             pointsize = 7, pointshape = 17, labelsize=6, ggtheme = theme_gray(),
             mean.point = FALSE, caption = "Source: PCA 2/3"
) +
  theme(text = element_text(size = 20),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 20))
dev.off()


#*#*#*#*#*#* VOLCANO PLOT *#*#*#*#*#*#*#
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

#Import Master file
Mstr <- read.delim("MSTR-ControlvsCre.txt", row.names = 1)
attach(Mstr)

#Set paramter values for plotting representation
FDR <- 0.05
FC <- 0.263
FC_Text.neg <- -1
FDR_text.neg <- 1e-5
FC_Text.pos <- 1
FDR_text.pos <- 1e-5

#Define significant and non-significant genes
Mstr$sig <- ifelse(FDR_CrevControl_Ret <= FDR & abs(Log2_FC_CrevControl_Ret) >= FC , "Sig", "Not Sig")

text <- ((Log2_FC_CrevControl_Ret <= FC_Text.neg & FDR_CrevControl_Ret <= FDR_text.neg) | (Log2_FC_CrevControl_Ret >= FC_Text.pos & FDR_CrevControl_Ret <= FDR_text.pos))

#Plotting Volcano Plot
Volcano_Plot <- ggplot(Mstr, aes(x = Log2_FC_CrevControl_Ret, y = -log10(FDR_CrevControl_Ret), size = 5)) +
  geom_point(aes(color = sig)) +
  scale_color_manual(values = c("green", "red")) +
  xlab("Log2 Fold Change") + ylab("-Log10 (FDR)")+
  theme_minimal() +
  theme_bw(base_size = 20) +
  scale_x_continuous(limits=c(-4,4), breaks=seq(-12,12,2))+
  theme(plot.title = element_text(lineheight=.8, face="bold", size=40))+
  theme(axis.title.x = element_text(face="bold", size=32),axis.text.x  = element_text(face="bold", vjust=0.5, size=25, colour="black")) +
  theme(axis.title.y = element_text(face="bold", size=32),axis.text.y  = element_text(face="bold", vjust=0.5, size=25, colour="black"))+
  geom_hline(yintercept=1.3,size=1, linetype = "solid", col="red" )+
  geom_vline(xintercept=-1:1,size=2, col="black", linetype = "dashed") +  
  geom_text_repel(
    data = subset(Mstr, text),
    aes(label = Gene.name),family = "Arial",
    size = 8,
    box.padding = unit(0.4, "lines"),
    point.padding = unit(0.6, "lines")
  )
tiff(filename="Volcano-CrevsControl-Retina.TIF", width=1024, height=1024, units="px", pointsize=14)
Volcano_Plot
dev.off()

#*#*#*#*#*#*#*#*#*HEATMAPS #*#*#*#*#*#*#*#*#
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

#Importing differentially expressed genes
data <- read.delim("DEG_Cre-Control-FC1.2-P0.05-TPM5.txt", header = TRUE)[,c(4,11,12,25:30)]

#Ordering the dataframe
data = select(data,Gene.name, Cre_1_TPM, Cre_2_TPM, Cre_3_TPM, Control_1_TPM, Control_2_TPM, Control_3_TPM, Mean_Brn3a_FPKM, Mean_Brn3b_FPKM)

data <- data[complete.cases(data),]
rownames(data) <- data[,1]
data[,1] <- NULL

#Using Z score to standardize by gene - we think it is more centred and normalized than using log
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}
zscore.data <- t(apply(data, 1, cal_z_score))
colnames(zscore.data)[1:8] <-c("Cre-1","Cre-2","Cre-3","Control-1","Control-2","Control-3", "Brn3a", "Brn3b")

#We plotted one way clustering and two way clustering

# Plotting one way clustering heatmap based on Z-score
jet <- brewer.pal(11, "Spectral")
tiff("2ways-HC-Cre-Control-DEG-FC2-F0.05-TPM5-z-score.tiff",width = 600, height = 5000, units = "px", pointsize = 10)
heatmap <- pheatmap(zscore.data, 
                    cluster_row = T, 
                    cluster_cols = F, 
                    color = rev(jet), 
                    fontsize = 15,
                    fontsize_row=12, 
                    legend = T,
                    gaps_col=c(3),
                    main = paste("Heatmap-Cre-Control-DEG-FC2-F0.05-TPM10-z.score"))
dev.off()

#Plotting two way clustering heatmap based on Z-score
jet <- brewer.pal(11, "Spectral")
tiff("2ways-HC-Cre-Control-DEG-FC2-F0.05-TPM5-z-score.tiff",width = 700, height = 5000, units = "px", pointsize = 10)
heatmap <- pheatmap(zscore.data, 
                    cluster_row = T, 
                    cluster_cols = T, 
                    color = rev(jet), 
                    fontsize = 15,
                    fontsize_row=12, 
                    legend = T,
                    gaps_col=c(3),
                    main = paste("2way-Heatmap-Cre-Control-DEG-FC2-F0.05-TPM10-z.score"))
dev.off()

#*#*#*#*#*#*#* VENN DIAGRAM #*#*#*#*#*#*
#*#*#*#*#*#*#*#*#**#*#*#*#*#*#*#*#*#*#**

#Import required files

DEG_Cre <- read.delim("DEG_unfilt_Cre-Control-FC1.2-F0.05-TPM5.txt", row.names = 1, header = TRUE)
Brn3ab <- read.delim("Expressed-genes-Brn3ab-FPKM51.2-P0.05-TPM5.txt", row.names = 1, header = TRUE)

#Create intersection table and unique tables for venn digram plotting
intersect <- subset(DEG_Cre, row.names(DEG_Cre) %in% row.names(Brn3ab))
write.table(data.frame("Ensembl_ID"=rownames(intersect),intersect),
            file = paste("Intersect-Deg-and-Brn3ab.txt", sep=""), quote=FALSE, sep = "\t", row.names = FALSE)


uniq_DEG_Cre <- subset(DEG_Cre, !(row.names(DEG_Cre) %in% row.names(Brn3ab)))
write.table(data.frame("Ensembl_ID"=rownames(uniq_DEG_Cre),uniq_DEG_Cre),
            file = paste("Uniq-DEG.txt", sep=""), quote=FALSE, sep = "\t", row.names = FALSE)


uniq_Brn3ab <-  subset(Brn3ab, !(row.names(Brn3ab) %in% row.names(DEG_Cre)))
write.table(data.frame("Ensembl_ID"=rownames(uniq_Brn3ab),uniq_Brn3ab),
            file = paste("Uniq-Brn3ab.txt", sep=""), quote=FALSE, sep = "\t", row.names = FALSE)

area1 <- nrow(DEG_Cre)
area2 <- nrow(Brn3ab)
cross <- nrow(intersect) 


#Plotting Venn Diagram
g <- draw.pairwise.venn(area1 = area1, area2 = area2, cross.area = cross, category = c("DEG_Cre", "Brn3ab"), 
                        scaled = FALSE, cex=2.5, lty = 3, lwd = 1.5, col = "#666666", 
                        cat.cex=3, cat.pos = c(0, 0), cat.dist = 0.05, cat.fontface = "bold",
                        fontfamily = ("Arial"), 
                        fill = c("yellow", "dodgerblue") ,  cat.fontfamily = rep("Arial", 2), 
                        print.mode = c("percent","raw"))
tiff(filename=paste("VD-DEG-Cre-and-Expressed.tiff", sep=""), width=1024, height=1024, units="px", pointsize=12) 
grid.arrange(gTree(children=g), top=textGrob("VennDiagram-DEG-Cre-and-Brn3ab",gp=gpar(fontsize=36,font=2)))
dev.off()


#Extracting gene names for pathway analysis
genes <-as.data.frame(intersect[,3])
genes
write.table(genes,"list_of_genes.txt",sep="\t",row.names=FALSE,col.names = FALSE)

