library(DESeq2)
library(tidyr)
library(dplyr)
library(data.table)
library(ggplot2)
library(vsn)
library(pheatmap)
library(RColorBrewer)

setwd("/Users/yunqingyu/Dropbox/postdoc/data/sh1 mutant/RNAseq/DEseq2_Hao")

####### input table into DESeq2 ######
StarMatrix <- as.matrix(read.csv("A10sh1_RNAseq_new.star.counts.output.colsorted.csv", row.names=1))
colnames(StarMatrix) <- gsub("(.*)_[A|G|C|T]{6}_trim", "\\1", colnames(StarMatrix))
StarMatrix_AZ <- StarMatrix[,c(1:28, 49:76)]
colnames(StarMatrix_AZ)


coldata <- read.csv("A10sh1_RNAseq.design.csv",header=T)
coldata$X <- gsub("(.*)_[A|G|C|T]{6}_trim", "\\1", coldata$X)
coldata_AZ <- coldata[coldata$X %in% colnames(StarMatrix_AZ),]
coldata_AZ$X == colnames(StarMatrix_AZ)
coldata_AZ <- coldata_AZ[,c("genotype", "stage", "tissue", "rep",	"group")]

dds_A10sh10_AZ <- DESeqDataSetFromMatrix(countData = StarMatrix_AZ, colData = coldata_AZ, 
                                          design = ~ group)

# Only keep genes with at least 3 samples with 10 reads or more
keep <- rowSums(counts(dds_A10sh10_AZ) >= 10) >= 3
dds_A10sh1_AZ <- dds_A10sh10_AZ[keep,]
# DE anlysis and get read counts
dds_A10sh1_AZ <- DESeq(dds_A10sh1_AZ, minReplicatesForReplace=Inf)
plotDispEsts(dds_A10sh1_AZ)
counts_A10sh1 <- as.data.frame(counts(dds_A10sh1_AZ, normalized = TRUE))
dim(counts_A10sh1)
colnames(counts_A10sh1)
counts_A10sh1$Sv_geneID <- rownames(counts_A10sh1)
write.table(counts_A10sh1, "A10sh1_counts.txt", sep = "\t")

# boxplot of the Cook’s distances
par(mar=c(4,3,1,3))
boxplot(log10(assays(dds_A10sh1_AZ)[["cooks"]]), range=0, las=2)



#extracting transformed values
vsd <- vst(dds_A10sh1_AZ, blind = FALSE)

head(assay(vsd))
meanSdPlot(assay(vsd), ranks = FALSE)
A10sh1_vsd <- as.data.frame(assay(vsd))
dim(A10sh1_vsd)
A10sh1_vsd$Sv_geneID <- rownames(A10sh1_vsd)
write.table(A10sh1_vsd, "A10sh1_vst.txt", sep = "\t")

#PCA graph
pcaData <- plotPCA(vsd, intgroup = c("genotype", "stage", "tissue"), returnData = TRUE)
pcaData
percentVar <- round(100*attr(pcaData, "percentVar"))
pdf("A10sh1_PCA.pdf", width = 6, height = 5)
ggplot(pcaData, aes(x = PC1, y = PC2, alpha = genotype)) + geom_point(aes(shape = tissue, fill = stage, color = stage)) +
  scale_alpha_discrete(range = c(0.5, 1)) +
  scale_shape_manual(values = c(21, 22, 23, 24)) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) + 
  ylab(paste0("PC2: ", percentVar[2], "% variance"))
dev.off()


### Get all data without filter
### sh1 vs. A10
res_sh1_21d_A.vs.A10_21d_A <- results(dds_A10sh1_AZ, 
                                      contrast = c("group", "sh1_21d_A", "A10_21d_A"))
summary(res_sh1_21d_A.vs.A10_21d_A)
dim(res_sh1_21d_A.vs.A10_21d_A)

res_sh1_31d_A.vs.A10_31d_A <- results(dds_A10sh1_AZ, 
                                      contrast = c("group", "sh1_31d_A", "A10_31d_A"))
summary(res_sh1_31d_A.vs.A10_31d_A)
dim(res_sh1_31d_A.vs.A10_31d_A)

res_sh1_38d_A.vs.A10_38d_A <- results(dds_A10sh1_AZ, 
                                      contrast = c("group", "sh1_38d_A", "A10_38d_A"))
summary(res_sh1_38d_A.vs.A10_38d_A)
dim(res_sh1_38d_A.vs.A10_38d_A)

res_sh1_21d_L.vs.A10_21d_L <- results(dds_A10sh1_AZ, 
                                      contrast = c("group", "sh1_21d_L", "A10_21d_L"))
summary(res_sh1_21d_L.vs.A10_21d_L)
dim(res_sh1_21d_L.vs.A10_21d_L)

res_sh1_31d_L.vs.A10_31d_L <- results(dds_A10sh1_AZ, 
                                      contrast = c("group", "sh1_31d_L", "A10_31d_L"))
summary(res_sh1_31d_L.vs.A10_31d_L)
dim(res_sh1_31d_L.vs.A10_31d_L)

res_sh1_21d_U.vs.A10_21d_U <- results(dds_A10sh1_AZ, 
                                      contrast = c("group", "sh1_21d_U", "A10_21d_U"))
summary(res_sh1_21d_U.vs.A10_21d_U)
dim(res_sh1_21d_U.vs.A10_21d_U)

res_sh1_31d_U.vs.A10_31d_U <- results(dds_A10sh1_AZ, 
                                      contrast = c("group", "sh1_31d_U", "A10_31d_U"))
summary(res_sh1_31d_U.vs.A10_31d_U)
dim(res_sh1_31d_U.vs.A10_31d_U)


### A v L / A v U in sh1
res_sh1_21d_A.vs.sh1_21d_L <- results(dds_A10sh1_AZ, 
                                      contrast = c("group", "sh1_21d_A", "sh1_21d_L"))
summary(res_sh1_21d_A.vs.sh1_21d_L)
dim(res_sh1_21d_A.vs.sh1_21d_L)

res_sh1_21d_A.vs.sh1_21d_U <- results(dds_A10sh1_AZ, 
                                      contrast = c("group", "sh1_21d_A", "sh1_21d_U"))
summary(res_sh1_21d_A.vs.sh1_21d_U)
dim(res_sh1_21d_A.vs.sh1_21d_U)

res_sh1_31d_A.vs.sh1_31d_L <- results(dds_A10sh1_AZ, 
                                      contrast = c("group", "sh1_31d_A", "sh1_31d_L"))
summary(res_sh1_31d_A.vs.sh1_31d_L)
dim(res_sh1_31d_A.vs.sh1_31d_L)

res_sh1_31d_A.vs.sh1_31d_U <- results(dds_A10sh1_AZ, 
                                      contrast = c("group", "sh1_31d_A", "sh1_31d_U"))
summary(res_sh1_31d_A.vs.sh1_31d_U)
dim(res_sh1_31d_A.vs.sh1_31d_U)

### A v L / A v U in A10
res_A10_21d_A.vs.A10_21d_L <- results(dds_A10sh1_AZ, 
                                      contrast = c("group", "A10_21d_A", "A10_21d_L"))
summary(res_A10_21d_A.vs.A10_21d_L)
dim(res_A10_21d_A.vs.A10_21d_L)

res_A10_21d_A.vs.A10_21d_U <- results(dds_A10sh1_AZ, 
                                      contrast = c("group", "A10_21d_A", "A10_21d_U"))
summary(res_A10_21d_A.vs.A10_21d_U)
dim(res_A10_21d_A.vs.A10_21d_U)


res_A10_31d_A.vs.A10_31d_L <- results(dds_A10sh1_AZ, 
                                      contrast = c("group", "A10_31d_A", "A10_31d_L"))
summary(res_A10_31d_A.vs.A10_31d_L)
dim(res_A10_31d_A.vs.A10_31d_L)

res_A10_31d_A.vs.A10_31d_U <- results(dds_A10sh1_AZ, 
                                      contrast = c("group", "A10_31d_A", "A10_31d_U"))
summary(res_A10_31d_A.vs.A10_31d_U)
dim(res_A10_31d_A.vs.A10_31d_U)

### 31/21d in sh1
res_sh1_31d_A.vs.sh1_21d_A <- results(dds_A10sh1_AZ, 
                                      contrast = c("group", "sh1_31d_A", "sh1_21d_A"))
summary(res_sh1_31d_A.vs.sh1_21d_A)
dim(res_sh1_31d_A.vs.sh1_21d_A)


res_sh1_31d_L.vs.sh1_21d_L <- results(dds_A10sh1_AZ, 
                                      contrast = c("group", "sh1_31d_L", "sh1_21d_L"))
summary(res_sh1_31d_L.vs.sh1_21d_L)
dim(res_sh1_31d_L.vs.sh1_21d_L)

res_sh1_31d_U.vs.sh1_21d_U <- results(dds_A10sh1_AZ, 
                                      contrast = c("group", "sh1_31d_U", "sh1_21d_U"))
summary(res_sh1_31d_U.vs.sh1_21d_U)
dim(res_sh1_31d_U.vs.sh1_21d_U)

### 31/21d in A10
res_A10_31d_A.vs.A10_21d_A <- results(dds_A10sh1_AZ, 
                                      contrast = c("group", "A10_31d_A", "A10_21d_A"))
summary(res_A10_31d_A.vs.A10_21d_A)
dim(res_A10_31d_A.vs.A10_21d_A)


res_A10_31d_L.vs.A10_21d_L <- results(dds_A10sh1_AZ, 
                                      contrast = c("group", "A10_31d_L", "A10_21d_L"))
summary(res_A10_31d_L.vs.A10_21d_L)
dim(res_A10_31d_L.vs.A10_21d_L)

res_A10_31d_U.vs.A10_21d_U <- results(dds_A10sh1_AZ, 
                                      contrast = c("group", "A10_31d_U", "A10_21d_U"))
summary(res_A10_31d_U.vs.A10_21d_U)
dim(res_A10_31d_U.vs.A10_21d_U)


### 38/31d A in sh1 and A10
res_sh1_38d_A.vs.sh1_31d_A <- results(dds_A10sh1_AZ, 
                                      contrast = c("group", "sh1_38d_A", "sh1_31d_A"))
summary(res_sh1_38d_A.vs.sh1_31d_A)
dim(res_sh1_38d_A.vs.sh1_31d_A)

res_A10_38d_A.vs.A10_31d_A <- results(dds_A10sh1_AZ, 
                                      contrast = c("group", "A10_38d_A", "A10_31d_A"))
summary(res_A10_38d_A.vs.A10_31d_A)
dim(res_A10_38d_A.vs.A10_31d_A)






####join the unfiltered res_tables together
##join the unfiltered DE table and gene count tables together
# Mutant effect
names(res_sh1_21d_A.vs.A10_21d_A) <-paste0(names(res_sh1_21d_A.vs.A10_21d_A),"_sh1_21d_A.vs.A10_21d_A")
names(res_sh1_31d_A.vs.A10_31d_A) <-paste0(names(res_sh1_31d_A.vs.A10_31d_A),"_sh1_31d_A.vs.A10_31d_A")
names(res_sh1_38d_A.vs.A10_38d_A) <-paste0(names(res_sh1_38d_A.vs.A10_38d_A),"_sh1_38d_A.vs.A10_38d_A")

names(res_sh1_21d_L.vs.A10_21d_L) <-paste0(names(res_sh1_21d_L.vs.A10_21d_L),"_sh1_21d_L.vs.A10_21d_L")
names(res_sh1_31d_L.vs.A10_31d_L) <-paste0(names(res_sh1_31d_L.vs.A10_31d_L),"_sh1_31d_L.vs.A10_31d_L")

names(res_sh1_21d_U.vs.A10_21d_U) <-paste0(names(res_sh1_21d_U.vs.A10_21d_U),"_sh1_21d_U.vs.A10_21d_U")
names(res_sh1_31d_U.vs.A10_31d_U) <-paste0(names(res_sh1_31d_U.vs.A10_31d_U),"_sh1_31d_U.vs.A10_31d_U")

# Tissue effect
names(res_sh1_21d_A.vs.sh1_21d_L) <-paste0(names(res_sh1_21d_A.vs.sh1_21d_L),"_sh1_21d_A.vs.sh1_21d_L")
names(res_sh1_21d_A.vs.sh1_21d_U) <-paste0(names(res_sh1_21d_A.vs.sh1_21d_U),"_sh1_21d_A.vs.sh1_21d_U")
names(res_sh1_31d_A.vs.sh1_31d_L) <-paste0(names(res_sh1_31d_A.vs.sh1_31d_L),"_sh1_31d_A.vs.sh1_31d_L")
names(res_sh1_31d_A.vs.sh1_31d_U) <-paste0(names(res_sh1_31d_A.vs.sh1_31d_U),"_sh1_31d_A.vs.sh1_31d_U")

names(res_A10_21d_A.vs.A10_21d_L) <-paste0(names(res_A10_21d_A.vs.A10_21d_L),"_A10_21d_A.vs.A10_21d_L")
names(res_A10_21d_A.vs.A10_21d_U) <-paste0(names(res_A10_21d_A.vs.A10_21d_U),"_A10_21d_A.vs.A10_21d_U")
names(res_A10_31d_A.vs.A10_31d_L) <-paste0(names(res_A10_31d_A.vs.A10_31d_L),"_A10_31d_A.vs.A10_31d_L")
names(res_A10_31d_A.vs.A10_31d_U) <-paste0(names(res_A10_31d_A.vs.A10_31d_U),"_A10_31d_A.vs.A10_31d_U")

# Stage effect
names(res_sh1_38d_A.vs.sh1_31d_A) <-paste0(names(res_sh1_38d_A.vs.sh1_31d_A),"_sh1_38d_A.vs.sh1_31d_A")

names(res_sh1_31d_A.vs.sh1_21d_A) <-paste0(names(res_sh1_31d_A.vs.sh1_21d_A),"_sh1_31d_A.vs.sh1_21d_A")
names(res_sh1_31d_L.vs.sh1_21d_L) <-paste0(names(res_sh1_31d_L.vs.sh1_21d_L),"_sh1_31d_L.vs.sh1_21d_L")
names(res_sh1_31d_U.vs.sh1_21d_U) <-paste0(names(res_sh1_31d_U.vs.sh1_21d_U),"_sh1_31d_U.vs.sh1_21d_U")

names(res_A10_31d_A.vs.A10_21d_A) <-paste0(names(res_A10_31d_A.vs.A10_21d_A),"_A10_31d_A.vs.A10_21d_A")
names(res_A10_31d_L.vs.A10_21d_L) <-paste0(names(res_A10_31d_L.vs.A10_21d_L),"_A10_31d_L.vs.A10_21d_L")
names(res_A10_31d_U.vs.A10_21d_U) <-paste0(names(res_A10_31d_U.vs.A10_21d_U),"_A10_31d_U.vs.A10_21d_U")


res_A10sh1 <- cbind(res_sh1_21d_A.vs.A10_21d_A, res_sh1_31d_A.vs.A10_31d_A, res_sh1_38d_A.vs.A10_38d_A,
                    res_sh1_21d_L.vs.A10_21d_L, res_sh1_31d_L.vs.A10_31d_L, 
                    res_sh1_21d_U.vs.A10_21d_U, res_sh1_31d_U.vs.A10_31d_U, 

                    res_sh1_21d_A.vs.sh1_21d_L, res_sh1_21d_A.vs.sh1_21d_U, 
                    res_sh1_31d_A.vs.sh1_31d_L, res_sh1_31d_A.vs.sh1_31d_U, 
                    res_A10_21d_A.vs.A10_21d_L, res_A10_21d_A.vs.A10_21d_U, 
                    res_A10_31d_A.vs.A10_31d_L, res_A10_31d_A.vs.A10_31d_U, 
                    
                    res_sh1_38d_A.vs.sh1_31d_A, res_sh1_31d_A.vs.sh1_21d_A, res_sh1_31d_L.vs.sh1_21d_L, res_sh1_31d_U.vs.sh1_21d_U, 
                    res_A10_38d_A.vs.A10_31d_A, res_A10_31d_A.vs.A10_21d_A, res_A10_31d_L.vs.A10_21d_L, res_A10_31d_U.vs.A10_21d_U, 
                    
                    counts_A10sh1)




res_A10sh1 <- as.data.frame(res_A10sh1) %>% mutate(A10_21d_L = rowMeans(select(., A10_21d_L_rep1:A10_21d_L_rep4)),
                                                   A10_21d_A = rowMeans(select(., A10_21d_A_rep1:A10_21d_A_rep4)),
                                                   A10_21d_U = rowMeans(select(., A10_21d_U_rep1:A10_21d_U_rep4)),
                                                   A10_31d_L = rowMeans(select(., A10_31d_L_rep1:A10_31d_L_rep4)),
                                                   A10_31d_A = rowMeans(select(., A10_31d_A_rep1:A10_31d_A_rep4)),
                                                   A10_31d_U = rowMeans(select(., A10_31d_U_rep1:A10_31d_U_rep4)),
                                                   A10_38d_A = rowMeans(select(., A10_38d_A_rep1:A10_38d_A_rep4)),                                         sh1_21d_L = rowMeans(select(., sh1_21d_L_rep1:sh1_21d_L_rep4)),
                                                   sh1_21d_A = rowMeans(select(., sh1_21d_A_rep1:sh1_21d_A_rep4)),
                                                   sh1_21d_U = rowMeans(select(., sh1_21d_U_rep1:sh1_21d_U_rep4)),
                                                   sh1_31d_L = rowMeans(select(., sh1_31d_L_rep1:sh1_31d_L_rep4)),
                                                   sh1_31d_A = rowMeans(select(., sh1_31d_A_rep1:sh1_31d_A_rep4)),
                                                   sh1_31d_U = rowMeans(select(., sh1_31d_U_rep1:sh1_31d_U_rep4)),
                                                   sh1_38d_A = rowMeans(select(., sh1_38d_A_rep1:sh1_38d_A_rep4)))
dim(res_A10sh1)
names(res_A10sh1)
##input Sv annotation
Sv_anno <- fread("/Users/yunqingyu/Dropbox/postdoc/data/LCM/rna seq/reference genome/PhytozomeV12/SviridisV2/annotation/Sviridis_500_v2.1.annotation_info.txt", 
                 header = T, sep = "\t", select = c(1:16))%>%tbl_df()
dim(Sv_anno)
names(Sv_anno)
Sv_anno <-Sv_anno %>%
  dplyr::select(geneID = locusName, Best_hit_arabi_name:rice_defline, Pfam:GO)
Sv_anno <-Sv_anno[!duplicated(Sv_anno$geneID),]
dim(Sv_anno)
names(Sv_anno) <- paste0("Sv_", names(Sv_anno))
head(Sv_anno)

##join the result table and annotation together
res_A10sh1_anno <- left_join(res_A10sh1, Sv_anno, by = "Sv_geneID")
dim(res_A10sh1_anno)
write.table(res_A10sh1_anno, "res_A10sh1_AZ_counts_anno_DEseq2.txt", row.names = F, sep = "\t")








library(limma)
#Create a dataframe first, include all pairwise compasrion, and set to 0
df_A10sh1_FDR0.05_FC1.5 <- data.frame(Sv_geneID = counts_A10sh1$Sv_geneID, 
                                      # Mutant effect AZ
                                      sh1_21d_A.vs.A10_21d_A = 0, 
                                      sh1_31d_A.vs.A10_31d_A = 0, 
                                      sh1_38d_A.vs.A10_38d_A = 0, 
                                      sh1_21d_L.vs.A10_21d_L = 0, 
                                      sh1_31d_L.vs.A10_31d_L = 0, 
                                      sh1_21d_U.vs.A10_21d_U = 0,
                                      sh1_31d_U.vs.A10_31d_U = 0,
                                      
                                      # tissue effect
                                      sh1_21d_A.vs.sh1_21d_L = 0, 
                                      sh1_21d_A.vs.sh1_21d_U = 0, 
                                      sh1_31d_A.vs.sh1_31d_L = 0, 
                                      sh1_31d_A.vs.sh1_31d_U = 0,
                                      A10_21d_A.vs.A10_21d_L = 0, 
                                      A10_21d_A.vs.A10_21d_U = 0, 
                                      A10_31d_A.vs.A10_31d_L = 0,
                                      A10_31d_A.vs.A10_31d_U = 0)


#get the up and down regulated genes
# Extract desired vaules from res_data set
# Get Mutant responses

FDR0.05_FC1.5_sh1_21d_A.vs.A10_21d_A_up <- rownames(subset(res_sh1_21d_A.vs.A10_21d_A, 
                                                           padj_sh1_21d_A.vs.A10_21d_A < 0.05 & log2FoldChange_sh1_21d_A.vs.A10_21d_A >= log2(1.5)))
FDR0.05_FC1.5_sh1_21d_A.vs.A10_21d_A_down <- rownames(subset(res_sh1_21d_A.vs.A10_21d_A, 
                                                             padj_sh1_21d_A.vs.A10_21d_A < 0.05 & log2FoldChange_sh1_21d_A.vs.A10_21d_A <= -log2(1.5)))

FDR0.05_FC1.5_sh1_31d_A.vs.A10_31d_A_up <- rownames(subset(res_sh1_31d_A.vs.A10_31d_A, 
                                                           padj_sh1_31d_A.vs.A10_31d_A < 0.05 & log2FoldChange_sh1_31d_A.vs.A10_31d_A >= log2(1.5)))
FDR0.05_FC1.5_sh1_31d_A.vs.A10_31d_A_down <- rownames(subset(res_sh1_31d_A.vs.A10_31d_A, 
                                                             padj_sh1_31d_A.vs.A10_31d_A < 0.05 & log2FoldChange_sh1_31d_A.vs.A10_31d_A <= -log2(1.5)))

FDR0.05_FC1.5_sh1_38d_A.vs.A10_38d_A_up <- rownames(subset(res_sh1_38d_A.vs.A10_38d_A, 
                                                           padj_sh1_38d_A.vs.A10_38d_A < 0.05 & log2FoldChange_sh1_38d_A.vs.A10_38d_A >= log2(1.5)))
FDR0.05_FC1.5_sh1_38d_A.vs.A10_38d_A_down <- rownames(subset(res_sh1_38d_A.vs.A10_38d_A, 
                                                             padj_sh1_38d_A.vs.A10_38d_A < 0.05 & log2FoldChange_sh1_38d_A.vs.A10_38d_A <= -log2(1.5)))

FDR0.05_FC1.5_sh1_21d_L.vs.A10_21d_L_up <- rownames(subset(res_sh1_21d_L.vs.A10_21d_L, 
                                                           padj_sh1_21d_L.vs.A10_21d_L < 0.05 & log2FoldChange_sh1_21d_L.vs.A10_21d_L >= log2(1.5)))
FDR0.05_FC1.5_sh1_21d_L.vs.A10_21d_L_down <- rownames(subset(res_sh1_21d_L.vs.A10_21d_L, 
                                                             padj_sh1_21d_L.vs.A10_21d_L < 0.05 & log2FoldChange_sh1_21d_L.vs.A10_21d_L <= -log2(1.5)))

FDR0.05_FC1.5_sh1_31d_L.vs.A10_31d_L_up <- rownames(subset(res_sh1_31d_L.vs.A10_31d_L, 
                                                           padj_sh1_31d_L.vs.A10_31d_L < 0.05 & log2FoldChange_sh1_31d_L.vs.A10_31d_L >= log2(1.5)))
FDR0.05_FC1.5_sh1_31d_L.vs.A10_31d_L_down <- rownames(subset(res_sh1_31d_L.vs.A10_31d_L, 
                                                             padj_sh1_31d_L.vs.A10_31d_L < 0.05 & log2FoldChange_sh1_31d_L.vs.A10_31d_L <= -log2(1.5)))

FDR0.05_FC1.5_sh1_21d_U.vs.A10_21d_U_up <- rownames(subset(res_sh1_21d_U.vs.A10_21d_U, 
                                                           padj_sh1_21d_U.vs.A10_21d_U < 0.05 & log2FoldChange_sh1_21d_U.vs.A10_21d_U >= log2(1.5)))
FDR0.05_FC1.5_sh1_21d_U.vs.A10_21d_U_down <- rownames(subset(res_sh1_21d_U.vs.A10_21d_U, 
                                                             padj_sh1_21d_U.vs.A10_21d_U < 0.05 & log2FoldChange_sh1_21d_U.vs.A10_21d_U <= -log2(1.5)))

FDR0.05_FC1.5_sh1_31d_U.vs.A10_31d_U_up <- rownames(subset(res_sh1_31d_U.vs.A10_31d_U, 
                                                           padj_sh1_31d_U.vs.A10_31d_U < 0.05 & log2FoldChange_sh1_31d_U.vs.A10_31d_U >= log2(1.5)))
FDR0.05_FC1.5_sh1_31d_U.vs.A10_31d_U_down <- rownames(subset(res_sh1_31d_U.vs.A10_31d_U, 
                                                             padj_sh1_31d_U.vs.A10_31d_U < 0.05 & log2FoldChange_sh1_31d_U.vs.A10_31d_U <= -log2(1.5)))


# tissue effect
FDR0.05_FC1.5_sh1_21d_A.vs.sh1_21d_L_up <- rownames(subset(res_sh1_21d_A.vs.sh1_21d_L, 
                                                           padj_sh1_21d_A.vs.sh1_21d_L < 0.05 & log2FoldChange_sh1_21d_A.vs.sh1_21d_L >= log2(1.5)))
FDR0.05_FC1.5_sh1_21d_A.vs.sh1_21d_L_down <- rownames(subset(res_sh1_21d_A.vs.sh1_21d_L, 
                                                             padj_sh1_21d_A.vs.sh1_21d_L < 0.05 & log2FoldChange_sh1_21d_A.vs.sh1_21d_L <= -log2(1.5)))

FDR0.05_FC1.5_sh1_21d_A.vs.sh1_21d_U_up <- rownames(subset(res_sh1_21d_A.vs.sh1_21d_U, 
                                                           padj_sh1_21d_A.vs.sh1_21d_U < 0.05 & log2FoldChange_sh1_21d_A.vs.sh1_21d_U >= log2(1.5)))
FDR0.05_FC1.5_sh1_21d_A.vs.sh1_21d_U_down <- rownames(subset(res_sh1_21d_A.vs.sh1_21d_U, 
                                                             padj_sh1_21d_A.vs.sh1_21d_U < 0.05 & log2FoldChange_sh1_21d_A.vs.sh1_21d_U <= -log2(1.5)))

FDR0.05_FC1.5_sh1_31d_A.vs.sh1_31d_L_up <- rownames(subset(res_sh1_31d_A.vs.sh1_31d_L, 
                                                           padj_sh1_31d_A.vs.sh1_31d_L < 0.05 & log2FoldChange_sh1_31d_A.vs.sh1_31d_L >= log2(1.5)))
FDR0.05_FC1.5_sh1_31d_A.vs.sh1_31d_L_down <- rownames(subset(res_sh1_31d_A.vs.sh1_31d_L, 
                                                             padj_sh1_31d_A.vs.sh1_31d_L < 0.05 & log2FoldChange_sh1_31d_A.vs.sh1_31d_L <= -log2(1.5)))

FDR0.05_FC1.5_sh1_31d_A.vs.sh1_31d_U_up <- rownames(subset(res_sh1_31d_A.vs.sh1_31d_U, 
                                                           padj_sh1_31d_A.vs.sh1_31d_U < 0.05 & log2FoldChange_sh1_31d_A.vs.sh1_31d_U >= log2(1.5)))
FDR0.05_FC1.5_sh1_31d_A.vs.sh1_31d_U_down <- rownames(subset(res_sh1_31d_A.vs.sh1_31d_U, 
                                                             padj_sh1_31d_A.vs.sh1_31d_U < 0.05 & log2FoldChange_sh1_31d_A.vs.sh1_31d_U <= -log2(1.5)))

FDR0.05_FC1.5_A10_21d_A.vs.A10_21d_L_up <- rownames(subset(res_A10_21d_A.vs.A10_21d_L, 
                                                           padj_A10_21d_A.vs.A10_21d_L < 0.05 & log2FoldChange_A10_21d_A.vs.A10_21d_L >= log2(1.5)))
FDR0.05_FC1.5_A10_21d_A.vs.A10_21d_L_down <- rownames(subset(res_A10_21d_A.vs.A10_21d_L, 
                                                             padj_A10_21d_A.vs.A10_21d_L < 0.05 & log2FoldChange_A10_21d_A.vs.A10_21d_L <= -log2(1.5)))

FDR0.05_FC1.5_A10_21d_A.vs.A10_21d_U_up <- rownames(subset(res_A10_21d_A.vs.A10_21d_U, 
                                                           padj_A10_21d_A.vs.A10_21d_U < 0.05 & log2FoldChange_A10_21d_A.vs.A10_21d_U >= log2(1.5)))
FDR0.05_FC1.5_A10_21d_A.vs.A10_21d_U_down <- rownames(subset(res_A10_21d_A.vs.A10_21d_U, 
                                                             padj_A10_21d_A.vs.A10_21d_U < 0.05 & log2FoldChange_A10_21d_A.vs.A10_21d_U <= -log2(1.5)))

FDR0.05_FC1.5_A10_31d_A.vs.A10_31d_L_up <- rownames(subset(res_A10_31d_A.vs.A10_31d_L, 
                                                           padj_A10_31d_A.vs.A10_31d_L < 0.05 & log2FoldChange_A10_31d_A.vs.A10_31d_L >= log2(1.5)))
FDR0.05_FC1.5_A10_31d_A.vs.A10_31d_L_down <- rownames(subset(res_A10_31d_A.vs.A10_31d_L, 
                                                             padj_A10_31d_A.vs.A10_31d_L < 0.05 & log2FoldChange_A10_31d_A.vs.A10_31d_L <= -log2(1.5)))

FDR0.05_FC1.5_A10_31d_A.vs.A10_31d_U_up <- rownames(subset(res_A10_31d_A.vs.A10_31d_U, 
                                                           padj_A10_31d_A.vs.A10_31d_U < 0.05 & log2FoldChange_A10_31d_A.vs.A10_31d_U >= log2(1.5)))
FDR0.05_FC1.5_A10_31d_A.vs.A10_31d_U_down <- rownames(subset(res_A10_31d_A.vs.A10_31d_U, 
                                                             padj_A10_31d_A.vs.A10_31d_U < 0.05 & log2FoldChange_A10_31d_A.vs.A10_31d_U <= -log2(1.5)))


#assign upregulated genes as 1 and downregulated genes as -1
# For some unknown reason, "-1" must be firstly assigned before "1"
# Downregulated genes
# Mutant effect
df_A10sh1_FDR0.05_FC1.5$sh1_21d_A.vs.A10_21d_A[df_A10sh1_FDR0.05_FC1.5$Sv_geneID%in%FDR0.05_FC1.5_sh1_21d_A.vs.A10_21d_A_down] <- -1
df_A10sh1_FDR0.05_FC1.5$sh1_31d_A.vs.A10_31d_A[df_A10sh1_FDR0.05_FC1.5$Sv_geneID%in%FDR0.05_FC1.5_sh1_31d_A.vs.A10_31d_A_down] <- -1
df_A10sh1_FDR0.05_FC1.5$sh1_38d_A.vs.A10_38d_A[df_A10sh1_FDR0.05_FC1.5$Sv_geneID%in%FDR0.05_FC1.5_sh1_38d_A.vs.A10_38d_A_down] <- -1
df_A10sh1_FDR0.05_FC1.5$sh1_21d_L.vs.A10_21d_L[df_A10sh1_FDR0.05_FC1.5$Sv_geneID%in%FDR0.05_FC1.5_sh1_21d_L.vs.A10_21d_L_down] <- -1
df_A10sh1_FDR0.05_FC1.5$sh1_31d_L.vs.A10_31d_L[df_A10sh1_FDR0.05_FC1.5$Sv_geneID%in%FDR0.05_FC1.5_sh1_31d_L.vs.A10_31d_L_down] <- -1
df_A10sh1_FDR0.05_FC1.5$sh1_21d_U.vs.A10_21d_U[df_A10sh1_FDR0.05_FC1.5$Sv_geneID%in%FDR0.05_FC1.5_sh1_21d_U.vs.A10_21d_U_down] <- -1
df_A10sh1_FDR0.05_FC1.5$sh1_31d_U.vs.A10_31d_U[df_A10sh1_FDR0.05_FC1.5$Sv_geneID%in%FDR0.05_FC1.5_sh1_31d_U.vs.A10_31d_U_down] <- -1

# Tissue effect
df_A10sh1_FDR0.05_FC1.5$sh1_21d_A.vs.sh1_21d_L[df_A10sh1_FDR0.05_FC1.5$Sv_geneID%in%FDR0.05_FC1.5_sh1_21d_A.vs.sh1_21d_L_down] <- -1
df_A10sh1_FDR0.05_FC1.5$sh1_21d_A.vs.sh1_21d_U[df_A10sh1_FDR0.05_FC1.5$Sv_geneID%in%FDR0.05_FC1.5_sh1_21d_A.vs.sh1_21d_U_down] <- -1
df_A10sh1_FDR0.05_FC1.5$sh1_31d_A.vs.sh1_31d_L[df_A10sh1_FDR0.05_FC1.5$Sv_geneID%in%FDR0.05_FC1.5_sh1_31d_A.vs.sh1_31d_L_down] <- -1
df_A10sh1_FDR0.05_FC1.5$sh1_31d_A.vs.sh1_31d_U[df_A10sh1_FDR0.05_FC1.5$Sv_geneID%in%FDR0.05_FC1.5_sh1_31d_A.vs.sh1_31d_U_down] <- -1
df_A10sh1_FDR0.05_FC1.5$A10_21d_A.vs.A10_21d_L[df_A10sh1_FDR0.05_FC1.5$Sv_geneID%in%FDR0.05_FC1.5_A10_21d_A.vs.A10_21d_L_down] <- -1
df_A10sh1_FDR0.05_FC1.5$A10_21d_A.vs.A10_21d_U[df_A10sh1_FDR0.05_FC1.5$Sv_geneID%in%FDR0.05_FC1.5_A10_21d_A.vs.A10_21d_U_down] <- -1
df_A10sh1_FDR0.05_FC1.5$A10_31d_A.vs.A10_31d_L[df_A10sh1_FDR0.05_FC1.5$Sv_geneID%in%FDR0.05_FC1.5_A10_31d_A.vs.A10_31d_L_down] <- -1
df_A10sh1_FDR0.05_FC1.5$A10_31d_A.vs.A10_31d_U[df_A10sh1_FDR0.05_FC1.5$Sv_geneID%in%FDR0.05_FC1.5_A10_31d_A.vs.A10_31d_U_down] <- -1

# upregulated genes
# Mutant effect
df_A10sh1_FDR0.05_FC1.5$sh1_21d_A.vs.A10_21d_A[df_A10sh1_FDR0.05_FC1.5$Sv_geneID%in%FDR0.05_FC1.5_sh1_21d_A.vs.A10_21d_A_up] <- 1
df_A10sh1_FDR0.05_FC1.5$sh1_31d_A.vs.A10_31d_A[df_A10sh1_FDR0.05_FC1.5$Sv_geneID%in%FDR0.05_FC1.5_sh1_31d_A.vs.A10_31d_A_up] <- 1
df_A10sh1_FDR0.05_FC1.5$sh1_38d_A.vs.A10_38d_A[df_A10sh1_FDR0.05_FC1.5$Sv_geneID%in%FDR0.05_FC1.5_sh1_38d_A.vs.A10_38d_A_up] <- 1
df_A10sh1_FDR0.05_FC1.5$sh1_21d_L.vs.A10_21d_L[df_A10sh1_FDR0.05_FC1.5$Sv_geneID%in%FDR0.05_FC1.5_sh1_21d_L.vs.A10_21d_L_up] <- 1
df_A10sh1_FDR0.05_FC1.5$sh1_31d_L.vs.A10_31d_L[df_A10sh1_FDR0.05_FC1.5$Sv_geneID%in%FDR0.05_FC1.5_sh1_31d_L.vs.A10_31d_L_up] <- 1
df_A10sh1_FDR0.05_FC1.5$sh1_21d_U.vs.A10_21d_U[df_A10sh1_FDR0.05_FC1.5$Sv_geneID%in%FDR0.05_FC1.5_sh1_21d_U.vs.A10_21d_U_up] <- 1
df_A10sh1_FDR0.05_FC1.5$sh1_31d_U.vs.A10_31d_U[df_A10sh1_FDR0.05_FC1.5$Sv_geneID%in%FDR0.05_FC1.5_sh1_31d_U.vs.A10_31d_U_up] <- 1

# Tissue effect
df_A10sh1_FDR0.05_FC1.5$sh1_21d_A.vs.sh1_21d_L[df_A10sh1_FDR0.05_FC1.5$Sv_geneID%in%FDR0.05_FC1.5_sh1_21d_A.vs.sh1_21d_L_up] <- 1
df_A10sh1_FDR0.05_FC1.5$sh1_21d_A.vs.sh1_21d_U[df_A10sh1_FDR0.05_FC1.5$Sv_geneID%in%FDR0.05_FC1.5_sh1_21d_A.vs.sh1_21d_U_up] <- 1
df_A10sh1_FDR0.05_FC1.5$sh1_31d_A.vs.sh1_31d_L[df_A10sh1_FDR0.05_FC1.5$Sv_geneID%in%FDR0.05_FC1.5_sh1_31d_A.vs.sh1_31d_L_up] <- 1
df_A10sh1_FDR0.05_FC1.5$sh1_31d_A.vs.sh1_31d_U[df_A10sh1_FDR0.05_FC1.5$Sv_geneID%in%FDR0.05_FC1.5_sh1_31d_A.vs.sh1_31d_U_up] <- 1
df_A10sh1_FDR0.05_FC1.5$A10_21d_A.vs.A10_21d_L[df_A10sh1_FDR0.05_FC1.5$Sv_geneID%in%FDR0.05_FC1.5_A10_21d_A.vs.A10_21d_L_up] <- 1
df_A10sh1_FDR0.05_FC1.5$A10_21d_A.vs.A10_21d_U[df_A10sh1_FDR0.05_FC1.5$Sv_geneID%in%FDR0.05_FC1.5_A10_21d_A.vs.A10_21d_U_up] <- 1
df_A10sh1_FDR0.05_FC1.5$A10_31d_A.vs.A10_31d_L[df_A10sh1_FDR0.05_FC1.5$Sv_geneID%in%FDR0.05_FC1.5_A10_31d_A.vs.A10_31d_L_up] <- 1
df_A10sh1_FDR0.05_FC1.5$A10_31d_A.vs.A10_31d_U[df_A10sh1_FDR0.05_FC1.5$Sv_geneID%in%FDR0.05_FC1.5_A10_31d_A.vs.A10_31d_U_up] <- 1

pdf("venn_sh1vsA10.pdf", width = 5, height = 8)
par(mfrow = c(3,1), mar = c(0, 0, 1, 0))
vennDiagram(df_A10sh1_FDR0.05_FC1.5[,2:4], include = c("up","down"), mar=rep(0, 4), 
            circle.col = c("red", "blue", "green"), cex = c(1, 2, 3))
title(main = "DEG_sh1_A.vs.A10_A (FC>=1.5, FDR<0.05)", cex.main = 1)

vennDiagram(df_A10sh1_FDR0.05_FC1.5[,5:6], include = c("up","down"), mar=rep(0, 4), 
            circle.col = c("red", "blue"), cex = c(1, 2, 3))
title(main = "DEG_sh1_L.vs.A10_L (FC>=1.5, FDR<0.05)", cex.main = 1)

vennDiagram(df_A10sh1_FDR0.05_FC1.5[,7:8], include = c("up","down"), mar=rep(0, 4), 
            circle.col = c("red", "blue"), cex = c(1, 2, 3))
title(main = "DEG_sh1_U.vs.A10_U (FC>=1.5, FDR<0.05)", cex.main = 1)
dev.off()



pdf("venn_AvsLU.pdf", width = 12, height = 5)
par(mfrow = c(2,2), mar = c(0, 0, 1, 0))
vennDiagram(df_A10sh1_FDR0.05_FC1.5[,14:15], include = c("up","down"), mar=rep(0.1, 4), 
            circle.col = c("red", "blue"), cex = c(1, 2, 3))
title(main = "DEG-21d-sh1_A/L.vs.A/U (FC>=1.5, FDR<0.05)", cex.main = 1)

vennDiagram(df_A10sh1_FDR0.05_FC1.5[,16:17], include = c("up","down"), mar=rep(0.1, 4), 
            circle.col = c("red", "blue"), cex = c(1, 2, 3))
title(main = "DEG-31d-sh1_A/L.vs.A/U (FC>=1.5, FDR<0.05)", cex.main = 1)

vennDiagram(df_A10sh1_FDR0.05_FC1.5[,18:19], include = c("up","down"), mar=rep(0.1, 4), 
            circle.col = c("red", "blue"), cex = c(1, 2, 3))
title(main = "DEG-21d-A10_A/L.vs.A/U (FC>=1.5, FDR<0.05)", cex.main = 1)

vennDiagram(df_A10sh1_FDR0.05_FC1.5[,20:21], include = c("up","down"), mar=rep(0.1, 4), 
            circle.col = c("red", "blue"), cex = c(1, 2, 3))
title(main = "DEG-31d-A10_A/L.vs.A/U (FC>=1.5, FDR<0.05)", cex.main = 1)
dev.off()



### Make a count table of up and down regulated genes
DEG_counttable <- data.frame(sh1vsA10 = c("up", "down"),
                             L_21d = 0, A_21d = 0, U_21d = 0, L_31d = 0, A_31d = 0, U_31d = 0, A_38d = 0)
DEG_counttable[1,2:8] <- c(length(FDR0.05_FC1.5_sh1_21d_L.vs.A10_21d_L_up),
                           length(FDR0.05_FC1.5_sh1_21d_A.vs.A10_21d_A_up),
                           length(FDR0.05_FC1.5_sh1_21d_U.vs.A10_21d_U_up),
                           length(FDR0.05_FC1.5_sh1_31d_L.vs.A10_31d_L_up),
                           length(FDR0.05_FC1.5_sh1_31d_A.vs.A10_31d_A_up),
                           length(FDR0.05_FC1.5_sh1_31d_U.vs.A10_31d_U_up),
                           length(FDR0.05_FC1.5_sh1_38d_A.vs.A10_38d_A_up))

DEG_counttable[2,2:8] <- c(length(FDR0.05_FC1.5_sh1_21d_L.vs.A10_21d_L_down),
                           length(FDR0.05_FC1.5_sh1_21d_A.vs.A10_21d_A_down),
                           length(FDR0.05_FC1.5_sh1_21d_U.vs.A10_21d_U_down),
                           length(FDR0.05_FC1.5_sh1_31d_L.vs.A10_31d_L_down),
                           length(FDR0.05_FC1.5_sh1_31d_A.vs.A10_31d_A_down),
                           length(FDR0.05_FC1.5_sh1_31d_U.vs.A10_31d_U_down),
                           length(FDR0.05_FC1.5_sh1_38d_A.vs.A10_38d_A_down))

DEG_counttable
write.table(DEG_counttable, "DEG_count_table.txt", row.names = F, sep = "\t", quote = F)


### matrix of normalized counts 
DEG_shvsA10 <- unique(c(FDR0.05_FC1.5_sh1_21d_L.vs.A10_21d_L_down,
                        FDR0.05_FC1.5_sh1_21d_A.vs.A10_21d_A_down,
                        FDR0.05_FC1.5_sh1_21d_U.vs.A10_21d_U_down,
                        FDR0.05_FC1.5_sh1_31d_L.vs.A10_31d_L_down,
                        FDR0.05_FC1.5_sh1_31d_A.vs.A10_31d_A_down,
                        FDR0.05_FC1.5_sh1_31d_U.vs.A10_31d_U_down,
                        FDR0.05_FC1.5_sh1_38d_A.vs.A10_38d_A_down,
                        
                        FDR0.05_FC1.5_sh1_21d_L.vs.A10_21d_L_up,
                        FDR0.05_FC1.5_sh1_21d_A.vs.A10_21d_A_up,
                        FDR0.05_FC1.5_sh1_21d_U.vs.A10_21d_U_up,
                        FDR0.05_FC1.5_sh1_31d_L.vs.A10_31d_L_up,
                        FDR0.05_FC1.5_sh1_31d_A.vs.A10_31d_A_up,
                        FDR0.05_FC1.5_sh1_31d_U.vs.A10_31d_U_up,
                        FDR0.05_FC1.5_sh1_38d_A.vs.A10_38d_A_up))
length(DEG_shvsA10)

counts_DEGsh1vsA10 <- counts_A10sh1[DEG_shvsA10,]
dim(counts_DEGsh1vsA10)
write.table(counts_DEGsh1vsA10, "DEG_sh1vsA10_counts.txt", sep = "\t")


res_DEGsh1vsA10 <- res[res$Sv_geneID %in% DEG_sh1vsA10_counts$Sv_geneID,]
dim(res_DEGsh1vsA10)
write_xlsx(res_DEGsh1vsA10, "res_DEGsh1vsA10.xlsx")




