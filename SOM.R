library(tidyr)
library(dplyr)
library(kohonen)
library(ggplot2)
library(writexl)

setwd("/Users/yunqingyu/Dropbox/postdoc/data/sh1 mutant/RNAseq/DEseq2_Hao/SOM")

## load count data from DEG (FDR0.05, FC1.5)
DEG_counts <- read.table("/Users/yunqingyu/Dropbox/postdoc/data/sh1 mutant/RNAseq/DEseq2_Hao/DEG_sh1vsA10_counts.txt", sep = "\t")
head(DEG_counts)

## Calculate mean of replicates
DEG_counts_mean0 <- DEG_counts %>% mutate(A10_21d_L = rowMeans(select(., A10_21d_L_rep1:A10_21d_L_rep4)),
                                          A10_21d_A = rowMeans(select(., A10_21d_A_rep1:A10_21d_A_rep4)),
                                          A10_21d_U = rowMeans(select(., A10_21d_U_rep1:A10_21d_U_rep4)),
                                          A10_31d_L = rowMeans(select(., A10_31d_L_rep1:A10_31d_L_rep4)),
                                          A10_31d_A = rowMeans(select(., A10_31d_A_rep1:A10_31d_A_rep4)),
                                          A10_31d_U = rowMeans(select(., A10_31d_U_rep1:A10_31d_U_rep4)),
                                          A10_38d_A = rowMeans(select(., A10_38d_A_rep1:A10_38d_A_rep4)),                                         
                                          sh1_21d_L = rowMeans(select(., sh1_21d_L_rep1:sh1_21d_L_rep4)),
                                          sh1_21d_A = rowMeans(select(., sh1_21d_A_rep1:sh1_21d_A_rep4)),
                                          sh1_21d_U = rowMeans(select(., sh1_21d_U_rep1:sh1_21d_U_rep4)),
                                          sh1_31d_L = rowMeans(select(., sh1_31d_L_rep1:sh1_31d_L_rep4)),
                                          sh1_31d_A = rowMeans(select(., sh1_31d_A_rep1:sh1_31d_A_rep4)),
                                          sh1_31d_U = rowMeans(select(., sh1_31d_U_rep1:sh1_31d_U_rep4)),
                                          sh1_38d_A = rowMeans(select(., sh1_38d_A_rep1:sh1_38d_A_rep4)))
names(DEG_counts_mean0)
DEG_counts_mean <- DEG_counts_mean0[,58:71]
head(DEG_counts_mean)

## Data needs to be scaled before SOM clustering
DEG_scale <- as.matrix(t(scale(t(DEG_counts_mean))))
head(DEG_scale)

## PCA map may show obvious clusters
pca <- prcomp(DEG_scale, scale = T)
summary(pca)
pca.scores <- data.frame(pca$x)

DEG_scale_pca <- cbind(DEG_scale, pca.scores)
head(DEG_scale_pca)

p <- ggplot(DEG_scale_pca, aes(PC1, PC2)) 
p + geom_point()

###### Run SOM #######

### super som: consider A10 and sh1 as different layers of data
set.seed(2)
DEG_scale_list <- list(DEG_scale[,1:7], DEG_scale[,8:14])
ssom <- supersom(data = DEG_scale_list, grid = somgrid(3, 3, "hexagonal"))
summary(ssom)

plot(ssom, type ="changes")
plot(ssom, type = "quality")
plot(ssom, type = "counts")
plot(ssom, type = "mapping")

pdf("ssom_fan_3x3.pdf", width = 6, height = 4)
par(mfrow = c(1,2), mar = c(0, 0, 0, 0))
plot(ssom, type = "codes")
par(mfrow = c(1,1))
dev.off()

DEG_scale_SSOMclass <- mutate(data.frame(DEG_scale), cluster = ssom$unit.classif, geneID = row.names(DEG_scale), geneNo = 0)
for (n in 1:max(ssom$unit.classif)){
  DEG_scale_SSOMclass$geneNo[DEG_scale_SSOMclass$cluster == n] <- sum(ssom$unit.classif == n)
}
DEG_scale_SSOMclass$cluster_geneNo <- paste0(DEG_scale_SSOMclass$cluster, " (", DEG_scale_SSOMclass$geneNo, ")")
head(DEG_scale_SSOMclass)

length(ssom$unit.classif)

DEG_scale_SSOMclass_long <- gather(DEG_scale_SSOMclass, key = "group", value = "scaledCounts", A10_21d_L:sh1_38d_A)%>%
  separate(group, into = c("genotype", "stage", "tissue"), sep = "_")%>%
  mutate(tissue_stage = paste0(tissue, "_", stage))%>%
  spread(key = genotype, value = scaledCounts)
head(DEG_scale_SSOMclass_long)
DEG_scale_SSOMclass_long$tissue_stage <- factor(DEG_scale_SSOMclass_long$tissue_stage,
                                                levels = c("L_21d", "A_21d", "U_21d", "L_31d", "A_31d", "U_31d", "A_38d"))

p <- ggplot(data = DEG_scale_SSOMclass_long, aes(x = tissue_stage)) +
  geom_line(aes(y = A10, group = geneID, color = "A10"), size = 0.2, alpha = 0.2) + 
  geom_line(aes(y = sh1, group = geneID, color = "sh1"), size = 0.2, alpha = 0.2) +  
  stat_summary(aes(y = A10, group = cluster_geneNo), geom="line", fun = "mean", size=0.8, color="blue", linetype="dashed") +
  stat_summary(aes(y = sh1, group = cluster_geneNo), geom="line", fun = "mean", size=0.8, color="red", linetype="dashed") +
  facet_wrap(.~cluster_geneNo) + 
  scale_color_manual(name = "genotype", values = c(A10 = "dodgerblue", sh1 = "salmon")) +
  ylab("scaled counts") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

p
pdf("A10sh1_SSOM_clusters_3x3.pdf", width = 7, height = 5)
p
dev.off()

## Output gene list from clusters
DEseq_res_anno <- read.table("/Users/yunqingyu/Dropbox/postdoc/data/sh1 mutant/RNAseq/DEseq2_Hao/res_A10sh1_counts_anno_DEseq2.txt", header = T, sep = "\t")
genelist_SSOM <- list()
reslist_SSOM <- list()
for (i in 1:max(ssom$unit.classif)){
  genelist_SSOM[[i]] <- DEG_scale_SSOMclass$geneID[DEG_scale_SSOMclass$cluster == i]
  write.table(genelist_SSOM[[i]], file = paste0("SSOM", i, "A10sh1.txt"), col.names = FALSE, row.names = FALSE, quote = FALSE, sep="\t")
  name <- paste0("SSOM",i)
  reslist_SSOM[[name]] <- DEseq_res_anno[DEseq_res_anno$Sv_geneID %in% genelist_SSOM[[i]],]
}
write_xlsx(reslist_SSOM, "res_A10sh1_SSOM_3x3.xlsx")  ## change file name when choosing different cluster size

### Look at some specific genes
DEG_scale_SSOMclass[DEG_scale_SSOMclass$geneID == "Sevir.9G153200",]  ## SH1 is in cluster 5
DEG_scale_SSOMclass[DEG_scale_SSOMclass$geneID == "Sevir.5G085400",]  ## LES1 is in cluster 7
