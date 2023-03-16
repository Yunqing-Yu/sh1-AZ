#Install SvOrgDb from local for the 1st run 
install.packages("/Users/yunqingyu/Dropbox/postdoc/data/sh1 mutant/RNAseq/DEseq2_Hao/clusterProfiler_A10sh1/Sviridis_500_v2.1_GOMAP_DB/org.Sviridis.eg.db",repos = NULL, type="source")

#Load clusterProfiler for GO analysis
library(clusterProfiler)
library(org.Sviridis.eg.db)
library(GOSemSim)
library(enrichplot)
library(writexl)

######### Prepare input for GO analysis #########
# Load gene list from GO_input
files <- list.files("/Users/yunqingyu/Dropbox/postdoc/data/sh1 mutant/RNAseq/DEseq2_Hao/SOM/SSOM3x3", pattern = ".txt")
files
gene_list <- list()
for (i in files){
  f <- read.table(paste0("/Users/yunqingyu/Dropbox/postdoc/data/sh1 mutant/RNAseq/DEseq2_Hao/SOM/SSOM3x3/", i), header = F)[,1]
  module <- gsub("(.*)A10sh1.txt", "\\1", i)
  gene_list[[module]] <- f
}
names(gene_list)
lapply(gene_list, length)

# Make GOsemsimDATA object used for enrichment analysis
SemData_A10sh1_BP <- godata(OrgDb = "org.Sviridis.eg.db", keytype = "GID", ont = "BP", computeIC = TRUE)
SemData_A10sh1_CC <- godata(OrgDb = "org.Sviridis.eg.db", keytype = "GID", ont = "CC", computeIC = TRUE)
SemData_A10sh1_MF <- godata(OrgDb = "org.Sviridis.eg.db", keytype = "GID", ont = "MF", computeIC = TRUE)



######################## Analyze different modules together (this is much easier than the individual ones, and essentially the same as doing separately) ####################3#####3
### Use compareCluster for comparing between gene lists
## BP
SSOM_BP <- compareCluster(gene_list, fun = "enrichGO",
                          OrgDb = "org.Sviridis.eg.db", keyType="GID", ont="BP", 
                          pAdjustMethod = "BH", qvalueCutoff = 0.05, 
                          minGSSize = 10, maxGSSize = 1000)
names(SSOM_BP@compareClusterResult)
GO_dup_BP <- SSOM_BP@compareClusterResult$ID[duplicated(SSOM_BP@compareClusterResult$geneID)]
length(GO_dup_BP)

# Drop the redundant GO with long text GO:0019288
SSOM_BP1 <- dropGO(SSOM_BP, term = "GO:0019288")

#Reduce GO redundancy with Rel
SSOM_BP_Rel <- simplify(SSOM_BP1, cutoff=0.7, by="p.adjust", measure = "Rel", semData = SemData_A10sh1_BP)
GO_dup_BP_Rel <- SSOM_BP_Rel@compareClusterResult$ID[duplicated(SSOM_BP_Rel@compareClusterResult$geneID)]
length(GO_dup_BP_Rel)

SSOM_BP_Rel_unique <- dropGO(SSOM_BP_Rel, term = GO_dup_BP_Rel)

pdf("dotplot_BP_Rel0.7_top10.pdf", width = 12, height = 15)
dotplot(SSOM_BP_Rel_unique, showCategory = 10, font.size =15, label_format = 60)
dev.off()



## CC
SSOM_CC <- compareCluster(gene_list, fun = "enrichGO",
                          OrgDb = "org.Sviridis.eg.db", keyType="GID", ont="CC", 
                          pAdjustMethod = "BH", qvalueCutoff = 0.05, 
                          minGSSize = 10, maxGSSize = 1000)
names(SSOM_CC@compareClusterResult)
GO_dup_CC <- SSOM_CC@compareClusterResult$ID[duplicated(SSOM_CC@compareClusterResult$geneID)]
length(GO_dup_CC)

#Reduce GO redundancy with Rel
SSOM_CC_Rel <- simplify(SSOM_CC, cutoff=0.7, by="p.adjust", measure = "Rel", semData = SemData_A10sh1_CC)
GO_dup_CC_Rel <- SSOM_CC_Rel@compareClusterResult$ID[duplicated(SSOM_CC_Rel@compareClusterResult$geneID)]
length(GO_dup_CC_Rel)

SSOM_CC_Rel_unique <- dropGO(SSOM_CC_Rel, term = GO_dup_CC_Rel)

pdf("dotplot/dotplot_CC_Rel0.7_top10.pdf", width = 9, height = 5)
dotplot(SSOM_CC_Rel_unique, showCategory = 10, font.size =15, label_format = 60)
dev.off()



## MF
SSOM_MF <- compareCluster(gene_list, fun = "enrichGO",
                          OrgDb = "org.Sviridis.eg.db", keyType="GID", ont="MF", 
                          pAdjustMethod = "BH", qvalueCutoff = 0.05, 
                          minGSSize = 10, maxGSSize = 1000)
names(SSOM_MF@compareClusterResult)
GO_dup_MF <- SSOM_MF@compareClusterResult$ID[duplicated(SSOM_MF@compareClusterResult$geneID)]
length(GO_dup_MF)

#Reduce GO redundancy with Rel
SSOM_MF_Rel <- simplify(SSOM_MF, cutoff=0.7, by="p.adjust", measure = "Rel", semData = SemData_A10sh1_MF)
GO_dup_MF_Rel <- SSOM_MF_Rel@compareClusterResult$ID[duplicated(SSOM_MF_Rel@compareClusterResult$geneID)]
length(GO_dup_MF_Rel)

SSOM_MF_Rel_unique <- dropGO(SSOM_MF_Rel, term = GO_dup_MF_Rel)

pdf("dotplot/dotplot_MF_Rel0.7_top10.pdf", width = 11, height = 5)
dotplot(SSOM_MF_Rel_unique, showCategory = 10, font.size =15, label_format = 60)
dev.off()


#### Focus on cluster 4-9
gene_list_49 <- gene_list[4:9]
  
  
## BP
SSOM_BP_49 <- compareCluster(gene_list_49, fun = "enrichGO",
                            OrgDb = "org.Sviridis.eg.db", keyType="GID", ont="BP", 
                            pAdjustMethod = "BH", qvalueCutoff = 0.05, 
                            minGSSize = 10, maxGSSize = 1000)
names(SSOM_BP_49@compareClusterResult)
GO_dup_BP_49 <- SSOM_BP_49@compareClusterResult$ID[duplicated(SSOM_BP_49@compareClusterResult$geneID)]
length(GO_dup_BP_49)

# Drop the redundant GO with long text GO:0019288
SSOM_BP1_49 <- dropGO(SSOM_BP_49, term = "GO:0019288")

#Reduce GO redundancy with Rel
SSOM_BP_49_Rel <- simplify(SSOM_BP1_49, cutoff=0.7, by="p.adjust", measure = "Rel", semData = SemData_A10sh1_BP)
GO_dup_BP_49_Rel <- SSOM_BP_49_Rel@compareClusterResult$ID[duplicated(SSOM_BP_49_Rel@compareClusterResult$geneID)]
length(GO_dup_BP_49_Rel)

SSOM_BP_49_Rel_unique <- dropGO(SSOM_BP_49_Rel, term = GO_dup_BP_49_Rel)

pdf("dotplot_SSOM3x3/dotplot_BP_49_Rel0.7_top10.pdf", width = 9, height = 13)
dotplot(SSOM_BP_49_Rel_unique, showCategory = 10, font.size =15, label_format = 60)
dev.off()

