#Loading Libraries
library(BiocParallel)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(ggplot2)
library(dplyr)
library(readr)
library(DESeq2)
library(readxl)
library(tidyverse)
library(limma)
library(pheatmap)
library(GEOquery)
library(EnhancedVolcano) 
library(clusterProfiler)
library(enrichplot)

#Importing Upregulated genes 
gse52202_upGene_C9 <- read.csv("C:/Users/LENOVO/Desktop/Workflow for internship/Group A (iPSC derived motor neuron)/GSE52202/Upregulated genes C9")
gse210969_upGene_SOD1 <- read.csv("C:/Users/LENOVO/Desktop/Workflow for internship/Group A (iPSC derived motor neuron)/GSE210969/upregulated gene SOD1")
gse210969_upGene_TDP43 <- read.csv("C:/Users/LENOVO/Desktop/Workflow for internship/Group A (iPSC derived motor neuron)/GSE210969/upregulated gene TDP43")
gse283507_upGene_TDP43 <- read.csv("C:/Users/LENOVO/Desktop/Workflow for internship/Group A (iPSC derived motor neuron)/GSE283507/upregulated genes TDP43")

#Importing downregulated genes
gse52202_downGene_C9 <- read.csv("C:/Users/LENOVO/Desktop/Workflow for internship/Group A (iPSC derived motor neuron)/GSE52202/Downregulated genes C9")
gse210969_downGene_SOD1 <- read.csv("C:/Users/LENOVO/Desktop/Workflow for internship/Group A (iPSC derived motor neuron)/GSE210969/downregulated genes SOD1")
gse210969_downGene_TDP43 <- read.csv("C:/Users/LENOVO/Desktop/Workflow for internship/Group A (iPSC derived motor neuron)/GSE210969/downregulated genes TDP43")
gse283507_downGene_TDP43 <- read.csv("C:/Users/LENOVO/Desktop/Workflow for internship/Group A (iPSC derived motor neuron)/GSE283507/downregulated genes TDP43")



#Gene IDs to Gene symbol
# Function to add Gene Symbol column
add_gene_symbol <- function(df) {
  
  df$SYMBOL <- mapIds(
    org.Hs.eg.db,
    keys = as.character(df$X),   #convert integer to character
    column = "SYMBOL",
    keytype = "ENTREZID",
    multiVals = "first"
  )
  
  return(df)
}

# Upregulated genes
gse52202_upGene_C9      <- add_gene_symbol(gse52202_upGene_C9)
gse210969_upGene_SOD1   <- add_gene_symbol(gse210969_upGene_SOD1)
gse210969_upGene_TDP43  <- add_gene_symbol(gse210969_upGene_TDP43)
gse283507_upGene_TDP43  <- add_gene_symbol(gse283507_upGene_TDP43)

# Downregulated genes
gse52202_downGene_C9      <- add_gene_symbol(gse52202_downGene_C9)
gse210969_downGene_SOD1   <- add_gene_symbol(gse210969_downGene_SOD1)
gse210969_downGene_TDP43  <- add_gene_symbol(gse210969_downGene_TDP43)
gse283507_downGene_TDP43  <- add_gene_symbol(gse283507_downGene_TDP43)

colnames(gse52202_upGene_C9)

str(gse52202_upGene_C9)

head(gse52202_upGene_C9)


#Removing genes that couldn't be mappped 
#UPREGULATED
gse52202_upGene_C9      <- gse52202_upGene_C9[!is.na(gse52202_upGene_C9$SYMBOL), ]
gse210969_upGene_SOD1   <- gse210969_upGene_SOD1[!is.na(gse210969_upGene_SOD1$SYMBOL), ]
gse210969_upGene_TDP43  <- gse210969_upGene_TDP43[!is.na(gse210969_upGene_TDP43$SYMBOL), ]
gse283507_upGene_TDP43  <- gse283507_upGene_TDP43[!is.na(gse283507_upGene_TDP43$SYMBOL), ]


#DOWNREGULATED
gse52202_downGene_C9      <- gse52202_downGene_C9[!is.na(gse52202_downGene_C9$SYMBOL), ]
gse210969_downGene_SOD1   <- gse210969_downGene_SOD1[!is.na(gse210969_downGene_SOD1$SYMBOL), ]
gse210969_downGene_TDP43  <- gse210969_downGene_TDP43[!is.na(gse210969_downGene_TDP43$SYMBOL), ]
gse283507_downGene_TDP43  <- gse283507_downGene_TDP43[!is.na(gse283507_downGene_TDP43$SYMBOL), ]



#Creating Gene Vector

#UPREGULATED
geneList_up_02_C9      <- as.character(gse52202_upGene_C9$X)
geneList_up_69_SOD1    <- as.character(gse210969_upGene_SOD1$X)
geneList_up_69_TDP43   <- as.character(gse210969_upGene_TDP43$X)
geneList_up_07_TDP43   <- as.character(gse283507_upGene_TDP43$X)


#DOWNREGULATED
geneList_down_02_C9      <- as.character(gse52202_downGene_C9$X)
geneList_down_69_SOD1    <- as.character(gse210969_downGene_SOD1$X)
geneList_down_69_TDP43   <- as.character(gse210969_downGene_TDP43$X)
geneList_down_07_TDP43   <- as.character(gse283507_downGene_TDP43$X)

#Checks
length(geneList_up_02_C9)
length(geneList_down_02_C9)

length(geneList_up_69_SOD1)
length(geneList_down_69_SOD1)

length(geneList_up_69_TDP43)
length(geneList_down_69_TDP43)

length(geneList_up_07_TDP43)
length(geneList_down_07_TDP43)


head(geneList_up_02_C9)
head(geneList_down_02_C9)

ego_BP_up_69_SOD1@gene
length(ego_BP_up_69_SOD1@gene)
head(geneList_up_69_SOD1)

AnnotationDbi::select(
  org.Hs.eg.db,
  keys = geneList_up_69_SOD1[1:10],
  columns = "SYMBOL",
  keytype = "ENTREZID"
)

ego_BP_up_69_SOD1
summary(ego_BP_up_69_SOD1)

ego_test <- enrichGO(
  gene = geneList_up_69_SOD1,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 1,
  qvalueCutoff = 1,
  readable = TRUE
)

nrow(as.data.frame(ego_test))


# GO annotations

#UPREGULATED GENES 
# GSE52202 C9
ego_BP_up_02_C9 <- enrichGO(
  gene = geneList_up_02_C9,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

go_BP_up_02_C9 <- as.data.frame(ego_BP_up_02_C9)
nrow(go_BP_up_02_C9)

write.csv(go_BP_up_02_C9,
          "GO_BP_GSE52202_C9_UP.csv",
          row.names = FALSE)



# GSE210969 SOD1
ego_BP_up_69_SOD1 <- enrichGO(
  gene = geneList_up_69_SOD1,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

go_BP_up_69_SOD1 <- as.data.frame(ego_BP_up_69_SOD1)
nrow(go_BP_up_69_SOD1)

write.csv(go_BP_up_69_SOD1,
          "GO_BP_GSE210969_SOD1_UP.csv",
          row.names = FALSE)



# GSE210969 TDP43
ego_BP_up_69_TDP43 <- enrichGO(
  gene = geneList_up_69_TDP43,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

go_BP_up_69_TDP43 <- as.data.frame(ego_BP_up_69_TDP43)
nrow(go_BP_up_69_TDP43)

write.csv(go_BP_up_69_TDP43,
          "GO_BP_GSE210969_TDP43_UP.csv",
          row.names = FALSE)



# GSE283507 TDP43
ego_BP_up_07_TDP43 <- enrichGO(
  gene = geneList_up_07_TDP43,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

go_BP_up_07_TDP43 <- as.data.frame(ego_BP_up_07_TDP43)

write.csv(go_BP_up_07_TDP43,
          "GO_BP_GSE283507_TDP43_UP.csv",
          row.names = FALSE)



#DOWNREGULATED GENES

# GSE52202 C9
ego_BP_down_02_C9 <- enrichGO(
  gene = geneList_down_02_C9,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

go_BP_down_02_C9 <- as.data.frame(ego_BP_down_02_C9)
nrow(go_BP_down_02_C9)

write.csv(go_BP_down_02_C9,
          "GO_BP_GSE52202_C9_DOWN.csv",
          row.names = FALSE)



# GSE210969 SOD1
ego_BP_down_69_SOD1 <- enrichGO(
  gene = geneList_down_69_SOD1,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

go_BP_down_69_SOD1 <- as.data.frame(ego_BP_down_69_SOD1)

write.csv(go_BP_down_69_SOD1,
          "GO_BP_GSE210969_SOD1_DOWN.csv",
          row.names = FALSE)



# GSE210969 TDP43
ego_BP_down_69_TDP43 <- enrichGO(
  gene = geneList_down_69_TDP43,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

go_BP_down_69_TDP43 <- as.data.frame(ego_BP_down_69_TDP43)

write.csv(go_BP_down_69_TDP43,
          "GO_BP_GSE210969_TDP43_DOWN.csv",
          row.names = FALSE)



# GSE283507 TDP43
ego_BP_down_07_TDP43 <- enrichGO(
  gene = geneList_down_07_TDP43,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

go_BP_down_07_TDP43 <- as.data.frame(ego_BP_down_07_TDP43)

write.csv(go_BP_down_07_TDP43,
          "GO_BP_GSE283507_TDP43_DOWN.csv",
          row.names = FALSE)

#Plots

# DOT PLOTS
# GSE52202 C9 (DOWN)
dotplot(
  ego_BP_down_02_C9,
  showCategory = 20,
  font.size = 12,
  title = "GO Biological Process - GSE52202 C9 (DOWN)"
)

# GSE210969 SOD1 (DOWN)
dotplot(
  ego_BP_down_69_SOD1,
  showCategory = 20,
  font.size = 12,
  title = "GO Biological Process - GSE210969 SOD1 (DOWN)"
)

# GSE210969 TDP43 (UP)
dotplot(
  ego_BP_up_69_TDP43,
  showCategory = 20,
  font.size = 12,
  title = "GO Biological Process - GSE210969 TDP43 (UP)"
)

# GSE283507 TDP43 (UP)
dotplot(
  ego_BP_up_07_TDP43,
  showCategory = 20,
  font.size = 12,
  title = "GO Biological Process - GSE283507 TDP43 (UP)"
)

# GSE283507 TDP43 (DOWN)
dotplot(
  ego_BP_down_07_TDP43,
  showCategory = 20,
  font.size = 12,
  title = "GO Biological Process - GSE283507 TDP43 (DOWN)"
)

#BAR PLOTS
# GSE52202 C9 (DOWN)
barplot(
  ego_BP_down_02_C9,
  showCategory = 20,
  title = "GO Biological Process - GSE52202 C9 (DOWN)"
)

# GSE210969 SOD1 (DOWN)
barplot(
  ego_BP_down_69_SOD1,
  showCategory = 20,
  title = "GO Biological Process - GSE210969 SOD1 (DOWN)"
)

# GSE210969 TDP43 (UP)
barplot(
  ego_BP_up_69_TDP43,
  showCategory = 20,
  title = "GO Biological Process - GSE210969 TDP43 (UP)"
)

# GSE283507 TDP43 (UP)
barplot(
  ego_BP_up_07_TDP43,
  showCategory = 20,
  title = "GO Biological Process - GSE283507 TDP43 (UP)"
)

# GSE283507 TDP43 (DOWN)
barplot(
  ego_BP_down_07_TDP43,
  showCategory = 20,
  title = "GO Biological Process - GSE283507 TDP43 (DOWN)"
)


#Tern similarity
ego_BP_down_02_C9   <- pairwise_termsim(ego_BP_down_02_C9)

ego_BP_down_69_SOD1 <- pairwise_termsim(ego_BP_down_69_SOD1)

ego_BP_up_69_TDP43  <- pairwise_termsim(ego_BP_up_69_TDP43)

ego_BP_up_07_TDP43  <- pairwise_termsim(ego_BP_up_07_TDP43)

ego_BP_down_07_TDP43 <- pairwise_termsim(ego_BP_down_07_TDP43)




#ENRICHMENT MAPS

emapplot(
  ego_BP_down_02_C9,
  showCategory = 20
)

emapplot(
  ego_BP_down_69_SOD1,
  showCategory = 20
)

emapplot(
  ego_BP_up_69_TDP43,
  showCategory = 20
)

emapplot(
  ego_BP_up_07_TDP43,
  showCategory = 20
)

emapplot(
  ego_BP_down_07_TDP43,
  showCategory = 20
)


#GO Hierarchy Plots

goplot(ego_BP_down_02_C9)

goplot(ego_BP_down_69_SOD1)

goplot(ego_BP_up_69_TDP43)

goplot(ego_BP_up_07_TDP43)

goplot(ego_BP_down_07_TDP43)


#Category–Gene Heatmaps

heatplot(
  ego_BP_down_02_C9,
  showCategory = 20
)

heatplot(
  ego_BP_down_69_SOD1,
  showCategory = 20
)

heatplot(
  ego_BP_up_69_TDP43,
  showCategory = 20
)

heatplot(
  ego_BP_up_07_TDP43,
  showCategory = 20
)

heatplot(
  ego_BP_down_07_TDP43,
  showCategory = 20
)

