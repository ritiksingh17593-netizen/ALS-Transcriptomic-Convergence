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


#Loading up the count matrix
gse210969 <- read.delim("C:/Users/LENOVO/Desktop/Workflow for internship/Group A (iPSC derived motor neuron)/GSE210969/Raw counts/GSE210969_raw_counts_GRCh38.p13_NCBI.tsv")


#Metadata
gse210969_met <- data.frame(
  SampleID = c(
    "GSM6443466",
    "GSM6443467",
    "GSM6443468",
    "GSM6443469",
    "GSM6443470",
    "GSM6443471",
    "GSM6443472",
    "GSM6443473"
  ),
  
  Condition = c(
    "Control",
    "Control",
    "Control",
    "SOD1",
    "SOD1",
    "TDP43",
    "TDP43",
    "TDP43"
  )
)
rownames(gse210969_met) <- gse210969_met$SampleID

rownames(gse210969) <- gse210969$GeneID
gse210969 <- gse210969[ ,-1]

#Does every sample name in the count matrix exist in the metadata?
all(colnames(gse210969)[-1] %in% rownames(gse210969_met))

#Are they not only present, but in exactly the same order?
all(colnames(gse210969) == rownames(gse210969_met))

#Designing the DESeq2 object
dds <- DESeqDataSetFromMatrix(
  countData = gse210969,
  colData = gse210969_met,
  design = ~ Condition
)
dds


#Setting up reference
dds$Condition <- relevel(dds$Condition , ref = "Control")
levels(dds$Condition)


#PCA
vsd <- vst(dds, blind = TRUE)
plotPCA(vsd, intgroup = "Condition")


#Sample to sample distance Heatmap
sampleDists <- dist(t(assay(vsd)))

sampleDistMatrix <- as.matrix(sampleDists)

pheatmap(
  sampleDistMatrix,
  clustering_distance_rows = sampleDists,
  clustering_distance_cols = sampleDists
)


#DESeq2 analysis
dds <- DESeq(dds)

resultsNames(dds)

#Extracting the results
res_SOD1 <- results(
  dds,
  contrast = c("Condition", "SOD1", "Control")
)
write.csv(res_SOD1, "Result for SOD1")
sum(res_SOD1$padj<0.05, na.rm = TRUE)

res_TDP43 <- results(
  dds,
  contrast = c("Condition", "TDP43", "Control")
)
write.csv(res_TDP43, "Result for TDP43")
sum(res_TDP43$padj<0.05, na.rm = TRUE)

#Significant genes
deg_SOD1 <- subset(
  as.data.frame(res_SOD1),
  padj < 0.05 & abs(log2FoldChange) > 1
)
nrow(deg_SOD1)
write.csv(deg_SOD1, "gse210969_sigGenes_SOD1")


deg_TDP43 <- subset(
  as.data.frame(res_TDP43),
  padj < 0.05 & abs(log2FoldChange) > 1
)
nrow(deg_TDP43)
write.csv(deg_TDP43, "gse210969_sigGenes_TDP43")


#upregulated gene
#SOD1 gene
upGenes_SOD1 <- as.data.frame(res_SOD1)
upGenes_SOD1 <- upGenes_SOD1[
  !is.na(upGenes_SOD1$padj) &
    upGenes_SOD1$padj < 0.05 &
    upGenes_SOD1$log2FoldChange > 1,
]
nrow(upGenes_SOD1)
write.csv(upGenes_SOD1, "upregulated gene SOD1")

#TDp43
upGene_TDP43 <- as.data.frame(res_TDP43)
upGene_TDP43 <- upGene_TDP43[
  !is.na(upGene_TDP43$padj) &
    upGene_TDP43$padj < 0.05 &
    upGene_TDP43$log2FoldChange > 1,
]
nrow(upGene_TDP43)
write.csv(upGene_TDP43, "upregulated gene TDP43")


#Down regulated genes
#SOD1 gene
downGenes_SOD1 <- as.data.frame(res_SOD1)
downGenes_SOD1 <- downGenes_SOD1[
  !is.na(downGenes_SOD1$padj) &
    downGenes_SOD1$padj < 0.05 &
    downGenes_SOD1$log2FoldChange < -1,
]
nrow(downGenes_SOD1)
write.csv(downGenes_SOD1, "downregulated genes SOD1")


#TDP43
downGene_TDP43 <- as.data.frame(res_TDP43)
downGene_TDP43 <- downGene_TDP43[
  !is.na(downGene_TDP43$padj) &
    downGene_TDP43$padj < 0.05 &
    downGene_TDP43$log2FoldChange < -1,
]
nrow(downGene_TDP43)
write.csv(downGene_TDP43, "downregulated genes TDP43")



#GeneID to symbol
res_SOD1_df <- as.data.frame(res_SOD1)

res_SOD1_df$Symbol <- mapIds(
  org.Hs.eg.db,
  keys = rownames(res_SOD1_df),
  column = "SYMBOL",
  keytype = "ENTREZID",
  multiVals = "first"
)

res_TDP43_df <- as.data.frame(res_TDP43)

res_TDP43_df$Symbol <- mapIds(
  org.Hs.eg.db,
  keys = rownames(res_TDP43_df),
  column = "SYMBOL",
  keytype = "ENTREZID",
  multiVals = "first"
)


#Volcano plot

#SOD1 vs Control
EnhancedVolcano(
  res_SOD1_df,
  lab = res_SOD1_df$Symbol,
  x = "log2FoldChange",
  y = "padj",
  title = "SOD1 vs Control",
  pCutoff = 0.05,
  FCcutoff = 1
)

#TDP43 vs Control
EnhancedVolcano(
  res_TDP43_df,
  lab = res_TDP43_df$Symbol,
  x = "log2FoldChange",
  y = "padj",
  title = "TDP43 vs Control",
  pCutoff = 0.05,
  FCcutoff = 1
)



#Top genes
top_genes <- head(
  rownames(res_SOD1[order(res_SOD1$padj), ]),
  50
)

#Heatmaps
mat <- assay(vsd)[top_genes, ]

pheatmap(
  mat,
  scale = "row",
  annotation_col = gse210969_met
)


#MA plot

#SOD1 vs control
ggplot(
  res_SOD1_df,
  aes(
    x = log10(baseMean + 1),
    y = log2FoldChange
  )
) +
  geom_point(
    aes(color = padj < 0.05),
    alpha = 0.5,
    size = 1
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  scale_color_manual(
    values = c("grey70", "red")
  ) +
  labs(
    title = "MA Plot: SOD1 vs Control",
    x = "log10 Mean Expression",
    y = "Log2 Fold Change"
  ) +
  theme_bw()


#TDP43 vs control
ggplot(
  res_TDP43_df,
  aes(
    x = log10(baseMean + 1),
    y = log2FoldChange
  )
) +
  geom_point(
    aes(color = padj < 0.05),
    alpha = 0.5,
    size = 1
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  scale_color_manual(
    values = c("grey70", "red")
  ) +
  labs(
    title = "MA Plot: TDP43 vs Control",
    x = "log10 Mean Expression",
    y = "Log2 Fold Change"
  ) +
  theme_bw()
