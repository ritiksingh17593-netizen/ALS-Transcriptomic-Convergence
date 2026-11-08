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
gse52202 <- read.delim("C:/Users/LENOVO/Desktop/ALS/Workflow for internship/Group A (iPSC derived motor neuron)/GSE52202/Raw counts/GSE52202_raw_counts_GRCh38.p13_NCBI.tsv")


#Making a simple metadata file
gse52202_met <- data.frame(
  SampleID = c(
    "GSM1261105",
    "GSM1261106",
    "GSM1261107",
    "GSM1261108",
    "GSM1261109",
    "GSM1261110",
    "GSM1261111",
    "GSM1261112"
  ),
  Condition = c(
    "Control",
    "Control",
    "Control",
    "Control",
    "C9ALS",
    "C9ALS",
    "C9ALS",
    "C9ALS"
  )
)



rownames(gse52202_met) <- gse52202_met$SampleID
rownames(gse52202_met)

rownames(gse52202) <- gse52202$GeneID
gse52202 <- gse52202[ ,-1]


#Does every sample name in the count matrix exist in the metadata?
all(colnames(gse52202)[-1] %in% rownames(gse52202_met))

#Are they not only present, but in exactly the same order?
all(colnames(gse52202) == rownames(gse52202_met))

#Designing the DESeq2 object
dds <- DESeqDataSetFromMatrix(
  countData = gse52202,
  colData = gse52202_met,
  design = ~ Condition
)
dds

#Removing low counts
keep <- rowSums(counts(dds)) >=10
dds <- dds[keep, ]

#Setting up reference
dds$Condition <- relevel(dds$Condition , ref = "Control")
resultsNames(deg)

#PCA
vsd <- vst(dds, blind = TRUE)

pcaData <- plotPCA(
  vsd,
  intgroup = "Condition",
  returnData = TRUE
)

ggplot(
  pcaData,
  aes(
    PC1,
    PC2,
    color = Condition,
    label = name
  )
) +
  geom_point(size = 4) +
  geom_text(vjust = -0.5)


#Running the analysis
register(SnowParam(8))
deg <- DESeq(dds, parallel = TRUE)

#normalized counts
sizeFactors(deg)
norm_counts <- counts(deg, normalized = TRUE)
write.csv(norm_counts, "Normalized conuts")


#Getting the results
res <- results(deg)
res
write.csv(res, "gse52202_result")

summary(res)

#Checks
sum(res$padj < 0.05, na.rm = TRUE)
sum(
  res$padj < 0.05 &
    abs(res$log2FoldChange) > 1,
  na.rm = TRUE
)

#Significant genes
sig_genes <- as.data.frame(res)

sig_genes <- sig_genes[
  !is.na(sig_genes$padj) &
    sig_genes$padj < 0.05,
]

sig_genes_df <- as.data.frame(sig_genes)
sig_genes_df
view(sig_genes_df)

write.csv(
  sig_genes_df,
  "GSE52202_Significant_Genes.csv"
)


up_genes <- sig_genes[
  sig_genes$log2FoldChange > 1,
]

down_genes <- sig_genes[
  sig_genes$log2FoldChange < -1,
]


#GeneID to Gene symbol
res_df <- as.data.frame(res)

res_df$symbol <- mapIds(
  org.Hs.eg.db,
  keys = rownames(res_df),
  column = "SYMBOL",
  keytype = "ENTREZID",
  multiVals = "first"
)


#Volcano plot
EnhancedVolcano(
  res_df,
  lab = res_df$symbol,
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 0.05,
  FCcutoff = 1,
  title = "GSE52202: C9ALS vs Control"
)


#Box plot
boxplot(
  log2(norm_counts + 1),
  las = 2
)


#MA plot
res_df$Significant <- ifelse(
  res_df$padj < 0.05,
  "Significant",
  "Not Significant"
)

ggplot(
  res_df,
  aes(
    x = log10(baseMean + 1),
    y = log2FoldChange,
    color = Significant
  )
) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_hline(yintercept = 0, linewidth = 0.5) +
  geom_hline(
    yintercept = c(-1, 1),
    linetype = "dashed"
  ) +
  scale_color_manual(
    values = c(
      "Not Significant" = "grey70",
      "Significant" = "red"
    )
  ) +
  labs(
    title = "MA Plot: GSE52202",
    x = "log10 Mean Normalized Expression",
    y = "Log2 Fold Change (C9ALS vs Control)"
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.title = element_blank(),
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    )
  )





