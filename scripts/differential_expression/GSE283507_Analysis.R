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

#Loading files
gse283507 <- read.csv("C:/Users/LENOVO/Desktop/ALS/GSE283507/Raw counts/GSE283507_raw_Count_FPKM_TPM.csv")
gse283507_clean <- gse283507 %>%
  dplyr::select(
    Gene_ID,
    RC802DMSO6H.1_Read_Count,
    RC802DMSO6H.2_Read_Count,
    RC802DMSO6H.3_Read_Count,
    TDP43DMSO6H.1_Read_Count,
    TDP43DMSO6H.2_Read_Count,
    TDP43DMSO6H.3_Read_Count
  )
colnames(gse283507_clean)
colnames(gse283507_clean) <- c(
  "Gene_ID",
  "GSM8664404",
  "GSM8664405",
  "GSM8664406",
  "GSM8664407",
  "GSM8664408",
  "GSM8664409"
)


#Metadata
gse283507_met <- data.frame(
  SampleID = c(
    "GSM8664404","GSM8664405","GSM8664406",
    "GSM8664407","GSM8664408","GSM8664409"
  ),
  Condition = c(
    "Control","Control","Control",
    "TDP43","TDP43","TDP43"
  )
)
rownames(gse283507_met) <- gse283507_met$SampleID


#Does every sample name in the count matrix exist in the metadata?
all(colnames(gse283507_clean)[-1] %in% rownames(gse283507_met))

#Are they not only present, but in exactly the same order?
all(colnames(gse283507_clean)[-1] == rownames(gse283507_met))


# Move Gene_ID into rownames
rownames(gse283507_clean) <- gse283507_clean$Gene_ID
gse283507_clean <- gse283507_clean[, -1]

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(
  countData = gse283507_clean,
  colData = gse283507_met,
  design = ~ Condition
)

# Remove low-count genes
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

# Set reference level
dds$Condition <- relevel(dds$Condition, ref = "Control")

#PCA
vsd <- vst(dds, blind = TRUE)
plotPCA(vsd, intgroup = "Condition")

# Run DESeq2
dds <- DESeq(dds)


#Extracting the results
resultsNames(dds)
res <- results(dds)

res_TDP43 <- results(
  dds,
  contrast = c("Condition", "TDP43", "Control")
)

res_TDP43_df <- as.data.frame(res_TDP43)
write.csv(res_TDP43_df, "Result with gene symbols")

#normalized counts
sizeFactors(dds)
norm_counts <- counts(dds, normalized = TRUE)
write.csv(norm_counts, "Normalized conuts")

#Checks
sum(res$padj < 0.05, na.rm = TRUE)
sum(
  res$padj < 0.05 &
    abs(res$log2FoldChange) > 1,
  na.rm = TRUE
)

#Significant genes
sig_genes <- as.data.frame(res_TDP43_df)

sig_genes <- sig_genes[
  !is.na(sig_genes$padj) &
    sig_genes$padj < 0.05,
]

write.csv(
  sig_genes,
  "GSE283507_Significant_Genes.csv"
)


up_genes <- sig_genes[
  sig_genes$log2FoldChange > 1,
]

down_genes <- sig_genes[
  sig_genes$log2FoldChange < -1,
]


#GeneID to Gene symbol
res_TDP43_df$symbol <- mapIds(
  org.Hs.eg.db,
  keys = rownames(res_TDP43_df),
  column = "SYMBOL",
  keytype = "ENTREZID",
  multiVals = "first"
)

#Volcano plot
EnhancedVolcano(
  res_TDP43_df,
  lab = res_TDP43_df$symbol,
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 0.05,
  FCcutoff = 1,
  title = "GSE283507: TDP43 vs Control"
)


#MA plot
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


## Top 50 most significant genes
top_genes <- head(
  rownames(res_TDP43[order(res_TDP43$padj), ]),
  50
)

# Expression matrix
mat <- assay(vsd)[top_genes, ]

# Heatmap
pheatmap(
  mat,
  scale = "row",
  annotation_col = gse283507_met
)


