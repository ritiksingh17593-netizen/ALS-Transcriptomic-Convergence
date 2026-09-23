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

#Loading up the count matrix
gse52202 <- read.delim("C:/Users/LENOVO/Desktop/ALS/Workflow for internship/Group A (iPSC derived motor neuron)/GSE52202/Raw counts/GSE52202_raw_counts_GRCh38.p13_NCBI.tsv")
gse210969 <- read.delim("C:/Users/LENOVO/Desktop/ALS/Workflow for internship/Group A (iPSC derived motor neuron)/GSE210969/Raw counts/GSE210969_raw_counts_GRCh38.p13_NCBI.tsv")
gse283507 <- read.csv("C:/Users/LENOVO/Desktop/ALS/Workflow for internship/Group A (iPSC derived motor neuron)/GSE283507/Raw counts/GSE283507_raw_Count_FPKM_TPM.csv")

#checking the dimentions of count matrix files
dim(gse52202)
dim(gse210969)
dim(gse283507)


#Some Important checks

#Finding how many genes are common in these 2 dataset
sum(gse210969$GeneID %in% gse283507$Gene_ID)   # length(gse969 ∩ gse507)

#Finding the intersection of genes in these 2 dataset 
length(intersect(
  gse210969$GeneID,
  gse52202$GeneID
))                                             # lenght(gse969 ∩ gse202)


#(A ∩ B) ∩ C, Genes common in all 3 datasets
common_genes <- intersect(
  intersect(gse210969$GeneID,
            gse52202$GeneID),
  gse283507$Gene_ID
)
length(common_genes)


#Selecting only relevant columns from gse283507
colnames(gse283507)
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
dim(gse283507_clean)


#Now creating the subset for all common genes
gse210969_sub <- gse210969[
  gse210969$GeneID %in% common_genes,
]

gse52202_sub <- gse52202[
  gse52202$GeneID %in% common_genes,
]

gse283507_sub <- gse283507_clean[
  gse283507_clean$Gene_ID %in% common_genes,
]

dim(gse210969_sub)
dim(gse52202_sub)
dim(gse283507_sub)

#Sorting the gene IDs
gse210969_sub <- gse210969_sub[order(gse210969_sub$GeneID), ]

gse283507_sub <- gse283507_sub[order(gse283507_sub$Gene_ID), ]

gse52202_sub <- gse52202_sub[order(gse52202_sub$GeneID), ]

#Checking the order of the gene IDs column
all(gse210969_sub$GeneID ==
      gse52202_sub$GeneID)

all(gse210969_sub$GeneID ==
      gse283507_sub$Gene_ID)


#extracting counts and mearging them
counts210969 <- gse210969_sub[, -1]
counts52202 <- gse52202_sub[, -1]
counts283507 <- gse283507_sub[, -1]


#Binding the columns together
merged_counts <- cbind(
  counts210969,
  counts52202,
  counts283507
)

#Assigning row names as gene IDs
rownames(merged_counts) <- gse210969_sub$GeneID
dim(merged_counts)
write.csv(merged_counts, "merged counts")

#Metadata loading
metadata <- read_excel("C:/Users/LENOVO/Desktop/ALS/Workflow for internship/Group A (iPSC derived motor neuron)/Master metadata/Master_metadata.xlsx")
dim(metadata)

head(colnames(merged_counts))
head(metadata$SampleID)

tail(colnames(merged_counts))
tail(metadata$SampleID)

all(metadata$SampleID %in% colnames(merged_counts))

metadata <- metadata[
  match(colnames(merged_counts),
        metadata$SampleID),
]

all(metadata$SampleID ==
      colnames(merged_counts))


rownames(metadata) <- metadata$SampleID


#Creating DESeq2 object
dds <- DESeqDataSetFromMatrix(
  countData = merged_counts,
  colData = metadata,
  design = ~ Batch + Condition
)
dds

#Design
table(metadata$Batch)
table(metadata$Condition)
table(metadata$Batch,
      metadata$Condition)


#PCA analysis

vsd <- vst(dds, blind = TRUE)

plotPCA(vsd, intgroup = "Batch")
plotPCA(vsd, intgroup = "Condition")

#Limma
vsd_corrected <- vst(dds, blind = FALSE)
assay(vsd_corrected) <- removeBatchEffect(
  assay(vsd_corrected),
  batch = vsd_corrected$Batch,
  design = model.matrix(~ Condition,
                        colData(vsd_corrected))
)
plotPCA(vsd_corrected,
        intgroup = "Condition")

plotPCA(vsd_corrected,
        intgroup = "Batch")



#Analysis
parallel::detectCores()

dds <- DESeq(
  dds,
  parallel = TRUE,
  BPPARAM = SnowParam(workers = 8)
)

#Result of the analysis
resultsNames(dds)
res <- results(
  dds,
  contrast = c(
    "Condition",
    "ALS",
    "Control"
  )
)
res
resultsNames(dds)
write.csv(res, "DESeq2_result")

#Checking levels
levels(dds$Condition)


#Turning ENTREZID to Gene Symbols
res_df <- as.data.frame(res)

res_df$ENTREZID <- rownames(res_df)

res_df$SYMBOL <- mapIds(
  org.Hs.eg.db,
  keys = res_df$ENTREZID,
  keytype = "ENTREZID",
  column = "SYMBOL",
  multiVals = "first"
)

head(res_df[, c(
  "ENTREZID",
  "SYMBOL",
  "log2FoldChange",
  "padj"
)])
head(res_df)

write_excel_csv(res_df, "result with gene symbols")


#Significant genes
sig_genes <- res_df %>%
  filter(
    padj < 0.05,
    abs(log2FoldChange) > 1
  ) %>%
  arrange(padj)
write.csv(sig_genes, "Significant genes")
nrow(sig_genes)
head(sig_genes)


als_genes <- c(
  "SOD1",
  "TARDBP",
  "FUS",
  "C9orf72",
  "TBK1",
  "NEK1",
  "OPTN",
  "SQSTM1",
  "VCP",
  "UBQLN2"
)

res_df[res_df$SYMBOL %in% als_genes, ]



#volcano plot
res_df$Significant <- ifelse(
  res_df$padj < 0.05 &
    abs(res_df$log2FoldChange) > 1,
  "Significant",
  "Not Significant"
)
table(res_df$Significant)

ggplot(
  res_df,
  aes(
    x = log2FoldChange,
    y = -log10(padj),
    color = Significant
  )
) +
  geom_point(alpha = 0.6)




#Heatmaps

sig_genes <- res_df %>%
  filter(
    padj < 0.05,
    abs(log2FoldChange) > 1
  ) %>%
  arrange(padj)
nrow(sig_genes)

mat <- assay(vsd)[top_genes, ]

gene_labels <- sig_genes$SYMBOL[1:30]

rownames(mat) <- gene_labels

mat.z <- t(scale(t(mat)))
annotation_col <- metadata[, c("Condition", "Batch")]
all(rownames(annotation_col) == colnames(mat.z))


pheatmap(
  mat.z,
  annotation_col = annotation_col,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  fontsize_row = 7,
  border_color = NA,
  main = "Top 30 DEGs",
  color = colorRampPalette(
    c("navy","white","firebrick3")
  )(100)
)


#MA plot
res_df <- res_df[!is.na(res_df$padj), ]

res_df$Significant <- res_df$padj < 0.05

ggplot(
  res_df,
  aes(
    x = baseMean,
    y = log2FoldChange,
    color = Significant
  )
) +
  geom_point(alpha = 0.5, size = 1) +
  scale_x_log10() +
  geom_hline(yintercept = 0,
             linetype = "dashed") +
  theme_minimal()


sum(res_df$padj < 0.05, na.rm = TRUE)






