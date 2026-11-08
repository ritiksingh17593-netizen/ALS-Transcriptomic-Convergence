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

gse52202_genes <- read.csv("gse52202_SigGenes")
gse210969_genes <- read.csv("gse210969_sigGenes_TDP43")
gse283507_genes <- read.csv("gse283507_sigGenes")

dim(gse52202_genes)
dim(gse283507_genes)
dim(gse210969_genes)


colnames(gse52202_genes)[1] <- "EntrezID"
colnames(gse210969_genes)[1] <- "EntrezID"
colnames(gse283507_genes)[1] <- "EntrezID"

genes_52202  <- gse52202_genes$EntrezID
genes_210969 <- gse210969_genes$EntrezID
genes_283507 <- gse283507_genes$EntrezID

#Pairwise overlap
common_52202_210969 <- intersect(genes52202, genes_210969)
length(common_52202_210969)

common_52202_283507 <- intersect(genes_52202, genes_283507)
length(common_52202_283507)

common_210969_283507 <- intersect(genes_210969, genes_283507)
length(common_210969_283507)


#Result file
gse52202_res <- read.csv("gse52202_DESeq2_result")
gse210969_res <- read.csv("Result for TDP43")
gse283507_res <- read.csv("Result with gene symbols")

common_res_genes <- Reduce(intersect, list(gse52202_res$X,
                                           gse210969_res$X,
                                           gse283507_res$X))
common_res_genes
length(common_res_genes)

common_res_genes_sigGenes <- as.data.frame(common_res_genes)

common_res_genes_sigGenes <- common_res_genes_sigGenes[
  !is.na(common_res_genes_sigGenes$padj) &
    common_res_genes_sigGenes$padj < 0.05 & abs(common_res_genes_sigGenes$log2FoldChange) > 1,
]















