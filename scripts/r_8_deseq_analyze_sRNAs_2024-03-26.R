# r_9_deseq_analyze_sRNAs_2024-03-26.R
# DESeq2

# Installed these packages using Packrat:

# install.packages("tidyverse")
# install.packages("vsn") ## not for this R version
# install.packages("ggrepel")
# install.packages("DEGreport") ## not for this R version
# install.packages("pheatmap")
# BiocManager::install(pattern="apeglm")
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("DESeq2")
#.libpaths()
#BiocManager::install("apeglm")

# Load packages via Packrat package manager
library(packrat)
pr_dir <- '/work/geisingerlab/Mark/rnaSeq/2024-01_rnaseq_pbpGlpsB/software'
setwd(pr_dir)
packrat::init()

# Load libraries
library(DESeq2)
library(ggplot2)
#library(vsn)
library(dplyr)
library(tidyverse)
library(ggrepel)
# library(DEGreport)
library(pheatmap)
library(apeglm)

wd <- "/work/geisingerlab/Mark/rnaSeq/stationary_phase_palethorpe_forNER_2024-03-04"
setwd(wd)

## Read in raw data and prepare feature count table
feature_count <- read.table("./data/featurecounts/srna_counts.txt", header=TRUE, row.names = 1)
# Grab the count data
data <- feature_count[,6:11]
# view(data)

## Set up column names in feature count data AND prep metadata table with strains and conditions
# Trim column names down to just the sample IDs
column_names <- colnames(data)
column_names <- sub("X.work.geisingerlab.Mark.rnaSeq.stationary_phase_palethorpe_forNER_2024.03.04.data.mapped..", "", column_names)
column_names <- sub("Aligned.sortedByCoord.out.bam", "", column_names)
column_names <- sub("SRR16949318", "17961-dbfmR_1", column_names)
column_names <- sub("SRR16949319", "17961-dbfmR_2", column_names)
column_names <- sub("SRR16949320", "17961-dbfmR_3", column_names)
column_names <- sub("SRR16949321", "17961-WT-1", column_names)
column_names <- sub("SRR16949322", "17961-WT-2", column_names)
column_names <- sub("SRR16949323", "17961-WT-3", column_names)
column_names
colnames(data) <- column_names
# Use regex to get condition (mutation) from strain IDs
# ".+?(?=_)":   .+?  means "match any character (.) any number of times (+?)"
# (?=_): a positive lookahead: find where the character - is matched
conditions <- str_extract(column_names, ".+?(?=_)")
# Use column_names and conditions to make metadata table
meta <- data.frame(column_names, conditions)
colnames(meta) <- c('id', 'condition')
meta  # Verify that IDs and conditions are as expected
# Write out metadata table (save!!)
#meta_filepath <- '/work/geisingerlab/Mark/rnaSeq/2024-01_rnaseq_pbpGlpsB/data'
#meta_file <- file.path(meta_filepath, 'metadata.txt')
#write.table(meta, meta_file, row.names=FALSE)  # TODO: Come back to this. Ok to have "id" label, or needs to be blank as in examples?  May want to change colnames before this.

## Load data, pre-filter low-count genes, and relevel to set WT as the reference
des_data <- DESeqDataSetFromMatrix(countData = data, colData = meta, design = ~ condition)
smallestGroupSize <- 3  # should be size of smallest group; I did 3 replicates
keep <- rowSums(counts(des_data) >= 10) >= smallestGroupSize  # keep data where count is >10 in all 3 samples
# relevel dds$condition to set WT as the reference
des_data$condition <- relevel(des_data$condition, ref = "WT")

dds <- DESeq(des_data)

# LFC shrinkage helps visualize and rank genes...
resultsNames(dds)  # get names of coefficients to shrink
## "Intercept" "condition_DlpsB_vs_WT" "condition_DpbpG_vs_WT"
# Use default shrinkage, which is apeglm method (Zhu, Ibrahim, and Love 2018)
resLFC_R <- lfcShrink(dds, coef="condition_dbfmR_vs_WT", type="apeglm")
resLFC_S <- lfcShrink(dds, coef="condition_dbfmS_vs_WT", type="apeglm")
resLFC_RS <- lfcShrink(dds, coef="condition_dbfmRS_vs_WT", type="apeglm")

# Clean up row names
#pbpG_rownames <- sub('gene-', '', row.names(resLFC_pbpG))
#lpsB_rownames <- sub('gene-', '', row.names(resLFC_lpsB))
#row.names(resLFC_pbpG) <- pbpG_rownames
#row.names(resLFC_lpsB) <- lpsB_rownames

# Count significant hits
sum(resLFC_R$padj < 0.1, na.rm=TRUE)
sum(resLFC_R$padj < 0.1 & abs(resLFC_R$log2FoldChange) > 1, na.rm=TRUE)

sum(resLFC_S$padj < 0.1, na.rm=TRUE)
sum(resLFC_S$padj < 0.1 & abs(resLFC_S$log2FoldChange) > 1, na.rm=TRUE)

sum(resLFC_RS$padj < 0.1, na.rm=TRUE)
sum(resLFC_RS$padj < 0.1 & abs(resLFC_RS$log2FoldChange) > 1, na.rm=TRUE)

# Check largest foldchanges from significant hits
resSig_R <- resLFC_R[which(resLFC_R$padj < 0.1), ]
resSig_S <- resLFC_S[which(resLFC_S$padj < 0.1), ]
resSig_RS <- resLFC_RS[which(resLFC_RS$padj < 0.1), ]

#downregulated
head(resSig_R[order(resSig_R$log2FoldChange),])
head(resSig_S[order(resSig_S$log2FoldChange),])
head(resSig_RS[order(resSig_RS$log2FoldChange),])

#upregulated
tail(resSig_R[order(resSig_R$log2FoldChange),])
tail(resSig_S[order(resSig_S$log2FoldChange),])
tail(resSig_RS[order(resSig_RS$log2FoldChange),])

# MA plot
#plotMA(resSig, ylim=c(-1, 1))

# Check PCA plot
#rld <- rlog(des, blind=TRUE)
#plotPCA(rld, intgroup="condition") + geom_text(aes(label=name))

# Plot the gene with the smallest value 
#plotCounts(dds, gene=which.min(res$padj), intgroup="condition")
# plotCounts(dds, gene="gene-ACX60_RS15100", intgroup="condition")

# Save data
#file_name = 'NAME/OF/FILE'
#write.table(<results/go/here>, file_name, quote=FALSE)
# Can also pre-filter for significant hits
#padj.cutoff <- 0.05 # False Discovery Rate cutoff
#significant_results <- res[which(res$padj < padj.cutoff),]
#write.table(significant_results, file_name, quote=FALSE)

#get session information, including package versions... can save to a variable, then write out.
#print(sessionInfo())


# Save current analyses
out_dir <- "/work/geisingerlab/Mark/rnaSeq/2024-03_eg_2020_rnaseq/data"
setwd(out_dir)

resOrdered_R <- resLFC_R[order(resLFC_R$pvalue),]
write.csv(as.data.frame(resOrdered_R), file="DES_d-bfmR-sRNAs_2024-03-07.csv")

resOrdered_S <- resLFC_S[order(resLFC_S$pvalue),]
write.csv(as.data.frame(resOrdered_S), file="DES_d-bfmS-sRNAs_2024-03-07.csv")

resOrdered_RS <- resLFC_RS[order(resLFC_RS$pvalue),]
write.csv(as.data.frame(resOrdered_RS), file="DES_d-bfmRS-sRNAs_2024-03-07.csv")
