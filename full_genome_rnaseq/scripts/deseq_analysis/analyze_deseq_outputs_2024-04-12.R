# analyze_deseq_outputs_2024-04-12.R

# install.packages('pheatmap')
# install DESeq if necessary
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# install and load package
# BiocManager::install("DESeq")


# Packages
library(pheatmap)
library(DESeq2)
library(stringr)
library(tidyverse)
library(RColorBrewer)
library(AnnotationDbi)
library(GO.db)

wd <- '~/Documents/geisinger_lab_research/bioinformatics_in_acinetobacter/rnaSeq/2024-04-03_pbpGlpsB_clean-redo'
setwd(wd)

pbpG_des_output <- read.table("./DES_pbpG_2024-04-04.csv", header=TRUE, sep = ',')

## Read in raw data and prepare feature count table
feature_count <- read.table("./counts.txt", header=TRUE, row.names = 1)
# Grab the count data
data <- feature_count[,6:14]
# view(data)

## Set up column names in feature count data AND prep metadata table with strains and conditions
# Trim column names down to just the sample IDs
column_names <- colnames(data)
column_names <- sub("X.work.geisingerlab.Mark.rnaSeq.2024.04.03_pbpGlpsB_clean.redo.data.mapped.", "", column_names)
column_names <- sub("Aligned.sortedByCoord.out.bam", "", column_names)
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

## Load data, pre-filter low-count genes, and relevel to set WT as the reference
des_data <- DESeqDataSetFromMatrix(countData = data, colData = meta, design = ~ condition)
smallestGroupSize <- 3  # should be size of smallest group; I did 3 replicates
keep <- rowSums(counts(des_data) >= 10) >= smallestGroupSize  # keep data where count is >10 in all 3 samples
des_data <- des_data[keep,]
# relevel dds$condition to set WT as the reference
des_data$condition <- relevel(des_data$condition, ref = "WT")

dds <- DESeq(des_data)
normalized_counts <- as.data.frame(counts(dds,normalized =TRUE))

#####
# Working on heatmap
# See this site: https://biostatsquid.com/step-by-step-heatmap-tutorial-with-pheatmap/
#####

# Prepare counts/DESeq data for heatmap
vsd <- vst(dds, blind=FALSE)

## Make a dataframe for significantly diff expressed genes with index = locus tag and a column "gene functions" containing GO Biological Process (BP) ontology term

# Start with pbpG. Get sig diff expressed UPREGULATED genes
pbpG_des_output <- pbpG_des_output[order(pbpG_des_output$padj),]
pbpG_sig_genes <- pbpG_des_output$X[pbpG_des_output$padj < 0.1 & pbpG_des_output$log2FoldChange > 0.5]

# Try KEGG
fxn_annos <- read.csv("/Users/mws/Documents/geisinger_lab_research/Abaum 17978-mff RefSeq Functional Annotation.csv")

# Construct dataframe
pbpG_up_fxnannos <- as.data.frame(pbpG_sig_genes)
colnames(pbpG_up_fxnannos) <- c("Gene")
pbpG_up_fxnannos <- left_join(pbpG_up_fxnannos, fxn_annos, by="Gene")
# # Count KEGG terms...
# pbpG_up_fxnannos %>% count(Category1)

# Make a dataframe that links gene with one KEGG annotation (from Category1 column of functional annotation spreadsheet)
gene_anno_df <- pbpG_up_fxnannos[, c('Gene', 'Category1')]
gene_anno_df$Gene <- paste("gene-", gene_anno_df$Gene, sep="")
gene_anno_df <- column_to_rownames(gene_anno_df, var = "Gene")
gene_anno_df$Category1 <- sub("^$", NA, gene_anno_df$Category1)
gene_anno_df[is.na(gene_anno_df)] <- "None"
colnames(gene_anno_df) <- c("KEGG term")

# select <- order(rowMeans(counts(dds,normalized=TRUE)),
#                 decreasing=TRUE)[1:40]
# pheatmap(assay(vsd)[select,], show_rownames=FALSE,)

## Heatmap; used VSD-transformed count data
# Set up colors
c25 <- c(
  "dodgerblue2", "#E31A1C", 
               "green4", 
               "#6A3D9A", 
               "#FF7F00", 
               "black", "gold1", 
               "skyblue2", "#FB9A99", 
               "palegreen2", 
               "#CAB2D6", 
               "#FDBF6F", 
               "gray70", "khaki2",
               "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
               "darkturquoise", "green1", "yellow4", "yellow3",
               "darkorange4", "brown"
               )

paste0(unique(gene_anno_df$`KEGG term`), '"', c25) # This needs tweaking

color_list <- list(
  `KEGG term` = c(
    None = "white", 
    `Cellular processes and signaling` = "#E31A1C",    
    `Membrane transport` = "green4",        
    `Various` = "#6A3D9A", 
    `Lipid metabolism` = "#FF7F00", 
    `Carbohydrate metabolism` = "black", 
    `Glycan biosynthesis and metabolism` = "gold1", 
    `Amino acid metabolism` = "skyblue2", 
    `Replication and repair`="#FB9A99", 
    `Transcription`="palegreen2", 
    `Metabolism`="#CAB2D6", 
    `Cell motility`="#FDBF6F", 
    `Poorly characterized`="gray70", 
    `Folding sorting and degradation`="khaki2", 
    `Cellular community - prokaryotes`="maroon", 
    `Xenobiotics biodegradation and metabolism`="orchid1", 
    `Enzyme families`="deeppink1", 
    `Translation`="blue1", 
    `Signal transduction`="steelblue4", 
    `Energy metabolism`="darkturquoise", 
    `Genetic information processing`="green1", 
    `Metabolism of cofactors and vitamins`="yellow4", 
    `Metabolism of other amino acids`="yellow3"
  ))

pbpG_sig_up_heatmap <- pheatmap(assay(vsd[rownames(gene_anno_df), ]), 
         annotation_row = gene_anno_df, 
         annotation_colors = color_list,
         show_rownames = F)
pbpG_sig_up_heatmap

tiff("./data/plot_heatmap_pbpG_sig_upregulated.tiff", width=1100, height=900, res = 150)
pbpG_sig_up_heatmap
dev.off()


# Function to save
save_pheatmap_tiff <- function(x, filename, width=1100, height=900, res = 300) {
  tiff(filename, width = width, height = height, res = res)
  dev.off()
}


save_pheatmap_tiff(pbpG_sig_up_heatmap, "./data/plot_heatmap_pbpG_sig_upregulated.tiff")

# Now get all diff expressed genes

# vsd <- vst(dds, blind=FALSE)
# fxn_annos <- read.csv("/Users/mws/Documents/geisinger_lab_research/Abaum 17978-mff RefSeq Functional Annotation.csv")

# Get sig diff expressed genes
pbpG_des_output <- pbpG_des_output[order(pbpG_des_output$padj),]
pbpG_sig_genes <- pbpG_des_output$X[pbpG_des_output$padj < 0.05 & abs(pbpG_des_output$log2FoldChange) > 1]
length(pbpG_sig_genes)  # 65

# Construct dataframe
pbpG_sig_genes <- pbpG_sig_genes[!is.na(pbpG_sig_genes)]  # drop rows with NA as locus tag
pbpG_sig_fxnannos <- as.data.frame(pbpG_sig_genes)
colnames(pbpG_sig_fxnannos) <- c("Gene")
pbpG_sig_fxnannos <- left_join(pbpG_sig_fxnannos, fxn_annos, by="Gene")
# # Count KEGG terms...
# pbpG_up_fxnannos %>% count(Category1)

# Make a dataframe that links gene with one KEGG annotation (from Category1 column of functional annotation spreadsheet)
gene_anno_df <- pbpG_sig_fxnannos[, c('Gene', 'Category1')]
gene_anno_df$Gene <- paste("gene-", gene_anno_df$Gene, sep="")
gene_anno_df <- column_to_rownames(gene_anno_df, var = "Gene")
gene_anno_df$Category1 <- sub("^$", NA, gene_anno_df$Category1)
gene_anno_df[is.na(gene_anno_df)] <- "None"
colnames(gene_anno_df) <- c("KEGG term")

# Set up colors
unique(gene_anno_df$`KEGG term`)  # 10 terms
# [1] "None"                                     
# [2] "Metabolism of cofactors and vitamins"     
# [3] "Amino acid metabolism"                    
# [4] "Glycan biosynthesis and metabolism"       
# [5] "Cellular processes and signaling"         
# [6] "Membrane transport"                       
# [7] "Various"                                  
# [8] "Xenobiotics biodegradation and metabolism"
# [9] "Cell motility"                            
# [10] "Carbohydrate metabolism"
color_list <- list(
  `KEGG term` = c(
    None = "white", 
    `Cellular processes and signaling` = "#E31A1C",    
    `Membrane transport` = "green4",        
    `Various` = "#6A3D9A", 
    `Carbohydrate metabolism` = "black", 
    `Glycan biosynthesis and metabolism` = "gold1", 
    `Amino acid metabolism` = "skyblue2", 
    `Cell motility`="#FDBF6F", 
    `Xenobiotics biodegradation and metabolism`="orchid1", 
    `Metabolism of cofactors and vitamins`="yellow4"
  ))

## Heatmap; used VSD-transformed count data
pbpG_sig_heatmap <- pheatmap(assay(vsd[rownames(gene_anno_df), ]), 
                             annotation_row = gene_anno_df,
                             annotation_colors = color_list, 
                             show_rownames = F,
                             main = "∆pbpG DEGs at padj < 0.05 and abs(l2FC) > 1")
pbpG_sig_heatmap

tiff("./data/plot_heatmap_pbpG_all_sig.tiff", width=1100, height=900, res = 150)
pbpG_sig_heatmap
dev.off()

# TODO: EG wants to see Z-score scaling.  Note that this may require pre-clustering: https://adairama.wordpress.com/2022/09/10/beware-of-using-the-scale-parameter-in-pheatmap-or-scaled-data-in-complexheatmap/ 


# #### GO annotation...
# 
# # Load term2gene with GO IDs and acx locus tags.  This was created from EG functional annotations spreadsheet.
# term2gene <- read.csv("/Users/mws/Documents/geisinger_lab_research/bioinformatics_in_acinetobacter/chip-seq_for_nicole/2024-02_nicole_chip_gsea_newtargetslist/term2files/term2gene_ab17978.csv")
# term2name <- read.csv("/Users/mws/Documents/geisinger_lab_research/bioinformatics_in_acinetobacter/chip-seq_for_nicole/2024-02_nicole_chip_gsea_newtargetslist/term2files/term2name.csv")
# 
# 
# # use gene_list to subset term2gene, grabbing all GO IDs that correspond to genes in gene_list.
# # Then, add go ontologies
# go_terms_in_genelist <- term2gene[term2gene$Gene %in% pbpG_sig_genes, ]
# go_terms_in_genelist <- go_terms_in_genelist[go_terms_in_genelist$GOID != "", ]
# # Get all GOID:ontology pairs (CC, BP, MF) from the gene list
# goid_to_ontology <- AnnotationDbi::select(GO.db, keys = go_terms_in_genelist$GOID, keytype="GOID", columns=c("ONTOLOGY"))
# # Apparently that contains a lot of duplicates... drop them
# goid_ontology_uniques <- goid_to_ontology %>% distinct(GOID, ONTOLOGY, .keep_all = TRUE)
# # Add ontologies to gene list with go terms.  Note that there are some "NA" ontologies somehow...
# go_table_ont <- left_join(go_terms_in_genelist, goid_ontology_uniques, by="GOID") 
# 
# # subset gene lists by filtering ontologies
# gene_table_bp <- go_table_ont[go_table_ont$ONTOLOGY == "BP",] %>% na.omit()
# gene_list_bp <- gene_table_bp$Gene
# gene_table_mf <- go_table_ont[go_table_ont$ONTOLOGY == "MF",] %>% na.omit()
# gene_list_mf <- gene_table_mf$Gene
# 
# # Add GO terms (e.g. "peptidoglycan turnover")
# gene_table_bp <- left_join(gene_table_bp, term2name, by="GOID")
# gene_table_mf <- left_join(gene_table_mf, term2name, by="GOID")
# 
# ### NOTE: I'm bailing here due to the large number of terms (65).  Would be ideal to narrow somehow, not sure.
# 
