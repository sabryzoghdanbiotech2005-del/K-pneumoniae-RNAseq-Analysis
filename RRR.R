install.packages(c("ggplot2", "pheatmap"))
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("DESeq2", "EnhancedVolcano"))
BiocManager::install("DESeq2", force = TRUE)

# ============================================
# Script 01: Load and Explore Data
# ============================================

# 1. Set the workflow path (make sure this path is exactly your project path)
getwd() 
setwd("C:/Users/SIRA/OneDrive/Desktop/kpneumoniae-rnaseq-analysis")

# 2. Importing libraries
library(DESeq2)
library(ggplot2)


# 3. Reading the Counts file
# (We'll ignore the first 3 columns Chr and Length because they're not numeric expressions)
counts_raw <- read.delim("data/GSE237490_raw_counts.txt", row.names = 1, check.names = FALSE)
counts <- counts_raw[, 4:ncol(counts_raw)]

# 4. Reading the Metadata File
coldata <- read.csv("data/metadata.csv", row.names = 1)

# 5. Make sure the sample names match
print("Do sample names match?")
print(all(colnames(counts) == rownames(coldata)))

# 6. Arrange the columns so they match the metadata
counts <- counts[, rownames(coldata)]

# 7. Preparing a DESeq2 Object
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = coldata,
                              design = ~ group)

# 8. Filtering out genes with very low expression (less than 10 in at least 3 samples)
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]
print(paste("Genes remaining after filtering:", nrow(dds)))

# 9. Transforming data using VST (so we can do PCA)
vsd <- vst(dds, blind = TRUE)

# 10. Extracting PCA data and plotting it
pcaData <- plotPCA(vsd, intgroup = "group", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

p <- ggplot(pcaData, aes(PC1, PC2, color = group)) +
  geom_point(size = 4) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_bw() +
  ggtitle("PCA - K. pneumoniae RNA-seq")

# 11. Save the drawing
ggsave("results/figures/01_PCA_QC.png", p, width = 8, height = 6, dpi = 300)
print("PCA plot saved successfully in results/figures/!")

colnames(counts)
rownames(coldata)

# Let's look at the columns shared between the two files
common_samples <- intersect(colnames(counts), rownames(coldata))

# Let's filter the two files so they end up in exactly the same order
counts <- counts[, common_samples]
coldata <- coldata[common_samples, , drop = FALSE]

# Let's make sure
print(all(colnames(counts) == rownames(coldata)))


# ============================================
# Script 02: Differential Expression Analysis
# ============================================

# 1. Run the differential analysis (basic model)
dds <- DESeq(dds)

# 2. Extracting the results for each group compared to the WT
# First comparison: CP17 vs WT
res_CP17 <- results(dds, contrast = c("group", "CP17", "WT"))
summary(res_CP17)

# Second Comparison: CP6 vs WT
res_CP6 <- results(dds, contrast = c("group", "CP6", "WT"))
summary(res_CP6)

# Third Comparison: P2 vs WT
res_P2 <- results(dds, contrast = c("group", "P2", "WT"))
summary(res_P2)

# Comparison Four: P13 vs WT
res_P13 <- results(dds, contrast = c("group", "P13", "WT"))
summary(res_P13)

#3. Saving the results in files so we can use them for plotting
write.csv(as.data.frame(res_CP17), file = "results/tables/DEGs_CP17_vs_WT.csv")
write.csv(as.data.frame(res_CP6), file = "results/tables/DEGs_CP6_vs_WT.csv")
write.csv(as.data.frame(res_P2), file = "results/tables/DEGs_P2_vs_WT.csv")
write.csv(as.data.frame(res_P13), file = "results/tables/DEGs_P13_vs_WT.csv")

print("Differential expression analysis complete and saved!")


# ============================================
# Script 02: Volcano Plot for DEGs
# ============================================

# 1. Importing Libraries
library(ggplot2)
library(EnhancedVolcano)

# 2. Reading the results file to compare CP17 vs WT
# (Make sure you're in the right project path)
res_CP17 <- read.csv("results/tables/DEGs_CP17_vs_WT.csv", row.names = 1)

# 3. Filtering important genes (|log2FC| > 1 and padj < 0.05)
# (We'll deal with padj values that are NA)
res_CP17_filtered <- res_CP17[!is.na(res_CP17$padj) & res_CP17$padj < 0.05 & abs(res_CP17$log2FoldChange) > 1, ]

# 4. Creating the Volcano Plot using EnhancedVolcano
# (This makes a very professional-looking plot)
p <- EnhancedVolcano(res_CP17,
                     lab = rownames(res_CP17),
                     x = 'log2FoldChange',
                     y = 'padj',
                     title = 'CP17 vs WT',
                     pCutoff = 0.05,
                     FCcutoff = 1.0,
                     pointSize = 3.0,
                     labSize = 3.0,
                     col=c('black', 'black', 'black', 'red3'),
                     colAlpha = 1)

# 5. Save the drawing in high quality
ggsave("results/figures/volcano_CP17_vs_WT.png", p, width = 10, height = 8, dpi = 300)

# 6. Printing the number of important genes
cat("Number of significant DEGs (|log2FC| > 1 & padj < 0.05):", nrow(res_CP17_filtered), "\n")


# ============================================
# Script 03: Heatmap for Top 50 DEGs
# ============================================

# 1. Calling the libraries
library(pheatmap)

# 2. Preparing the VST data (the one we did before)
# (If it's not available in the Environment, we'll recalculate it quickly)
if (!exists("vsd")) {
  library(DESeq2)
  dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ group)
  dds <- dds[rowSums(counts(dds) >= 10) >= 3, ]
  vsd <- vst(dds, blind = TRUE)
}

# 3. Extract the top 50 genes from the CP17 results
# (we sort by p-value and take the top 50)
res_sorted <- res_CP17[order(res_CP17$padj), ]
top50_genes <- rownames(res_sorted)[1:50]

# 4. Extracting expression data for the top 50 genes
mat <- assay(vsd)[top50_genes, ]

# 5. Scaling data (Z-score) so the comparison is clearer in the Heatmap
mat_scaled <- t(scale(t(mat)))

# 6. Creating the Heatmap
pheatmap(mat_scaled,
         cluster_rows = TRUE,      # we group similar genes
         cluster_cols = TRUE,      # we group similar samples
         show_rownames = TRUE,     # we show gene names
         show_colnames = TRUE,     # we show sample names
annotation_col = coldata["group"], # we color samples by group
main = "Top 50 DEGs: CP17 vs WT"
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50), # professional color gradient
         fontsize_row = 8,
         fontsize_col = 10,
         filename = "results/figures/heatmap_CP17_vs_WT.png", # direct save
         width = 8, height = 10)
cat("Heatmap saved successfully in results/figures/heatmap_CP17_vs_WT.png\n")


# ============================================
# Script 04: Volcano Plots for other groups
# ============================================

# List of files we want to analyze
groups <- c("CP6", "P2", "P13")

for (grp in groups) {
# Reading the file
  res_file <- paste0("results/tables/DEGs_", grp, "_vs_WT.csv")
  
  if (file.exists(res_file)) {
    res_data <- read.csv(res_file, row.names = 1)
    
    # Making the drawing
    p_volcano <- EnhancedVolcano(res_data,
                                 lab = rownames(res_data),
                                 x = 'log2FoldChange',
                                 y = 'padj',
                                 title = paste(grp, "vs WT"),
                                 pCutoff = 0.05,
                                 FCcutoff = 1.0,
                                 pointSize = 3.0,
                                 labSize = 3.0,
                                 col=c('black', 'black', 'black', 'red3'),
                                 colAlpha = 1)
    
    # Save the drawing
    ggsave(paste0("results/figures/volcano_", grp, "_vs_WT.png"), p_volcano, width = 10, height = 8, dpi = 300)
    
    # Print the number of different genes
    sig_count <- sum(!is.na(res_data$padj) & res_data$padj < 0.05 & abs(res_data$log2FoldChange) > 1)
    cat("Done for", grp, "- Significant DEGs:", sig_count, "\n")
  }
}



# ============================================
# Script 05: Extract Top Genes per Group
# ============================================

# List of groups
groups <- c("CP17", "CP6", "P2", "P13")

for (grp in groups) {
  
# Reading the file
  res_file <- paste0("results/tables/DEGs_", grp, "_vs_WT.csv")
  
  if (file.exists(res_file)) {
    res_data <- read.csv(res_file, row.names = 1)
    
# Filtering important genes
    sig_genes <- res_data[!is.na(res_data$padj) & res_data$padj < 0.05 & abs(res_data$log2FoldChange) > 2, ]
    
# Sorting by log2FoldChange (to highlight the strongest changes)
    sig_genes <- sig_genes[order(sig_genes$log2FoldChange, decreasing = TRUE), ]
    
# Save the table
    output_file <- paste0("results/tables/Top_DEGs_", grp, ".csv")
    write.csv(sig_genes, output_file)
    
    cat("Done for", grp, "- Top DEGs saved to", output_file, "\n")
  }
}


#####################################################################
# Reading and Displaying Important Genes
top_cp17 <- read.csv("results/tables/Top_DEGs_CP17.csv")
top_cp6  <- read.csv("results/tables/Top_DEGs_CP6.csv")
top_p2   <- read.csv("results/tables/Top_DEGs_P2.csv")
top_p13  <- read.csv("results/tables/Top_DEGs_P13.csv")

cat("=== Top 10 Genes CP17 ===\n")
print(head(top_cp17[, c("X", "log2FoldChange", "padj")], 10))

cat("\n=== Top 10 Genes CP6 ===\n")
print(head(top_cp6[, c("X", "log2FoldChange", "padj")], 10))

cat("\n=== Top 5 Genes P2 ===\n")
print(head(top_p2[, c("X", "log2FoldChange", "padj")], 5))

cat("\n=== Top 5 Genes P13 ===\n")
print(head(top_p13[, c("X", "log2FoldChange", "padj")], 5))
