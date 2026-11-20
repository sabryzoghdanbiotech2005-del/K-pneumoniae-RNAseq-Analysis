# ==============================================================================
# 1. Setup & Libraries


# setwd("C:/Users/SIRA/OneDrive/Desktop/kpneumoniae-rnaseq-analysis")

library(DESeq2)
library(ggplot2)
library(EnhancedVolcano)
library(pheatmap)

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# 2. Data Import, Sample Matching & Reference Level
# ==============================================================================
counts_raw <- read.delim("data/GSE237490_raw_counts.txt", row.names = 1, check.names = FALSE)
counts <- counts_raw[, 4:ncol(counts_raw)]

coldata <- read.csv("data/metadata.csv", row.names = 1)

common_samples <- intersect(colnames(counts), rownames(coldata))
counts <- counts[, common_samples]
coldata <- coldata[common_samples, , drop = FALSE]

stopifnot(all(colnames(counts) == rownames(coldata)))

coldata$group <- relevel(factor(coldata$group), ref = "WT")
print(table(coldata$group))

# ==============================================================================
# 3. DESeq2 Filtering, Modeling & Contrasts (vs. WT)
# ==============================================================================
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = coldata,
  design = ~ group
)

keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]

cat("Samples:", ncol(dds), "\n")
cat("Genes remaining after filtering:", nrow(dds), "\n")

dds <- DESeq(dds)

res_CP17 <- results(dds, contrast = c("group", "CP17", "WT"), alpha = 0.05)
res_CP6  <- results(dds, contrast = c("group", "CP6",  "WT"), alpha = 0.05)
res_P13  <- results(dds, contrast = c("group", "P13",  "WT"), alpha = 0.05)
res_P2   <- results(dds, contrast = c("group", "P2",   "WT"), alpha = 0.05)

write.csv(as.data.frame(res_CP17), file = "results/tables/DEGs_CP17_vs_WT.csv")
write.csv(as.data.frame(res_CP6),  file = "results/tables/DEGs_CP6_vs_WT.csv")
write.csv(as.data.frame(res_P13),  file = "results/tables/DEGs_P13_vs_WT.csv")
write.csv(as.data.frame(res_P2),   file = "results/tables/DEGs_P2_vs_WT.csv")

count_DEGs <- function(res) {
  res_df <- as.data.frame(res)
  significant <- res_df[
    !is.na(res_df$padj) &
      res_df$padj < 0.05 &
      abs(res_df$log2FoldChange) > 1,
  ]
  c(
    Upregulated = sum(significant$log2FoldChange > 0),
    Downregulated = sum(significant$log2FoldChange < 0),
    Total = nrow(significant)
  )
}

results_table <- rbind(
  CP17 = count_DEGs(res_CP17),
  CP6  = count_DEGs(res_CP6),
  P13  = count_DEGs(res_P13),
  P2   = count_DEGs(res_P2)
)
print(results_table)
write.csv(results_table, file = "results/tables/DEG_counts_summary.csv")

key_genes <- c(
  "VK055_0479",
  "VK055_1763", "VK055_1764", "VK055_1765", "VK055_1766",
  "VK055_4636",
  "VK055_2558"
)

extract_key <- function(res, grp_name) {
  df <- as.data.frame(res)
  present <- intersect(key_genes, rownames(df))
  out <- df[present, c("log2FoldChange", "pvalue", "padj"), drop = FALSE]
  out$Comparison <- grp_name
  out$Gene <- rownames(out)
  out[, c("Comparison", "Gene", "log2FoldChange", "pvalue", "padj")]
}

key_summary <- rbind(
  extract_key(res_CP17, "CP17 vs WT"),
  extract_key(res_CP6,  "CP6 vs WT"),
  extract_key(res_P13,  "P13 vs WT"),
  extract_key(res_P2,   "P2 vs WT")
)
print(key_summary, row.names = FALSE)
write.csv(key_summary, file = "results/tables/key_genes_summary.csv", row.names = FALSE)

# ==============================================================================
# 4. Visualization: PCA, Volcano Plots & Top 50 Heatmap
# ==============================================================================
vsd <- vst(dds, blind = TRUE)

pcaData <- plotPCA(vsd, intgroup = "group", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

p_pca <- ggplot(pcaData, aes(PC1, PC2, color = group)) +
  geom_point(size = 4) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_bw() +
  ggtitle("PCA - K. pneumoniae RNA-seq")

ggsave("results/figures/01_PCA_QC.png", p_pca, width = 8, height = 6, dpi = 300)

res_list <- list(
  CP17 = res_CP17,
  CP6  = res_CP6,
  P13  = res_P13,
  P2   = res_P2
)

for (grp in names(res_list)) {
  res_obj <- as.data.frame(res_list[[grp]])
  p_vol <- EnhancedVolcano(
    res_obj,
    lab = rownames(res_obj),
    x = "log2FoldChange",
    y = "padj",
    title = paste(grp, "vs WT"),
    pCutoff = 0.05,
    FCcutoff = 1.0,
    pointSize = 3.0,
    labSize = 3.0,
    col = c("black", "black", "black", "red3"),
    colAlpha = 1
  )
  ggsave(
    paste0("results/figures/volcano_", grp, "_vs_WT.png"),
    p_vol, width = 10, height = 8, dpi = 300
  )
}

res_cp17_df <- as.data.frame(res_CP17)
res_valid <- res_cp17_df[
  !is.na(res_cp17_df$padj) &
    res_cp17_df$padj < 0.05 &
    abs(res_cp17_df$log2FoldChange) > 1,
]
res_sorted <- res_valid[order(res_valid$padj), ]
top50_genes <- head(rownames(res_sorted), 50)

mat <- assay(vsd)[top50_genes, , drop = FALSE]
mat_scaled <- t(scale(t(mat)))
mat_scaled <- mat_scaled[complete.cases(mat_scaled), , drop = FALSE]

pheatmap(
  mat_scaled,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  annotation_col = coldata["group", drop = FALSE],
  main = "Top 50 DEGs: CP17 vs WT",
  color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
  fontsize_row = 8,
  fontsize_col = 10,
  filename = "results/figures/heatmap_CP17_vs_WT.png",
  width = 8,
  height = 10
)

print("All analyses and plots completed successfully!")