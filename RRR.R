install.packages(c("ggplot2", "pheatmap"))
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("DESeq2", "EnhancedVolcano"))
BiocManager::install("DESeq2", force = TRUE)

# ============================================
# Script 01: Load and Explore Data
# ============================================

# 1. ضبط مسار العمل (تأكد إن المسار ده هو مسار مشروعك بالظبط)
getwd() 
setwd("C:/Users/SIRA/OneDrive/Desktop/kpneumoniae-rnaseq-analysis")

# 2. استدعاء المكتبات
library(DESeq2)
library(ggplot2)


# 3. قراءة ملف الـ Counts
# (هنتجاهل أول 3 أعمدة Chr و Length عشان هما مش أرقام تعبير)
counts_raw <- read.delim("data/GSE237490_raw_counts.txt", row.names = 1, check.names = FALSE)
counts <- counts_raw[, 4:ncol(counts_raw)]

# 4. قراءة ملف الـ Metadata
coldata <- read.csv("data/metadata.csv", row.names = 1)

# 5. التأكد إن أسماء العينات متطابقة
print("Do sample names match?")
print(all(colnames(counts) == rownames(coldata)))

# 6. ترتيب الأعمدة عشان تتطابق مع الميتاداتا
counts <- counts[, rownames(coldata)]

# 7. تجهيز كائن DESeq2
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = coldata,
                              design = ~ group)

# 8. فلترة الجينات اللي تعبيرها قليل جداً (أقل من 10 في 3 عينات على الأقل)
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]
print(paste("Genes remaining after filtering:", nrow(dds)))

# 9. تحويل البيانات باستخدام VST (عشان نعمل PCA)
vsd <- vst(dds, blind = TRUE)

# 10. استخراج بيانات الـ PCA ورسمها
pcaData <- plotPCA(vsd, intgroup = "group", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

p <- ggplot(pcaData, aes(PC1, PC2, color = group)) +
  geom_point(size = 4) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_bw() +
  ggtitle("PCA - K. pneumoniae RNA-seq")

# 11. حفظ الرسمة
ggsave("results/figures/01_PCA_QC.png", p, width = 8, height = 6, dpi = 300)
print("PCA plot saved successfully in results/figures/!")

colnames(counts)
rownames(coldata)

# بنشوف الأعمدة المشتركة بين الملفين
common_samples <- intersect(colnames(counts), rownames(coldata))

# نفلتر الملفين عشان يبقوا في نفس الترتيب بالظبط
counts <- counts[, common_samples]
coldata <- coldata[common_samples, , drop = FALSE]

# نتأكد
print(all(colnames(counts) == rownames(coldata)))


# ============================================
# Script 02: Differential Expression Analysis
# ============================================

# 1. تشغيل التحليل التفاضلي (الموديل الأساسي)
dds <- DESeq(dds)

# 2. استخراج النتائج لكل مجموعة مقارنة بالـ WT
# المقارنة الأولى: CP17 vs WT
res_CP17 <- results(dds, contrast = c("group", "CP17", "WT"))
summary(res_CP17)

# المقارنة الثانية: CP6 vs WT
res_CP6 <- results(dds, contrast = c("group", "CP6", "WT"))
summary(res_CP6)

# المقارنة الثالثة: P2 vs WT
res_P2 <- results(dds, contrast = c("group", "P2", "WT"))
summary(res_P2)

# المقارنة الرابعة: P13 vs WT
res_P13 <- results(dds, contrast = c("group", "P13", "WT"))
summary(res_P13)

# 3. حفظ النتائج في ملفات عشان نستخدمها في الرسم
write.csv(as.data.frame(res_CP17), file = "results/tables/DEGs_CP17_vs_WT.csv")
write.csv(as.data.frame(res_CP6), file = "results/tables/DEGs_CP6_vs_WT.csv")
write.csv(as.data.frame(res_P2), file = "results/tables/DEGs_P2_vs_WT.csv")
write.csv(as.data.frame(res_P13), file = "results/tables/DEGs_P13_vs_WT.csv")

print("Differential expression analysis complete and saved!")


# ============================================
# Script 02: Volcano Plot for DEGs
# ============================================

# 1. استدعاء المكتبات
library(ggplot2)
library(EnhancedVolcano)

# 2. قراءة ملف النتائج لمقارنة CP17 vs WT
# (تأكد إنك في المسار الصح للمشروع)
res_CP17 <- read.csv("results/tables/DEGs_CP17_vs_WT.csv", row.names = 1)

# 3. فلترة الجينات المهمة (|log2FC| > 1 و padj < 0.05)
# (هنتعامل مع قيم padj اللي فيها NA)
res_CP17_filtered <- res_CP17[!is.na(res_CP17$padj) & res_CP17$padj < 0.05 & abs(res_CP17$log2FoldChange) > 1, ]

# 4. عمل الـ Volcano Plot باستخدام EnhancedVolcano
# (ده بيعمل رسمة احترافية جداً)
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

# 5. حفظ الرسمة بجودة عالية
ggsave("results/figures/volcano_CP17_vs_WT.png", p, width = 10, height = 8, dpi = 300)

# 6. طباعة عدد الجينات المهمة
cat("Number of significant DEGs (|log2FC| > 1 & padj < 0.05):", nrow(res_CP17_filtered), "\n")


# ============================================
# Script 03: Heatmap for Top 50 DEGs
# ============================================

# 1. استدعاء المكتبات
library(pheatmap)

# 2. تجهيز بيانات الـ VST (اللي عملناها قبل كده)
# (لو مش موجودة في الـ Environment، هنعيد حسابها بسرعة)
if (!exists("vsd")) {
  library(DESeq2)
  dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ group)
  dds <- dds[rowSums(counts(dds) >= 10) >= 3, ]
  vsd <- vst(dds, blind = TRUE)
}

# 3. استخراج أعلى 50 جين من نتائج CP17
# (بنرتب حسب الـ p-value وناخد أول 50)
res_sorted <- res_CP17[order(res_CP17$padj), ]
top50_genes <- rownames(res_sorted)[1:50]

# 4. استخراج بيانات التعبير لأفضل 50 جين
mat <- assay(vsd)[top50_genes, ]

# 5. تحجيم البيانات (Z-score) عشان المقارنة تبقى أوضح في الـ Heatmap
mat_scaled <- t(scale(t(mat)))

# 6. عمل الـ Heatmap
pheatmap(mat_scaled,
         cluster_rows = TRUE,      # بنجمع الجينات المتشابهة
         cluster_cols = TRUE,      # بنجمع العينات المتشابهة
         show_rownames = TRUE,     # نعرض أسماء الجينات
         show_colnames = TRUE,     # نعرض أسماء العينات
         annotation_col = coldata["group"], # نلون العينات حسب المجموعة
         main = "Top 50 DEGs: CP17 vs WT",
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50), # تدرج لوني احترافي
         fontsize_row = 8,
         fontsize_col = 10,
         filename = "results/figures/heatmap_CP17_vs_WT.png", # حفظ مباشر
         width = 8, height = 10)

cat("Heatmap saved successfully in results/figures/heatmap_CP17_vs_WT.png\n")


# ============================================
# Script 04: Volcano Plots for other groups
# ============================================

# قائمة بالملفات اللي عايزين نعمللها تحليل
groups <- c("CP6", "P2", "P13")

for (grp in groups) {
  # قراءة الملف
  res_file <- paste0("results/tables/DEGs_", grp, "_vs_WT.csv")
  
  if (file.exists(res_file)) {
    res_data <- read.csv(res_file, row.names = 1)
    
    # عمل الرسمة
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
    
    # حفظ الرسمة
    ggsave(paste0("results/figures/volcano_", grp, "_vs_WT.png"), p_volcano, width = 10, height = 8, dpi = 300)
    
    # طباعة عدد الجينات المختلفة
    sig_count <- sum(!is.na(res_data$padj) & res_data$padj < 0.05 & abs(res_data$log2FoldChange) > 1)
    cat("Done for", grp, "- Significant DEGs:", sig_count, "\n")
  }
}



# ============================================
# Script 05: Extract Top Genes per Group
# ============================================

# قائمة بالمجموعات
groups <- c("CP17", "CP6", "P2", "P13")

for (grp in groups) {
  # قراءة الملف
  res_file <- paste0("results/tables/DEGs_", grp, "_vs_WT.csv")
  
  if (file.exists(res_file)) {
    res_data <- read.csv(res_file, row.names = 1)
    
    # فلترة الجينات المهمة
    sig_genes <- res_data[!is.na(res_data$padj) & res_data$padj < 0.05 & abs(res_data$log2FoldChange) > 2, ]
    
    # ترتيب حسب الـ log2FoldChange (لجذب أقوى التغييرات)
    sig_genes <- sig_genes[order(sig_genes$log2FoldChange, decreasing = TRUE), ]
    
    # حفظ الجدول
    output_file <- paste0("results/tables/Top_DEGs_", grp, ".csv")
    write.csv(sig_genes, output_file)
    
    cat("Done for", grp, "- Top DEGs saved to", output_file, "\n")
  }
}


#####################################################################
# قراءة وعرض الجينات المهمة
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