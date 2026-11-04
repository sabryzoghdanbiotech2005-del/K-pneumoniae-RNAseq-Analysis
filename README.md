# Transcriptomic Analysis of Antibiotic Resistance in *Klebsiella pneumoniae*

## 📌 Project Overview
This project analyzes publicly available RNA-seq data (GSE237490) to investigate the transcriptomic response of *K. pneumoniae* to different stress conditions (Antibiotics/Phage). The goal is to identify key genes and pathways associated with antimicrobial resistance.

## 📊 Key Findings
- **Differential Expression:** Dual-stress groups (CP17 & CP6) showed massive transcriptomic changes (1,832 and 840 DEGs respectively).
- **Efflux Pump:** The gene `VK055_0479` was consistently upregulated across resistant strains (CP17, CP6, P13), suggesting a shared resistance mechanism.
- **Capsule Biosynthesis:** Galactose metabolism genes (`VK055_1763-1766`) were highly upregulated in CP6.

## 🛠️ Tools & Technologies
- **Language:** R
- **Packages:** DESeq2, EnhancedVolcano, pheatmap, ggplot2, knitr.
- **Environment:** RStudio, Linux (Ubuntu).

## 📂 Repository Structure
- `data/`: Raw counts and metadata.
- `scripts/`: R scripts for the analysis pipeline.
- `results/`: Generated figures (PCA, Volcano, Heatmap) and DEG tables.
- `report_final.pdf`: Comprehensive PDF report of the analysis.

## 🚀 How to Reproduce
1. Clone the repo: `git clone https://github.com/sabryzoghdanbiotech2005-del/K-pneumoniae-RNAseq-Analysis.git`
2. Open the `.Rmd` file in RStudio.
3. Knit to PDF.

---
*Author: Sabry Ali Zoghdan*
