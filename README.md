
# Transcriptomic Analysis of Antibiotic Resistance in *Klebsiella pneumoniae*

## 📌 Project Overview
This project analyzes publicly available RNA-seq data (**GSE237490**) to investigate the transcriptomic response of *K. pneumoniae* to different stress conditions (Antibiotics/Phage). The goal is to identify key genes and pathways associated with antimicrobial resistance and stress adaptation.

## 🛠️ Tools & Technologies
![R](https://img.shields.io/badge/R-276DC3?style=for-the-badge&logo=r&logoColor=white)
![RStudio](https://img.shields.io/badge/RStudio-75AADB?style=for-the-badge&logo=rstudio&logoColor=white)
![Linux](https://img.shields.io/badge/Linux-FCC624?style=for-the-badge&logo=linux&logoColor=black)

- **Language:** R
- **Packages:** DESeq2, EnhancedVolcano, pheatmap, ggplot2, knitr.
- **Environment:** RStudio, Linux (Ubuntu).

## 🔄 Workflow
1. **Data Retrieval:** Downloaded raw counts and metadata for *K. pneumoniae* RNA-seq dataset (GSE237490).
2. **Quality Control & Preprocessing:** Filtered low-count genes and normalized data using DESeq2.
3. **Differential Expression Analysis:** Identified Differentially Expressed Genes (DEGs) between control and stress conditions (CP17, CP6, P13).
4. **Functional Analysis:** Investigated genes related to efflux pumps and capsule biosynthesis.
5. **Visualization:** Generated PCA plots, Volcano plots, and Heatmaps.

## 📊 Key Findings
- **Differential Expression:** Dual-stress groups (CP17 & CP6) showed massive transcriptomic changes (**1,832** and **840** DEGs respectively).
- **Efflux Pump Mechanism:** The gene `VK055_0479` was consistently upregulated across resistant strains (CP17, CP6, P13), suggesting a shared resistance mechanism.
- **Capsule Biosynthesis:** Galactose metabolism genes (`VK055_1763-1766`) were highly upregulated in CP6.

## 📂 Repository Structure
- `data/` : Raw counts and metadata.
- `scripts/` : R scripts (`.Rmd`) for the analysis pipeline.
- `results/` : Generated figures (PCA, Volcano, Heatmap) and DEG tables.
- `report_final.pdf` : Comprehensive PDF report of the analysis.

## 🚀 How to Reproduce
1. Clone the repository:
   ```bash
   git clone https://github.com/sabryzoghdanbiotech2005-del/K-pneumoniae-RNAseq-Analysis.git
Open the .Rmd file in RStudio.

Knit to PDF to generate the full report.

👤 Author
Sabry Ali Zoghdan

GitHub: sabryzoghdanbiotech2005-del

Email: sabryzoghdan.biotech2005@gmail.com
