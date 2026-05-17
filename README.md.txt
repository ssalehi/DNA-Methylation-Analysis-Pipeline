# 🧬 DNA Methylation Analysis Pipeline (Illumina 450K)

## 📌 Overview
This repository contains a complete bioinformatics pipeline developed in **R** for processing, normalizing, and analyzing DNA methylation data from **Illumina HumanMethylation450 BeadChip** arrays. The workflow heavily utilizes the Bioconductor `minfi` package to move from raw fluorescence intensities (IDAT files) to the discovery of Differentially Methylated Probes (DMPs) and advanced genomic visualizations.

This project was implemented as part of the DRD 2025 Bioinformatics workflow, focusing on comparing a Control group (CTRL) with a Disease group (DIS).

## 🚀 Pipeline Workflow & Key Features
The pipeline is structured into standard epigenetic data processing steps:

1. **Raw Data Importation:** Reading `.idat` files and metadata via `read.metharray.exp()`.
2. **Signal Extraction:** Extracting Red/Green fluorescence signals for targeted probe evaluation.
3. **Quality Control (QC):** 
   - Evaluation of negative control probes.
   - Calculation of detection p-values.
   - **Custom Threshold:** Probes with a detection $p > 0.01$ were strictly filtered out to ensure high data reliability.
4. **Data Normalization:** 
   - Applied the **SWAN (Subset-quantile Within Array Normalization)** method (`preprocessSWAN`).
   - This effectively corrected technical biases between Type I and Type II probe chemistries, adjusting their distributions to be directly comparable.
5. **Dimensionality Reduction:** Conducted **Principal Component Analysis (PCA)** on normalized $\beta$-values to assess sample clustering based on biological groups (CTRL vs. DIS) and batch effects.
6. **Differential Methylation Analysis:** 
   - Applied the non-parametric **Mann-Whitney U test** across the first 50,000 probes.
   - Applied Multiple Testing Corrections (Benjamini-Hochberg FDR and Bonferroni).
7. **Advanced Visualization:** 
   - **Volcano Plots:** Visualizing significance vs. $\Delta\beta$ magnitude.
   - **Manhattan Plots:** Mapping p-values across genomic positions using the `IlluminaHumanMethylation450kanno.ilmn12.hg19` database.
   - **Heatmaps:** Unsupervised hierarchical clustering (using `gplots`) of the Top 100 differentially methylated CpGs.

## 🛠️ Tools & Libraries Used
* **Language:** R
* **Core Packages:** `minfi`, `minfiData`
* **Annotation & Manifest:** `IlluminaHumanMethylation450kmanifest`, `IlluminaHumanMethylation450kanno.ilmn12.hg19`
* **Visualization:** `gplots`, Base R Graphics

## 📁 Repository Structure
* `project_codes.R`: The main R script containing the entire step-by-step pipeline.
* `DRD_2025_Final_Report.pdf`: Comprehensive project report containing biological interpretations, exploratory data analysis, PCA plots, Volcano/Manhattan plots, and the final clustered Heatmap.
*(Note: Raw `.idat` and `.RData` environment files are excluded from this repository due to size limits.)*

## 💡 Results
The pipeline successfully normalized the probe chemistries and identified distinct methylation signatures separating the CTRL and DIS groups, visualized perfectly in the final hierarchical clustering heatmap. For full graphical outputs and biological interpretation, please refer to the attached PDF report.

--------------------------------------------------------------------------------