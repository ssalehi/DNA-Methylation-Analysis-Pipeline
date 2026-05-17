
--------------------------------------------------------------------------------
# 🧬 DNA Methylation Analysis Pipeline (Illumina 450K)

## 📌 Overview
This repository contains a complete bioinformatics pipeline developed in R for processing, normalizing, and analyzing DNA methylation data from Illumina HumanMethylation450 BeadChip arrays. The workflow heavily utilizes the Bioconductor `minfi` package to move from raw fluorescence intensities (IDAT files) to the discovery of Differentially Methylated Probes (DMPs) and advanced genomic visualizations.

This project was implemented as part of the DRD 2025 Bioinformatics workflow, focusing on comparing a Control group (CTRL) with a Disease group (DIS).

## 🚀 Pipeline Workflow & Key Features
The pipeline is structured into standard epigenetic data processing steps:

*   **Raw Data Importation:** Reading `.idat` files and metadata via `read.metharray.exp()`.
*   **Signal Extraction:** Extracting Red/Green fluorescence signals for targeted probe evaluation.
*   **Quality Control (QC):** 
    *   Evaluation of negative control probes.
    *   Calculation of detection p-values.
    *   **Custom Threshold:** Probes with a detection p > 0.01 were strictly filtered out to ensure high data reliability.
*   **Data Normalization:** 
    *   Applied the SWAN (Subset-quantile Within Array Normalization) method (`preprocessSWAN`).
    *   This effectively corrected technical biases between Type I and Type II probe chemistries, adjusting their distributions to be directly comparable.
*   **Dimensionality Reduction:** Conducted Principal Component Analysis (PCA) on normalized β-values to assess sample clustering based on biological groups (CTRL vs. DIS), sex, and batch effects.
*   **Differential Methylation Analysis:** 
    *   Applied the non-parametric Mann-Whitney U test across the first 50,000 probes.
    *   Applied Multiple Testing Corrections (Benjamini-Hochberg FDR and Bonferroni).
*   **Advanced Visualization (Automated Exports):** 
    *   **Density & PCA Plots:** Visualizing data distribution pre/post-normalization and PCA clustering.
    *   **Volcano Plots:** Visualizing significance vs. Δβ magnitude.
    *   **Manhattan Plots:** Mapping p-values across genomic positions using the `IlluminaHumanMethylation450kanno.ilmn12.hg19` database.
    *   **Heatmaps:** Unsupervised hierarchical clustering (using `gplots`) of the Top 100 differentially methylated CpGs.

## 🛠️ Tools & Libraries Used
*   **Language:** R
*   **Core Packages:** `minfi`, `minfiData`
*   **Annotation & Manifest:** `IlluminaHumanMethylation450kmanifest`, `IlluminaHumanMethylation450kanno.ilmn12.hg19`
*   **Visualization:** `gplots`, Base R Graphics

## 📁 Repository Structure
*   `project_codes.R`: The main R script containing the entire step-by-step pipeline. Running this script automatically generates the `Results/` directory.
*   `Results/`: An auto-generated directory containing all the final outputs of the pipeline, including:
    *   **Plots:** High-quality PNG files of all analyses (`Density_6Panels.png`, `PCA_Group.png`, `PCA_Sex.png`, `PCA_Batch.png`, `Volcano_Plot.png`, `Manhattan_Plot.png`, `Heatmap_Top100.png`).
    *   **Statistical Tables:** Exported CSV results (e.g., `Differential_Methylation_Results_Top50k.csv`).
    *   `🧬 DRD 2025 – Final Report.pdf`: Comprehensive project report containing biological interpretations, exploratory data analysis, and the final clustered Heatmap.

*(Note: Raw `.idat` files and intermediate `.RData` environments are excluded from this repository via `.gitignore` to maintain a clean structure and comply with size limits.)*

## 💡 Results
The pipeline successfully normalized the probe chemistries and identified distinct methylation signatures separating the CTRL and DIS groups, visualized perfectly in the final hierarchical clustering heatmap. For full graphical outputs and biological interpretation, please refer to the attached PDF report inside the `Results/` folder.

## 📥 Reproducibility & How to Run
Due to repository size limitations, the raw `.idat` files and temporary `.RData` environments are not included here.

To completely reproduce this pipeline on your local machine:
1. Create a folder named `Input_Data/` in the root directory.
2. Download the required raw Illumina 450K `.idat` files (GEO accessions starting with GSM5319592, etc.) and place them inside the `Input_Data/` folder along with the metadata file.
3. Open `project_codes.R` and run it. The code is already configured to automatically read from `./Input_Data/`, process the files seamlessly, and output all figures and tables directly into a newly created `Results/` folder.

--------------------------------------------------------------------------------