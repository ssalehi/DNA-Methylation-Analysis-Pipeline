
# Title: "🧬 DRD 2025 – Final Report"
# date: "2025-06-27"

## **Group Members:**

#| Name                | Student ID   |
#|---------------------|--------------|
#| Sareh Salehi        | 0001164160   |
#| Azam Bakhshandeh    | 0001164039   |
#| Shaghayegh Shirani  | 0001164225   |
#| Delnia Khezragha    | 0001167984   |
#| Maryam Kazemi       | 0001168170   |
#| Milad Arvand        | 0001164178   |

#--------------------------------------- Step 0

## ⚙️ Step 0 – Initialization (Environment Setup)

# Create Results directory if it doesn't exist
if (!dir.exists("Results")) {
  dir.create("Results")
}

# Clear the environment
rm(list = ls())

# Set working directory
# setwd("...")

# Install required packages (if not already installed)
#install.packages("rmarkdown")
#install.packages("knitr")
#install.packages("BiocManager")
#BiocManager::install("minfi")
#BiocManager::install("minfiData")

# Load libraries
library(rmarkdown)
library(knitr)
library(minfi)
library(minfiData)

# (Optional) View package documentation
# vignette("minfi")

# Define input directory
data_directory <- "./Input_Data/"
print(data_directory)

# List files in the input folder
list.files("./Input_Data")

# Read the sample sheet
SampleSheet <- read.csv("./Input_Data/SampleSheet_Report_II.csv", header = TRUE)

head(SampleSheet)

#--------------------------------------- Step 1

## 🧬 Step 1 – Importing Raw IDAT Data

# Load raw IDAT data using the sample sheet
RGset <- read.metharray.exp(base = data_directory)

# Save object for later steps
save(RGset,file="RGset.RData")

# Check the object class
class(RGset)

#--------------------------------------- Step 2

## 🔴🟢 Step 2 – Extracting Red and Green Fluorescence Signals

# Extract fluorescence signals from RGset
Red <- data.frame(getRed(RGset))
Green <- data.frame(getGreen(RGset))

# Check structure of the data
dim(Red); head(Red)
dim(Green); head(Green)

# Extract intensities for the assigned address (Group 2)
Red[rownames(Red) == "71773431", ]
Green[rownames(Green) == "71773431", ]

#--------------------------------------- Step 3

## 🔬 Step 3 – Extracting Signal for the Assigned Probe Address

# Build a data frame to extract signal intensities for the probe ID
df <- data.frame(
  Sample = colnames(Red),
  Red_fluor = as.numeric(Red["71773431", ]),
  Green_fluor = as.numeric(Green["71773431", ]),
  Type = rep(NA, ncol(Red)),     # Not found in manifest
  Color = rep(NA, ncol(Red))     # Not found in manifest
)

# Display the table
print(df)

#--------------------------------------- Step 3--Optional


### 🧪 Step 3 (Optional) – Probe Address Investigation

manifest <- getManifest(RGset)
probe_info <- getProbeInfo(manifest)
head(probe_info)

probe_info[probe_info$AddressA == "71773431" | probe_info$AddressB == "71773431", ]
probe_info[probe_info$AddressA == 71773431 | probe_info$AddressB == 71773431, ]

annotation <- getAnnotation(RGset)
annotation[annotation$AddressA_ID == "71773431" | annotation$AddressB_ID == "71773431", ]
annotation[annotation$AddressA_ID == 71773431 | annotation$AddressB_ID == 71773431, ]

nrow(getProbeInfo(RGset, type = "I"))
nrow(getProbeInfo(RGset, type = "II"))
nrow(getProbeInfo(RGset, type = "Control"))

df_I <- data.frame(getProbeInfo(RGset, type = "I"))
df_II <- data.frame(getProbeInfo(RGset, type = "II"))
df_control <- data.frame(getProbeInfo(RGset, type = "Control"))

df_I[df_I$AddressA == 71773431 | df_I$AddressB == 71773431, ]
df_II[df_II$AddressA == 71773431 | df_II$AddressB == 71773431, ]
df_control[df_control$AddressA == 71773431 | df_control$AddressB == 71773431, ]

df_I[df_I$AddressA == "71773431" | df_I$AddressB == "71773431", ]
df_II[df_II$AddressA == "71773431" | df_II$AddressB == "71773431", ]
df_control[df_control$AddressA == "71773431" | df_control$AddressB == "71773431", ]

#--------------------------------------- Step 4

## 🧬 Step 4 – Creating the Initial `MethylSet` Object (`MSet.raw`)

# Create the raw MethylSet object from RGset
MSet.raw <- preprocessRaw(RGset)

# Save the MSet object for future use
save(MSet.raw, file = "MSet_raw.RData")

# Extract methylated and unmethylated signals as matrices
Meth <- as.matrix(getMeth(MSet.raw))
Unmeth <- as.matrix(getUnmeth(MSet.raw))

# (Optional) Inspect structure and some rows of the matrices
str(Meth); head(Meth)
str(Unmeth); head(Unmeth)

#--------------------------------------- Step 5

## 🧬 Step 5 – Quality Control Assessment

qc <- getQC(MSet.raw)           # Extract QC metrics from raw MethylSet
plotQC(qc)                      # Visualize QC scores for all samples

controlStripPlot(RGset, controls = "NEGATIVE")

# Step 1: Calculate detection p-values
detP <- detectionP(RGset)  # This step might take a few minutes
save(detP, file = "detP.RData")  # Save the result for future use

# Step 2: Load detection p-values (if previously saved)
load("detP.RData")  # Alternatively: detP <- readRDS("detP.rds")

# Step 3: Count failed probes per sample (p-value > 0.01)
threshold <- 0.01  # Specific threshold used for Group 2
failedPositions <- colSums(detP > threshold)

# Step 4: Summarize the result
table(failedPositions)
summary(failedPositions)

# Step 5: Create result table
resultTable <- data.frame(
  Sample = names(failedPositions),
  `n° Failed positions` = as.vector(failedPositions),
  row.names = names(failedPositions)
)

resultTable  # Show the final table

#--------------------------------------- Step 6

## 🧬 Step 6 – B-Value and M-Value Analysis

# Load pre-processed methylation data
load("MSet_raw.RData")

# Extract beta and M values
beta <- getBeta(MSet.raw)
M <- getM(MSet.raw)

# Inspect distributions
colnames(beta)
summary(beta)
summary(M)


# 📥 Read sample metadata
SampleSheet <- read.csv("./Input_Data/SampleSheet_Report_II.csv", 
                        colClasses = c("Sentrix_ID" = "character"), 
                        header = TRUE)

# 🧬 Create unique sample names for matching
SampleSheet$Sample_Name <- paste(SampleSheet$SampleID,
                                 SampleSheet$Sentrix_ID,
                                 SampleSheet$Sentrix_Position,
                                 sep = "_")

# 🔄 Reorder metadata to match beta matrix columns
SampleSheet_ordered <- SampleSheet[match(colnames(beta), SampleSheet$Sample_Name), ]
group <- SampleSheet_ordered$Group  # Extract group labels
sum(is.na(group))

# Split the samples for β
beta_CTRL <- beta[, group == "CTRL"]
beta_DIS  <- beta[, group == "DIS"]

# Calculate the average β per CpG for each group
mean_beta_CTRL <- apply(beta_CTRL, 1, mean, na.rm = TRUE)
mean_beta_DIS  <- apply(beta_DIS, 1, mean, na.rm = TRUE)
mean_beta_CTRL <- na.omit(mean_beta_CTRL)
mean_beta_DIS <- na.omit(mean_beta_DIS)

# Split the samples for M
M_CTRL <- M[, group == "CTRL"]
M_DIS  <- M[, group == "DIS"]

# Calculate the average M per CpG for each group
mean_M_CTRL <- apply(M_CTRL, 1, mean, na.rm = TRUE)
mean_M_DIS  <- apply(M_DIS, 1, mean, na.rm = TRUE)
mean_M_CTRL <- na.omit(mean_M_CTRL)
mean_M_DIS <- na.omit(mean_M_DIS)

cat("▶ Mean Beta values (CTRL group):\n")
head(mean_beta_CTRL)

cat("\n▶ Mean Beta values (DIS group):\n")
head(mean_beta_DIS)

cat("\n▶ Mean M values (CTRL group):\n")
head(mean_M_CTRL)

cat("\n▶ Mean M values (DIS group):\n")
head(mean_M_DIS)


# Density plot for mean Beta values - CTRL group
d_mean_beta_CTRL <- density(mean_beta_CTRL)
plot(d_mean_beta_CTRL ,main="Density of Beta Values for group CTRL",col="orange",lwd = 2)

# Density plot for mean M values - CTRL group
d_mean_M_CTRL <- density(mean_M_CTRL)
plot(d_mean_M_CTRL,main="Density of M Values for group CTRL",col="purple",lwd = 2)

# Density plot for mean Beta values - DIS group
d_mean_beta_DIS <- density(mean_beta_DIS)
plot(d_mean_beta_DIS ,main="Density of Beta Values for group DIS",col="orange", lwd = 2)

# Density plot for mean M values - DIS group
d_mean_M_DIS <- density(mean_M_DIS)
plot(d_mean_M_DIS, main="Density of M Values for group DIS",col="purple", lwd = 2)


# Comparison of mean Beta values across groups
plot(density(mean_beta_CTRL), col = "orange", lwd = 2,
     main = "Comparison of Mean Beta Values (CTRL vs DIS)",
     xlab = "Mean Beta Value")
lines(density(mean_beta_DIS), col = "blue", lwd = 2)
legend("topright", legend = c("CTRL", "DIS"),
       col = c("orange", "blue"), lwd = 2)

# Comparison of mean M values across groups
plot(density(mean_M_CTRL), col = "orange", lwd = 2,
     main = "Comparison of Mean M Values (CTRL vs DIS)",
     xlab = "Mean M Value")
lines(density(mean_M_DIS), col = "blue", lwd = 2)
legend("topright", legend = c("CTRL", "DIS"),
       col = c("orange", "blue"), lwd = 2)

#--------------------------------------- Step 7

## 🧬 Step 7 – Normalize the Data using `preprocessSWAN`

MSet.raw <- preprocessRaw(RGset)
MSet.norm <- preprocessSWAN(RGset, MSet.raw)

# Extract Beta values
beta_raw <- getBeta(MSet.raw)
beta_norm <- getBeta(MSet.norm)

# Determine probe types
probeType <- getProbeType(RGset)
typeI <- probeType == "I"
typeII <- probeType == "II"


# Raw
mean_raw_I <- rowMeans(beta_raw[typeI, ], na.rm = TRUE)
mean_raw_II <- rowMeans(beta_raw[typeII, ], na.rm = TRUE)
sd_raw_I <- apply(beta_raw[typeI, ], 1, sd, na.rm = TRUE)
sd_raw_II <- apply(beta_raw[typeII, ], 1, sd, na.rm = TRUE)

# Normalized
mean_norm_I <- rowMeans(beta_norm[typeI, ], na.rm = TRUE)
mean_norm_II <- rowMeans(beta_norm[typeII, ], na.rm = TRUE)
sd_norm_I <- apply(beta_norm[typeI, ], 1, sd, na.rm = TRUE)
sd_norm_II <- apply(beta_norm[typeII, ], 1, sd, na.rm = TRUE)


# Label setup
short_labels_raw <- substr(colnames(beta_raw), 1, 8)
short_labels_norm <- substr(colnames(beta_norm), 1, 8)

# 6 plots layout
png("Results/Density_6Panels.png", width = 1200, height = 800, res = 120)
par(mfrow = c(2, 3))

# 1. Mean of RAW
plot(density(mean_raw_I, na.rm = TRUE), col = "black", 
     main = "Mean of raw Beta values", xlab = "Beta mean")
lines(density(mean_raw_II, na.rm = TRUE), col = "blue")
legend("topright", legend = c("Type I", "Type II"), col = c("black", "blue"), lty = 1)

# 2. SD of RAW
plot(density(sd_raw_I, na.rm = TRUE), col = "black", 
     main = "SD of raw Beta values", xlab = "Beta SD")
lines(density(sd_raw_II, na.rm = TRUE), col = "blue")
legend("topright", legend = c("Type I", "Type II"), col = c("black", "blue"), lty = 1)

# 3. Boxplot of RAW
boxplot(beta_raw, main = "Raw Beta values", col = c("red", "green", "blue"),
        names = short_labels_raw)

# 4. Mean of NORM
plot(density(mean_norm_I, na.rm = TRUE), col = "black", 
     main = "Mean of normalized Beta values", xlab = "Beta mean")
lines(density(mean_norm_II, na.rm = TRUE), col = "blue")
legend("topright", legend = c("Type I", "Type II"), col = c("black", "blue"), lty = 1)

# 5. SD of NORM
plot(density(sd_norm_I, na.rm = TRUE), col = "black", 
     main = "SD of normalized Beta values", xlab = "Beta SD")
lines(density(sd_norm_II, na.rm = TRUE), col = "blue")
legend("topright", legend = c("Type I", "Type II"), col = c("black", "blue"), lty = 1)

# 6. Boxplot of NORM
boxplot(beta_norm, main = "Normalized Beta values", col = c("red", "green", "blue"),
        names = short_labels_norm)

dev.off()

#--------------------------------------- Step 7--ptional

group <- c("WT", "WT", "WT", "MUT", "MUT", "MUT", "WT", "MUT")
pData(RGset)$Group <- group

group_colors <- ifelse(group == "WT", "darkgreen", "darkred")
par(mfrow = c(1, 2))

# Raw
boxplot(beta_raw, col = group_colors,
        main = "Raw Beta values by group",
        las = 2, names = substr(colnames(beta_raw), 1, 8))

# Normalized
boxplot(beta_norm, col = group_colors,
        main = "Normalized Beta values by group",
        las = 2, names = substr(colnames(beta_norm), 1, 8))

#--------------------------------------- Step 8

## 🧬 Step 8 – PCA Analysis on Normalized Beta Values

meta <- read.csv("./Input_Data/SampleSheet_Report_II.csv", stringsAsFactors = FALSE)
meta$CombinedID <- paste(meta$SampleID, meta$Sentrix_ID, meta$Sentrix_Position, sep = "_")

meta <- meta[match(colnames(beta_norm), meta$CombinedID), ]
all(meta$CombinedID == colnames(beta_norm))

# Remove probes that have full NA
beta_pca <- beta_norm[complete.cases(beta_norm), ]
beta_t <- t(beta_pca)

# Run PCA
pca_result <- prcomp(beta_t, center = TRUE, scale. = TRUE)

# Colorization based on CTRL / DIS group
group_colors <- ifelse(meta$Group == "CTRL", "darkgreen", "firebrick")

# Plot PCA by group
png("Results/PCA_Group.png", width = 800, height = 600, res = 120)
plot(pca_result$x[,1:2],
     col = group_colors,
     pch = 19,
     xlab = "PC1", ylab = "PC2",
     main = "PCA – Colored by Group")
grid(col = "gray85", lty = "dotted")
legend("topright", legend = c("CTRL", "DIS"),
       col = c("darkgreen", "firebrick"), pch = 19)

dev.off()

# Color coding based on gender
sex_colors <- ifelse(meta$Sex == "Male", "navy", "orange")

# Plot PCA based on gender
png("Results/PCA_Sex.png", width = 800, height = 600, res = 120)
plot(pca_result$x[,1:2],
     col = sex_colors,
     pch = 19,
     xlab = "PC1", ylab = "PC2",
     main = "PCA – Colored by Sex")
grid(col = "gray85", lty = "dotted")
legend("topright", legend = c("Male", "Female"),
       col = c("navy", "orange"), pch = 19)

dev.off() 

# Define colors based on Batch (Sentrix_ID)
png("Results/PCA_Batch.png", width = 800, height = 600, res = 120)
batch_colors <- as.factor(meta$Sentrix_ID)

plot(pca_result$x[,1:2],
     col = batch_colors,
     pch = 19,
     xlab = "PC1", ylab = "PC2",
     main = "PCA – Colored by Batch")
grid(col = "gray85", lty = "dotted")
legend("topright", legend = levels(batch_colors),
       col = 1:length(levels(batch_colors)), pch = 19)
dev.off()

#--------------------------------------- Step 9

## 🧪 Step 9 – Identification of Differentially Methylated Probes

beta_norm
first50k_beta_norm <- beta_norm[1:50000, ]

My_mannwhitney_function <- function(x) {
  wilcox <- wilcox.test(x ~ meta$Group)
  return(wilcox$p.value)
}

pValues_wilcox_first50k <- apply(first50k_beta_norm, 1, My_mannwhitney_function)

final_wilcox_first50k <- data.frame(first50k_beta_norm, pValues_wilcox_first50k)
final_wilcox_first50k <- final_wilcox_first50k[order(final_wilcox_first50k$pValues_wilcox_first50k), ]
final_wilcox_first50k_0.01 <- final_wilcox_first50k[final_wilcox_first50k$pValues_wilcox_first50k <= 0.01, ]
final_wilcox_first50k_0.05 <- final_wilcox_first50k[final_wilcox_first50k$pValues_wilcox_first50k <= 0.05, ]


dim(final_wilcox_first50k_0.01)
dim(final_wilcox_first50k_0.05)

#--------------------------------------- Step 10

## 🧬 Step 10 – Multiple Testing Correction


# Output from previous Wilcoxon test on first 50,000 probes
corrected_pValues_BH <- p.adjust(final_wilcox_first50k$pValues_wilcox_first50k, method = "BH")
corrected_pValues_Bonf <- p.adjust(final_wilcox_first50k$pValues_wilcox_first50k, method = "bonferroni")

final_wilcox_first50k_corrected <- data.frame(final_wilcox_first50k, corrected_pValues_BH, corrected_pValues_Bonf)

# test based on P-Value (nominal)
dim(final_wilcox_first50k_corrected[final_wilcox_first50k_corrected$pValues_wilcox_first50k <= 0.05, ])

# after BH
dim(final_wilcox_first50k_corrected[final_wilcox_first50k_corrected$corrected_pValues_BH <= 0.05, ])

# after Bonferroni
dim(final_wilcox_first50k_corrected[final_wilcox_first50k_corrected$corrected_pValues_Bonf <= 0.05, ])


#--------------------------------------- Step 11

## 🔬 Step 11 – Visualization of Differential Methylation Results

# Extract the beta matrix for the first 50,000 probes
beta_first50k <- final_wilcox_first50k_corrected[, 1:8]


# Construct submatrices for each group
beta_first50k_CTRL <- beta_first50k[, meta$Group == "CTRL"]
beta_first50k_DIS  <- beta_first50k[, meta$Group == "DIS"]

# Calculate the average beta in each group
mean_beta_CTRL <- apply(beta_first50k_CTRL, 1, mean)
mean_beta_DIS  <- apply(beta_first50k_DIS, 1, mean)

# Calculate Δβ
delta_beta <- mean_beta_DIS - mean_beta_CTRL
delta_beta

# Create a data frame for Volcano Plot
toVolcano <- data.frame(
  delta_beta = delta_beta,
  negLog10_p = -log10(final_wilcox_first50k_corrected$pValues_wilcox_first50k)
)

# Draw the initial chart
plot(toVolcano$delta_beta, toVolcano$negLog10_p, pch=16, cex=0.5,
     xlab = "Δβ (DIS - CTRL)", ylab = "-log10(p-value)", main = "Volcano Plot")

toHighlight <- toVolcano[abs(toVolcano$delta_beta) > 0.05 & toVolcano$negLog10_p > -log10(0.05), ]
nrow(toHighlight)

# Redraw with coloring
plot(toVolcano$delta_beta, toVolcano$negLog10_p, pch=16, cex=0.5,
     xlab="Δβ (DIS - CTRL)", ylab="-log10(p-value)", main="Volcano Plot")

# Add threshold line
abline(h = -log10(0.01), col = "red", lty = 2)

# Add important points with different colors
points(toHighlight$delta_beta, toHighlight$negLog10_p, pch=16, cex=0.7, col="orange")



if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("IlluminaHumanMethylation450kmanifest")
BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19")

library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(minfi)

annotation_full <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
colnames(annotation_full)

manifest_clean <- data.frame(
  Name = annotation_full$Name,
  CHR = annotation_full$chr,
  MAPINFO = annotation_full$pos
)

save(manifest_clean, file = "Illumina450Manifest_clean.RData")

# Load genomic annotation
load("Illumina450Manifest_clean.RData") 

# Create a data frame containing p-value and probe IDs
toManhattan <- data.frame(
  Name = rownames(final_wilcox_first50k_corrected),
  pval = final_wilcox_first50k_corrected$pValues_wilcox_first50k
)

# Merge with annotation to capture chromosome and position information
toManhattan_annotated <- merge(toManhattan, manifest_clean[, c("Name", "CHR", "MAPINFO")],
                               by = "Name")
# Remove incomplete values​(NA)
toManhattan_annotated <- na.omit(toManhattan_annotated)

# Drawing Manhattan Plot
png("Results/Manhattan_Plot.png", width = 1000, height = 600, res = 120)
plot(toManhattan_annotated$MAPINFO, -log10(toManhattan_annotated$pval),
     pch = 16, cex = 0.5, col = "darkblue",
     xlab = "Genomic Position", ylab = "-log10(p-value)",
     main = "Manhattan Plot (first 50k probes)",
     ylim = c(0, 5))

abline(h = -log10(0.01), col = "red", lty = 2, lwd = 2)
dev.off() 

#--------------------------------------- Step 12

## 🔥 Step 12 – Heatmap of Top 100 Differentially Methylated Probes

# Load required package
library(gplots)

# Extract beta-values of top 100 CpGs (ranked by significance)
colnames(final_wilcox_first50k_corrected)
input_heatmap <- as.matrix(final_wilcox_first50k_corrected[1:100, 1:8])

# Load phenotype data to assign group labels
pheno <- read.csv("./Input_Data/SampleSheet_Report_II.csv", header = TRUE, stringsAsFactors = TRUE)
pheno$Group

# Assign color based on group
colorbar <- ifelse(pheno$Group == "CTRL", "green", "orange")

# Generate heatmap with clustering
heatmap.2(input_heatmap,
          col = terrain.colors(100),
          Rowv = TRUE,
          Colv = TRUE,
          dendrogram = "both",
          trace = "none",
          ColSideColors = colorbar,
          key = TRUE,
          scale = "none",
          cexRow = 0.5,
          main = "Heatmap of Top 100 Differentially Methylated Probes")





# Export Differential Methylation Results
write.csv(final_wilcox_first50k_corrected, 
          file = "Results/Differential_Methylation_Results_Top50k.csv", 
          row.names = TRUE)


# Export Volcano Plot
png("Results/Volcano_Plot.png", width = 800, height = 600, res = 120)

plot(toVolcano$delta_beta, toVolcano$negLog10_p, pch=16, cex=0.5, 
     xlab="Δβ (DIS - CTRL)", ylab="-log10(p-value)", main="Volcano Plot")
abline(h = -log10(0.01), col = "red", lty = 2)
points(toHighlight$delta_beta, toHighlight$negLog10_p, pch=16, cex=0.7, col="orange")

dev.off()



# Export Heatmap
png("Results/Heatmap_Top100.png", width = 800, height = 800, res = 120)

heatmap.2(input_heatmap,
          col = terrain.colors(100),
          Rowv = TRUE,
          Colv = TRUE,
          dendrogram = "both",
          trace = "none",
          ColSideColors = colorbar,
          key = TRUE,
          scale = "none",
          cexRow = 0.5,
          main = "Heatmap of Top 100 Differentially Methylated Probes")

dev.off()
