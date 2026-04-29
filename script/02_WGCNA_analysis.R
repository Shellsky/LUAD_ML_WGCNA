############################################################
# TCGA-LUAD WGCNA Analysis Pipeline
#
# Description:
# This script performs weighted gene co-expression network
# analysis (WGCNA) using TCGA-LUAD gene expression data.
#
# Workflow:
# 1. Load LUAD expression and clinical data
# 2. Preprocess expression matrix and clinical phenotype
# 3. Select highly variable genes
# 4. Detect and remove outlier samples
# 5. Construct WGCNA co-expression modules
# 6. Correlate module eigengenes with survival status
# 7. Identify module membership and gene significance
# 8. Export module networks for Cytoscape
############################################################


############################################################
# 1. Load required packages
############################################################

library(tidyverse)
library(org.Hs.eg.db)
library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)
library(tibble)
library(WGCNA)
library(stringr)


############################################################
# 2. Clear environment and set working directory
############################################################

rm(list = ls())

# Set working directory for LUAD analysis
setwd("LUAD")


############################################################
# 3. Load LUAD gene expression data
############################################################

# Read TCGA-LUAD gene expression matrix
exp <- read_csv(".csv")

# Convert tibble to data frame
exp <- as.data.frame(exp)

# Use the first column as gene names
rownames(exp) <- as.character(exp[, 1])

# Remove the first column after assigning row names
exp <- exp[, -1]

# Store original sample names
get_exp_name <- colnames(exp)


############################################################
# 4. Load LUAD clinical data
############################################################

# Read LUAD clinical data
lc <- read_csv("LUAD_Clinical.csv")

# Convert TCGA sample barcode to patient-level barcode
colnames(exp) <- substr(colnames(exp), 1, 12)

exp_name <- substr(colnames(exp), 1, 12)


############################################################
# 5. Load cleaned LUAD clinical data
############################################################


# Select relevant clinical variables
lc <- lc %>%
  select(
    vital_status,
    ajcc_pathologic_m,
    ajcc_pathologic_n,
    ajcc_pathologic_t,
    age_at_index,
    cigarettes_per_day,
    bcr_patient_barcode
  )

# Remove missing values
data <- lc %>% 
  na.omit()

# Filter invalid TNM staging values and keep survival status
lc <- lc %>%
  filter(
    ajcc_pathologic_m != "MX",
    ajcc_pathologic_t != "TX"
  ) %>%
  dplyr::select(
    bcr_patient_barcode,
    vital_status
  )

# Remove missing records
lc <- lc %>%
  na.omit()


############################################################
# 6. Remove duplicated patient barcodes
############################################################

duplicated_positions <- duplicated(exp_name)

if (any(duplicated_positions)) {
  duplicate_indices <- which(duplicated_positions)
  cat("Duplicated sample positions:", paste(duplicate_indices, collapse = ", "), "\n")
} else {
  cat("No duplicated sample positions found.\n")
}

# Extract duplicated samples for checking
test <- exp[, duplicated_positions]
ch <- data.frame(colnames(test))

# Remove duplicated samples
exp <- exp[, !duplicated_positions]

# Store filtered expression sample names
get_exp_name <- colnames(exp)

get_exp_name <- data.frame(
  bcr_patient_barcode = get_exp_name
)


############################################################
# 7. Match expression data with clinical data
############################################################

# Merge expression sample names with clinical survival data
data_frame <- left_join(
  get_exp_name,
  lc,
  by = "bcr_patient_barcode"
)

# Remove unmatched or missing samples
data_frame <- data_frame %>%
  na.omit()


############################################################
# 8. Check sample type
############################################################

# TCGA sample code:
# 01-09 = tumor samples
# 10-19 = normal samples
check_sample <- ifelse(
  substring(colnames(exp), 14, 15) == "11",
  "normal",
  "cancer"
)

table(check_sample)


############################################################
# 9. Select top 5,000 highly variable genes
############################################################

# Select genes with the highest median absolute deviation (MAD)
# Rows = genes, columns = samples
# Transpose matrix so that rows = samples and columns = genes
data.mat <- t(
  exp[
    order(apply(exp, 1, mad), decreasing = TRUE)[1:5000],
  ]
)

dim(data.mat)


############################################################
# 10. Check good samples and genes
############################################################

gsg <- goodSamplesGenes(data.mat, verbose = 3)

# Remove poor-quality genes or samples if detected
if (!gsg$allOK) {
  
  if (sum(!gsg$goodGenes) > 0) {
    printFlush(
      paste(
        "Removing genes:",
        paste(names(data.mat)[!gsg$goodGenes], collapse = ",")
      )
    )
  }
  
  if (sum(!gsg$goodSamples) > 0) {
    printFlush(
      paste(
        "Removing samples:",
        paste(rownames(data.mat)[!gsg$goodSamples], collapse = ",")
      )
    )
  }
  
  data.mat <- data.mat[gsg$goodSamples, gsg$goodGenes]
}


############################################################
# 11. Sample clustering for outlier detection
############################################################

sampleTree <- hclust(
  dist(data.mat),
  method = "average"
)

sizeGrWindow(15, 12)

par(cex = 0.8)
par(mar = c(0, 5, 2, 0))

plot(
  sampleTree,
  main = "Sample clustering to detect outliers",
  sub = "",
  xlab = "",
  cex.lab = 1.5,
  cex.axis = 1.5,
  cex.main = 2
)

abline(h = 120000, col = "red")

sampleTree$height


############################################################
# 12. Remove outlier samples based on clustering
############################################################

# Cut the sample tree to remove outlier clusters
clusters <- cutree(sampleTree, h = 150000)

table(clusters)

# Keep samples in the major cluster
data.mat_filtered <- data.mat[clusters == 1, ]

# Re-cluster filtered samples
sampleTree3 <- hclust(
  dist(data.mat_filtered),
  method = "average"
)

sizeGrWindow(15, 12)

par(cex = 0.8)
par(mar = c(0, 5, 2, 0))

plot(
  sampleTree3,
  main = "Sample clustering after outlier removal",
  sub = "",
  xlab = "",
  cex.lab = 1.5,
  cex.axis = 1.5,
  cex.main = 2
)

# Replace original matrix with filtered matrix
data.mat <- data.mat_filtered


############################################################
# 13. Select soft-thresholding power
############################################################

# Enable multi-threading for WGCNA
allowWGCNAThreads()

# Define network type
type <- "unsigned"

# Candidate soft-thresholding powers
powers <- c(1:10, seq(from = 12, to = 30, by = 2))

# Estimate soft-thresholding power
sft <- pickSoftThreshold(
  data.mat,
  powerVector = powers,
  networkType = type,
  verbose = 3
)


############################################################
# 14. Plot scale-free topology and mean connectivity
############################################################

par(mfrow = c(1, 2))

cex1 <- 0.9

# Scale-free topology model fit
plot(
  sft$fitIndices[, 1],
  -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
  xlab = "Soft Threshold (power)",
  ylab = "Scale Free Topology Model Fit, signed R^2",
  type = "n",
  main = "Scale independence"
)

text(
  sft$fitIndices[, 1],
  -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
  labels = powers,
  cex = cex1,
  col = "red"
)

abline(h = 0.85, col = "red")

# Mean connectivity
plot(
  sft$fitIndices[, 1],
  sft$fitIndices[, 5],
  xlab = "Soft Threshold (power)",
  ylab = "Mean Connectivity",
  type = "n",
  main = "Mean connectivity"
)

text(
  sft$fitIndices[, 1],
  sft$fitIndices[, 5],
  labels = powers,
  cex = cex1,
  col = "red"
)

# Display estimated soft-thresholding power
print(sft$powerEstimate)


############################################################
# 15. Define soft-thresholding power
############################################################

# Use estimated power if available; otherwise use a default rule
power <- sft$powerEstimate

nSamples <- nrow(data.mat)

if (is.na(power)) {
  power <- ifelse(
    nSamples < 20,
    ifelse(type == "unsigned", 9, 18),
    ifelse(
      nSamples < 30,
      ifelse(type == "unsigned", 8, 16),
      ifelse(
        nSamples < 40,
        ifelse(type == "unsigned", 7, 14),
        ifelse(type == "unsigned", 6, 12)
      )
    )
  )
}

print(power)


############################################################
# 16. Construct weighted gene co-expression network
############################################################

net <- blockwiseModules(
  data.mat,
  power = power,
  maxBlockSize = 5000,
  TOMType = type,
  minModuleSize = 100,
  reassignThreshold = 0.25,
  mergeCutHeight = 0.05,
  numericLabels = TRUE,
  pamRespectsDendro = FALSE,
  saveTOMs = TRUE,
  verbose = 3
)

# Show number of genes in each module
table(net$colors)

# Convert numeric labels to color labels
moduleColors <- labels2colors(net$colors)


############################################################
# 17. Plot gene dendrogram and module colors
############################################################

plotDendroAndColors(
  net$dendrograms[[1]],
  moduleColors[net$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05
)


############################################################
# 18. Module eigengene network analysis
############################################################

# Extract module eigengenes
MEs_col <- net$MEs

# Convert numeric module names to color names
colnames(MEs_col) <- paste0(
  "ME",
  labels2colors(
    as.numeric(str_replace_all(colnames(MEs_col), "ME", ""))
  )
)

# Order module eigengenes
MEs_col <- orderMEs(MEs_col)

# Plot eigengene adjacency heatmap
plotEigengeneNetworks(
  MEs_col,
  "Eigengene adjacency heatmap",
  marDendro = c(3, 3, 2, 4),
  marHeatmap = c(3, 4, 2, 2),
  plotDendrograms = TRUE,
  xLabelsAngle = 90
)


############################################################
# 19. Prepare clinical trait matrix
############################################################

# Prepare phenotype table
test_datafrmae <- data.frame(
  ID = data_frame$bcr_patient_barcode,
  session = data_frame$vital_status
)

# Create design matrix for survival status
design <- model.matrix(~0 + test_datafrmae$session)

dimnames(design) <- list(test_datafrmae$ID)

# Match clinical trait data with WGCNA samples
design <- design[rownames(MEs_col), ]

# Rename phenotype columns
colnames(design) <- c("Alive", "Dead")


############################################################
# 20. Module-trait correlation analysis
############################################################

# Pearson correlation between module eigengenes and clinical traits
modTraitCor <- cor(MEs_col, design, use = "p")

modTraitP <- corPvalueStudent(
  modTraitCor,
  dim(test_datafrmae)[1]
)

# Robust bicorrelation analysis
modTraitCorP <- bicorAndPvalue(MEs_col, design)

modTraitCor <- modTraitCorP$bicor
modTraitP <- modTraitCorP$p

# Create text matrix for heatmap
textMatrix <- paste0(
  signif(modTraitCor, 2),
  "\n(",
  signif(modTraitP, 1),
  ")"
)

dim(textMatrix) <- dim(modTraitCor)

dev.off()

par(mar = c(2, 10, 2, 0))

# Plot module-trait relationship heatmap
labeledHeatmap(
  Matrix = modTraitCor,
  xLabels = colnames(design),
  yLabels = colnames(MEs_col),
  cex.lab = 0.5,
  ySymbols = colnames(MEs_col),
  colorLabels = TRUE,
  colors = blueWhiteRed(100),
  textMatrix = textMatrix,
  setStdMargins = FALSE,
  cex.text = 0.8,
  zlim = c(-1, 1),
  main = "Module-trait relationships"
)


############################################################
# 21. Extract genes from the key module
############################################################

# Select the module of interest
module <- "blue"

# Extract genes belonging to the selected module
moduleGenes <- names(net$colors)[which(moduleColors == module)]

print(moduleGenes)


############################################################
# 22. TOM similarity matrix and network heatmap
############################################################

# Load saved TOM matrix
load(net$TOMFiles[1], verbose = TRUE)

# Recalculate TOM if needed
TOM <- TOMsimilarityFromExpr(
  data.mat,
  power = power,
  networkType = type
)

TOM <- as.matrix(TOM)

# Convert TOM to dissimilarity matrix
dissTOM <- 1 - TOM

plotTOM <- dissTOM^7
diag(plotTOM) <- NA

TOMplot(
  plotTOM,
  net$dendrograms[[1]],
  moduleColors[net$blockGenes[[1]]],
  main = "Network heatmap plot, all genes"
)


############################################################
# 23. Recalculate module eigengenes
############################################################

nGenes <- ncol(data.mat)
nSamples <- nrow(data.mat)

Color <- labels2colors(net$colors)

MEs0 <- moduleEigengenes(data.mat, Color)$eigengenes
MEs <- orderMEs(MEs0)


############################################################
# 24. Calculate module membership and gene significance
############################################################

nSamples <- dim(data.mat)[1]

# Module membership: correlation between genes and module eigengenes
geneModuleMembership <- cor(
  data.mat,
  MEs_col,
  use = "p",
  method = "spearman"
)

MMPvalue <- corPvalueStudent(
  geneModuleMembership,
  nSamples
)

# Gene significance: correlation between genes and clinical traits
geneSignificanceCor <- cor(
  data.mat,
  design,
  use = "p",
  method = "spearman"
)

geneSignificanceP <- corPvalueStudent(
  geneSignificanceCor,
  nSamples
)


############################################################
# 25. Evaluate genes in the key module
############################################################

module <- "blue"
column <- paste0("ME", module)

moduleGenes <- names(net$colors)[which(moduleColors == module)]

# Keep only valid module genes present in correlation matrices
valid_genes <- moduleGenes[moduleGenes %in% rownames(geneModuleMembership)]

# Module membership values
MM <- abs(geneModuleMembership[, column])
MM <- MM[rownames(geneModuleMembership) %in% valid_genes]

# Gene significance values for the first phenotype column
GS <- abs(geneSignificanceCor[, 1])
GS <- GS[rownames(geneSignificanceCor) %in% valid_genes]


############################################################
# 26. Scatter plot: Module membership vs Gene significance
############################################################

dev.off()

verboseScatterplot(
  MM,
  GS,
  xlab = paste("Module Membership in", module, "module"),
  ylab = "Gene significance for survival status",
  main = "Module membership vs. gene significance",
  abline = TRUE,
  pch = 21,
  cex.main = 1.2,
  cex.lab = 1.2,
  cex.axis = 1.2,
  col = "black",
  bg = module
)


############################################################
# 27. Save GS and MM results
############################################################

setwd("")

write.csv(GS, file = "LUAD_GS.csv")
write.csv(MM, file = "LUAD_MM.csv")


############################################################
# 28. Identify candidate hub genes
############################################################

# Candidate genes with high gene significance and module membership
candidate_hub_genes <- moduleGenes[
  (GS > 0.5 & MM > 0.3)
]

print(candidate_hub_genes)


############################################################
# 29. Export selected module network to Cytoscape
############################################################

# Recalculate TOM for module export
TOM <- TOMsimilarityFromExpr(
  data.mat,
  power = power,
  networkType = type
)

dimnames(TOM) <- list(
  colnames(data.mat),
  colnames(data.mat)
)

# Extract TOM for the selected module
modTOM <- TOM[moduleGenes, moduleGenes]

# Output file prefix
file <- "G:/HOSPITAL/TCGAData/WGCNA/LUAD.blue.module"

# Export module network
cyt <- exportNetworkToCytoscape(
  modTOM,
  edgeFile = paste0(file, ".edges.txt"),
  nodeFile = paste0(file, ".nodes.txt"),
  weighted = TRUE,
  threshold = 0,
  nodeNames = moduleGenes,
  nodeAttr = module
)


############################################################
# 30. Alternative TOM heatmap color scheme
############################################################

mycolor <- gplots::colorpanel(
  250,
  "red",
  "orange",
  "lemonchiffon"
)

TOMplot(
  plotTOM,
  net$dendrograms[[1]],
  moduleColors[net$blockGenes[[1]]],
  main = "Network heatmap plot, all genes",
  col = mycolor
)


############################################################
# 31. Export whole WGCNA network to Cytoscape
############################################################

file <- "G:/HOSPITAL/TCGAData/WGCNA/LUAD.net"

genes <- names(net$colors[net$blockGenes[[1]]])

dimnames(TOM) <- list(genes, genes)

cyt <- exportNetworkToCytoscape(
  TOM,
  edgeFile = paste0(file, ".edges.txt"),
  nodeFile = paste0(file, ".nodes.txt"),
  weighted = TRUE,
  threshold = 0,
  nodeNames = genes,
  nodeAttr = moduleColors[net$blockGenes[[1]]]
)