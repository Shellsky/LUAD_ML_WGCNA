############################################################
# External Validation Survival Analysis: GSE3141
# Gene of Interest: SF3B4
#
# Description:
# This script performs Kaplan-Meier survival analysis using
# the GSE3141 validation cohort. Patients are divided into
# high- and low-risk groups according to the median expression
# level of SF3B4.
############################################################


############################################################
# 1. Load required packages
############################################################

library(DESeq2)
library(pheatmap)
library(ggplot2)
library(tidyverse)
library(SummarizedExperiment)
library(dplyr)
library(tibble)
library(stringr)
library(tidyr)
library(limma)
library(ggpubr)
library(RColorBrewer)
library(pROC)
library(survival)
library(survminer)


############################################################
# 2. Load clinical and expression data
############################################################

# Load clinical survival data
clinical_raw <- read_csv(
  "D:/HOSPITAL/張醫生資料/GSE資料庫/GSE3141/GSE3141_Clinical_ch.csv"
)

# Load gene expression matrix
expression_raw <- read_csv(
  "D:/HOSPITAL/張醫生資料/GSE資料庫/GSE3141/GSE3141.csv"
)


############################################################
# 3. Prepare clinical data
############################################################

# Select sample ID, overall survival time, and survival status
clinical <- clinical_raw[c("geo_accession", "OS_TIME", "OS")]

# Rename sample ID column
colnames(clinical)[1] <- "ID"

# Convert survival time to numeric format
clinical$OS_TIME <- as.numeric(clinical$OS_TIME)


############################################################
# 4. Prepare expression matrix
############################################################

# Rename the first column as gene identifier
colnames(expression_raw)[1] <- "columns"

# Make duplicated gene names unique
expression_raw$columns <- make.unique(expression_raw$columns)

# Remove missing values
expression_raw <- na.omit(expression_raw)

# Set gene names as row names
expression_raw <- column_to_rownames(
  expression_raw,
  var = "columns"
)

# Keep only samples with available clinical data
data_ID <- clinical$ID
expression_matched <- expression_raw[, data_ID]

# Transpose expression matrix:
# rows = samples, columns = genes
expression_data <- t(expression_matched)


############################################################
# 5. Extract gene of interest
############################################################

# Define target gene
target_gene <- "SF3B4"

# Check whether target gene exists in the dataset
gene_id <- colnames(expression_data) %in% target_gene

# Extract target gene expression
filtered_data <- expression_data[, gene_id, drop = FALSE]


############################################################
# 6. Combine gene expression with clinical data
############################################################

# Convert expression matrix to data frame
risk_scores <- as.data.frame(filtered_data)

# Combine expression values and clinical survival data
validation_data <- bind_cols(risk_scores, clinical)

# Create risk score using SF3B4 expression
validation_data$sum <- validation_data[[target_gene]]


############################################################
# 7. Define high- and low-risk groups
############################################################

# Use median SF3B4 expression as cutoff
median_cutoff <- median(validation_data$sum, na.rm = TRUE)

# Assign patients to high- or low-risk groups
validation_data$risk_group <- ifelse(
  validation_data$sum > median_cutoff,
  "High",
  "Low"
)


############################################################
# 8. Restrict analysis to 5-year overall survival
############################################################

# Keep patients with survival time less than 5 years
validation_data <- validation_data %>%
  filter(OS_TIME < 1825)


############################################################
# 9. Calculate median OS time by risk group
############################################################

median_os_high <- median(
  validation_data$OS_TIME[validation_data$risk_group == "High"],
  na.rm = TRUE
)

median_os_low <- median(
  validation_data$OS_TIME[validation_data$risk_group == "Low"],
  na.rm = TRUE
)

print(paste(
  "High Risk Median OS Time:",
  round(median_os_high, 1),
  "Days"
))

print(paste(
  "Low Risk Median OS Time:",
  round(median_os_low, 1),
  "Days"
))


############################################################
# 10. Kaplan-Meier survival analysis
############################################################

surv_fit <- survfit(
  Surv(OS_TIME, OS) ~ risk_group,
  data = validation_data
)


############################################################
# 11. Plot Kaplan-Meier survival curve
############################################################

ggsurvplot(
  surv_fit,
  data = validation_data,
  pval = TRUE,
  risk.table = TRUE,
  ggtheme = theme_minimal(),
  title = "GSE3141"
)