############################################################
# LUAD Survival Analysis and Gene Signature Evaluation
# Description:
# This script performs survival analysis using TCGA-LUAD
# expression data and clinical survival information.
# The workflow includes:
# 1. Loading expression and clinical data
# 2. Selecting tumor samples
# 3. Matching expression profiles with clinical records
# 4. Constructing overall survival variables
# 5. Calculating gene-based risk scores
# 6. Performing Kaplan-Meier survival analysis
# 7. Performing Cox proportional hazards regression
############################################################


############################################################
# 1. Load required packages
############################################################

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
library(readr)
library(glmnet)
library(survival)
library(survminer)
library(readxl)


############################################################
# 2. Clear environment and load LUAD expression data
############################################################

rm(list = ls())

# Load LUAD TPM expression data
load(file = "LUAD_TPM.Rdata")


############################################################
# 3. Load gene list and identify common genes
############################################################

# Read gene list from Excel file
gene_df <- read_excel("Gene.xlsx")

# Convert each column into a gene list and remove missing values
gene_lists <- lapply(gene_df, function(col) na.omit(as.character(col)))

# Identify genes shared across all gene lists
common_all <- Reduce(intersect, gene_lists)


############################################################
# 4. Prepare clinical survival information
############################################################

# Extract patient barcode and vital status
lc <- lc_data %>%
  dplyr::select(bcr_patient_barcode, vital_status)

# Convert vital status into binary format
# Dead = 1, Alive = 0
family <- lc %>%
  mutate(vital_status = ifelse(vital_status == "Dead", 1, 0))


############################################################
# 5. Select TCGA tumor samples
############################################################

# TCGA barcode positions 14-15 indicate sample type
# Values < 10 represent tumor samples
tumor <- colnames(exp)[as.integer(substr(colnames(exp), 14, 15)) < 10]

# Subset expression matrix to tumor samples only
tumor_sample <- exp[, tumor]

# Create a tumor sample annotation table
Tumor <- data.frame(ID = tumor, Value = 1)

# Replace expression matrix with tumor-only expression data
exp <- tumor_sample

# Shorten TCGA barcode to the first 12 characters for patient-level matching
colnames(exp) <- substr(colnames(exp), 1, 12)


############################################################
# 6. Match expression data with clinical data
############################################################

# Identify patients with both expression and clinical data
matching_barcodes <- lc[lc$bcr_patient_barcode %in% colnames(exp), ]

# Keep only expression samples with matched clinical data
exp <- exp[, colnames(exp) %in% matching_barcodes$bcr_patient_barcode]

# Store sample barcodes
barcode <- colnames(exp)


############################################################
# 7. Create expression-vital status mapping table
############################################################

exp_vital_status <- data.frame(
  barcode = colnames(exp),
  vital_status = sapply(barcode, function(barcode) {
    status <- lc$vital_status[lc$bcr_patient_barcode == substr(barcode, 1, 12)]
    
    if (length(status) == 0) {
      message("No vital_status found for barcode: ", barcode)
      return(NA)
    } else {
      return(status[1])
    }
  })
)

print(exp_vital_status)

# Convert survival status into binary session variable
test_datafrmae <- exp_vital_status %>%
  mutate(session = ifelse(vital_status == "Dead", 1, 0))


############################################################
# 8. Function: Survival analysis based on mutation-associated gene
############################################################

run_survival_analysis <- function(gene_name, data, filtered_data, merged_data) {
  
  # Extract mutation data for the selected gene
  filtered_gene_data <- data %>%
    select(Hugo_Symbol, Tumor_Sample_Barcode) %>%
    filter(Hugo_Symbol == gene_name)
  
  # Convert tumor sample barcode to patient-level barcode
  filtered_gene_data$Tumor_Sample_Barcode <- substr(
    filtered_gene_data$Tumor_Sample_Barcode, 1, 12
  )
  
  # Remove duplicated patient barcodes
  filtered_gene_data <- filtered_gene_data %>%
    distinct(Tumor_Sample_Barcode, .keep_all = TRUE)
  
  # Identify barcodes with both mutation and expression data
  matching_barcodes <- filtered_gene_data$Tumor_Sample_Barcode[
    filtered_gene_data$Tumor_Sample_Barcode %in% colnames(filtered_data)
  ]
  
  filtered_data_matching <- filtered_data %>%
    select(all_of(matching_barcodes))
  
  # Filter survival data for matched patients
  merged_data_os <- merged_data %>%
    filter(bcr_patient_barcode %in% matching_barcodes)
  
  print(gene_name)
  
  # Display vital status distribution
  ostable <- table(merged_data_os$vital_status)
  
  # Identify common patients between datasets
  common_barcodes <- intersect(
    merged_data_test$bcr_patient_barcode,
    merged_data_os$bcr_patient_barcode
  )
  
  # Prepare survival dataset
  merged_data_test_filtered <- merged_data_test %>%
    filter(bcr_patient_barcode %in% common_barcodes) %>%
    mutate(OS_TIME = OS_TIME / 30.4375)
  
  # Check whether survival data are available
  if (nrow(merged_data_test_filtered) == 0) {
    warning(paste("No matching data for gene:", gene_name))
    return(NULL)
  }
  
  # Check missing survival information
  if (any(is.na(merged_data_test_filtered$OS_TIME)) ||
      any(is.na(merged_data_test_filtered$vital_status))) {
    warning(paste("Missing OS_TIME or vital_status for gene:", gene_name))
    return(NULL)
  }
  
  # Fit Kaplan-Meier survival curve
  surv_fit <- survfit(
    Surv(OS_TIME, vital_status) ~ risk_group,
    data = merged_data_test_filtered
  )
  
  # Save survival model
  surv_fit_list[[gene_name]] <- surv_fit
  
  # Perform log-rank test
  surv_diff <- survdiff(
    Surv(OS_TIME, vital_status) ~ risk_group,
    data = merged_data_test_filtered
  )
  
  # Calculate p-value
  p_value <- 1 - pchisq(surv_diff$chisq, length(surv_diff$n) - 1)
  
  # Save p-value
  p_value_list[[gene_name]] <- p_value
  
  # Calculate median overall survival time by risk group
  median_os_high <- median(
    merged_data_test_filtered$OS_TIME[
      merged_data_test_filtered$risk_group == 1
    ],
    na.rm = TRUE
  )
  
  median_os_low <- median(
    merged_data_test_filtered$OS_TIME[
      merged_data_test_filtered$risk_group == 0
    ],
    na.rm = TRUE
  )
  
  print(paste("Gene:", gene_name))
  print(paste("High Risk Median OS Time:", round(median_os_high, 1), "Month"))
  print(paste("Low Risk Median OS Time:", round(median_os_low, 1), "Month"))
  
  # Plot Kaplan-Meier survival curve
  print(
    ggsurvplot(
      surv_fit,
      data = merged_data_test_filtered,
      pval = TRUE,
      risk.table = TRUE,
      ggtheme = theme_minimal(),
      title = gene_name
    )
  )
  
  return(list(p_value = p_value, ostable = ostable))
}


############################################################
# 9. Function: Survival analysis based on single-gene expression
############################################################

run_survival_analysis_people <- function(
    gene_name,
    merged_data,
    merged_data_test,
    median_sum
) {
  
  # Extract expression value of the selected gene
  merged_data_test_gene <- as.data.frame(merged_data_test[[gene_name]])
  colnames(merged_data_test_gene) <- gene_name
  
  print(gene_name)
  
  # Extract survival outcome variables
  merged_y <- merged_data_test[, (ncol(merged_data_test) - 1):ncol(merged_data_test)]
  
  # Create risk group according to median cutoff
  risk_scores_df <- as.data.frame(merged_data_test_gene)
  median_value <- median_sum
  
  risk_scores_df$risk_group <- ifelse(
    risk_scores_df[[gene_name]] > median_value,
    "High",
    "Low"
  )
  
  # Combine gene expression, survival information, and risk group
  merged_data_filtered <- cbind(
    merged_data_test_gene,
    merged_y,
    risk_scores_df
  )
  
  # Fit Kaplan-Meier survival curve
  surv_fit <- survfit(
    Surv(OS_TIME, vital_status) ~ risk_group,
    data = merged_data_filtered
  )
  
  # Save survival model
  surv_fit_list[[gene_name]] <- surv_fit
  
  # Perform log-rank test
  surv_diff <- survdiff(
    Surv(OS_TIME, vital_status) ~ risk_group,
    data = merged_data_filtered
  )
  
  # Calculate p-value
  p_value <- 1 - pchisq(surv_diff$chisq, length(surv_diff$n) - 1)
  
  # Save p-value
  p_value_list[[gene_name]] <- p_value
  
  # Plot Kaplan-Meier survival curve
  print(
    ggsurvplot(
      surv_fit,
      data = merged_data_filtered,
      pval = TRUE,
      risk.table = TRUE,
      ggtheme = theme_minimal(),
      title = gene_name
    )
  )
  
  return(p_value)
}


############################################################
# 10. Construct overall survival dataset
############################################################

cox_lc <- lc_data %>%
  dplyr::select(
    bcr_patient_barcode,
    vital_status,
    days_to_last_follow_up,
    days_to_death
  )

# Convert vital status into binary outcome
cox_lc <- cox_lc %>%
  mutate(vital_status = ifelse(vital_status == "Dead", 1, 0))

# Define overall survival time
Final <- cox_lc %>%
  mutate(OS_TIME = case_when(
    vital_status == 0 ~ days_to_last_follow_up,
    vital_status == 1 ~ days_to_death,
    TRUE ~ NA_real_
  )) %>%
  select(bcr_patient_barcode, vital_status, OS_TIME)


############################################################
# 11. Initialize result containers
############################################################

surv_fit_list <- list()
p_value_list <- list()

results_df <- data.frame(
  n = integer(),
  seednumber = integer(),
  gene = character(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)


############################################################
# 12. Define candidate genes
############################################################

common_genes <- c(
  "BPIFB1", "CTNNA1", "HNRNPF", "FN1", "CD14", "IGHA1", "CCND1",
  "CMPK1", "IGKV1-16", "COL3A1", "EIF4G2", "GPRC5A", "FOS", "CES1",
  "COL1A1", "IGHM", "SEPTIN2", "EIF4G1", "PDIA4", "TIMP3", "HSP90B1",
  "COX5B", "SF3B4", "APOL1", "MAGED2", "FLNA", "AEBP1", "IGHG1",
  "XRCC6", "IGKV1-17", "IGHV1-2", "EGR1", "MORF4L2", "IGHG4", "TPI1",
  "CDH1", "MSN", "COL1A2", "CANX", "ATP1B1", "EPAS1", "SLC34A2",
  "HBB", "NCL", "ANXA5", "RPLP2", "A2M", "P4HB", "SFTPA2", "SEC61A1",
  "CYC1", "MYH9", "RPL6", "PLXNB2", "GANAB", "EEF1A1", "C4BPA", "COL6A1",
  "ZFP36L1"
)


############################################################
# 13. Filter expression matrix by selected genes
############################################################

filtered_data <- exp[rownames(exp) %in% common_genes, ]

# Transpose expression matrix: rows = patients, columns = genes
x <- t(filtered_data)

x <- as.data.frame(x)

# Add patient barcode as a column
x <- x %>%
  rownames_to_column(var = "bcr_patient_barcode")

# Merge gene expression data with survival data
merged_data <- x %>%
  left_join(Final, by = "bcr_patient_barcode") %>%
  drop_na()

# Extract gene expression matrix
x <- merged_data[, -c(1, (ncol(merged_data) - 1):ncol(merged_data))]

# Extract survival status
y <- merged_data$vital_status
y <- as.numeric(y)

# Convert expression data to matrix
x <- as.matrix(x)


############################################################
# 14. Calculate unweighted gene signature risk score
############################################################

# Select genes for risk score calculation
selected_genes <- common_genes

# Assign equal coefficient to all selected genes
weights <- data.frame(
  Gene = selected_genes,
  Coefficient = rep(1, length(selected_genes))
)

rownames(weights) <- selected_genes

# Calculate gene-wise weighted expression score
risk_scores <- sweep(
  x[, rownames(weights)],
  2,
  weights$Coefficient,
  `*`
)

risk_scores_df <- as.data.frame(risk_scores)

# Calculate total risk score by summing all gene scores
risk_scores_sum <- rowSums(risk_scores_df)

risk_scores_df$sum <- risk_scores_sum

# Use median risk score as cutoff
median_sum <- median(risk_scores_df$sum, na.rm = TRUE)

# Define high-risk and low-risk groups
risk_scores_df$risk_group <- ifelse(
  risk_scores_df$sum > median_sum,
  "High",
  "Low"
)


############################################################
# 15. Merge risk score with survival data
############################################################

merged_data_lasso <- merged_data %>%
  select(bcr_patient_barcode, vital_status, OS_TIME)

merged_data_test <- bind_cols(risk_scores_df, merged_data_lasso)

# Restrict analysis to patients with OS time less than 5 years
merged_data_test_5year <- merged_data_test %>%
  filter(OS_TIME < 1825)


############################################################
# 16. Kaplan-Meier survival analysis for total risk score
############################################################

surv_fit <- survfit(
  Surv(OS_TIME, vital_status) ~ risk_group,
  data = merged_data_test_5year
)

ggsurvplot(
  surv_fit,
  data = merged_data_test_5year,
  pval = TRUE,
  risk.table = TRUE,
  ggtheme = theme_minimal(),
  title = "No Threshold"
)


############################################################
# 17. Function: Cox proportional hazards analysis for each gene
############################################################

COXPH_Peolpe <- function(gene_name, merged_data) {
  
  # Extract selected gene expression
  merged_data_gene <- as.data.frame(merged_data[[gene_name]])
  colnames(merged_data_gene) <- gene_name
  
  # Extract survival outcome variables
  merged_y <- merged_data[, (ncol(merged_data) - 1):ncol(merged_data)]
  
  # Create median-based risk group
  risk_scores_df <- as.data.frame(merged_data_gene)
  median_value <- median(merged_data_gene[[gene_name]])
  
  risk_scores_df$risk_group <- ifelse(
    risk_scores_df[[gene_name]] > median_value,
    "High",
    "Low"
  )
  
  # Combine gene expression, survival outcome, and risk group
  merged_data_filtered <- cbind(
    merged_data_gene,
    merged_y,
    risk_scores_df
  )
  
  merged_data_filtered$risk_group <- as.factor(merged_data_filtered$risk_group)
  
  # Fit Cox proportional hazards model
  cox_model <- coxph(
    Surv(OS_TIME, vital_status) ~ risk_group,
    data = merged_data_filtered
  )
  
  # Extract HR and p-value
  cox_summary <- summary(cox_model)
  hr_values <- cox_summary$coefficients[, "exp(coef)"]
  p_values <- cox_summary$coefficients[, "Pr(>|z|)"]
  
  print(
    data.frame(
      gene = gene_name,
      Variable = rownames(cox_summary$coefficients),
      HR = hr_values,
      P_value = p_values
    )
  )
  
  return(list(hr_values, p_values))
}


############################################################
# 18. Run CoxPH analysis for all candidate genes
############################################################

results_coxph <- data.frame(
  gene = character(),
  HR = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

genes <- common_genes

for (gene in genes) {
  tryCatch({
    
    result_cox <- COXPH_Peolpe(gene, merged_data_test_5year)
    
    result_cox_df <- data.frame(
      gene = gene,
      HR = result_cox[1],
      p_value = result_cox[2],
      stringsAsFactors = FALSE
    )
    
    colnames(result_cox_df) <- c("gene", "HR", "p_value")
    
    results_coxph <- rbind(results_coxph, result_cox_df)
    
  }, error = function(e) {
    message(sprintf(
      "Error in processing gene: %s, error: %s",
      gene,
      e$message
    ))
  })
}


############################################################
# 19. Multivariable Cox proportional hazards model
############################################################

cox_model <- coxph(
  Surv(OS_TIME, vital_status) ~ FLNC + KRT6A + TNNC2,
  data = merged_data
)

summary(cox_model)

# Extract model coefficients
coefficients <- coef(cox_model)


############################################################
# 20. Check zero-value distribution across variables
############################################################

# Count zero values in each column
zero_counts <- colSums(merged_data_test == 0)

# Total number of rows
total_counts <- nrow(merged_data_test)

# Calculate zero ratio
zero_ratios <- zero_counts / total_counts

# Summarize zero-value distribution
zero_summary <- data.frame(
  Column = names(zero_counts),
  Zero_Count = zero_counts,
  Zero_Ratio = zero_ratios
)

print(zero_summary)


############################################################
# 21. Run survival analysis for selected genes
############################################################

genes <- c("SF3B4", "SLC34A2")

for (gene in genes) {
  
  print(gene)
  
  tryCatch({
    
    result <- run_survival_analysis(
      gene,
      LUAD_EGFR,
      filtered_data,
      merged_data_test_5year
    )
    
    p <- result[[1]]
    ostable <- result[[2]]
    
    print(ostable)
    
    results <- rbind(
      results_df,
      data.frame(
        n = n,
        seednumber = seednumbers,
        gene = gene,
        p_value = p
      )
    )
    
    p <- run_survival_analysis_people(
      gene,
      merged_data,
      merged_data_test,
      median_sum
    )
    
    print(p)
    
  }, error = function(e) {
    message(sprintf(
      "Error in processing gene: %s, error: %s",
      gene,
      e$message
    ))
  })
}


############################################################
# 22. Single-gene risk score analysis: SF3B4
############################################################

# Extract coefficient of SF3B4
sf3b4_coef <- weights["SF3B4", "Coefficient"]

# Calculate SF3B4-based risk score
risk_scores <- x[, "SF3B4"] * sf3b4_coef

# Convert risk score into data frame
risk_scores_df <- data.frame(SF3B4_score = risk_scores)

# Since only one gene is used, the total score equals SF3B4 score
risk_scores_df$sum <- risk_scores_df$SF3B4_score

# Define median cutoff
median_sum <- median(risk_scores_df$sum, na.rm = TRUE)

# Define high- and low-risk groups
risk_scores_df$risk_group <- ifelse(
  risk_scores_df$sum > median_sum,
  "High",
  "Low"
)

# Merge SF3B4 risk score with survival data
merged_data_lasso <- merged_data %>%
  select(bcr_patient_barcode, vital_status, OS_TIME)

merged_data_test <- bind_cols(risk_scores_df, merged_data_lasso)

# Restrict analysis to patients with follow-up less than 5 years
merged_data_test_5year <- merged_data_test %>%
  filter(OS_TIME < 1825)


############################################################
# 23. Kaplan-Meier survival analysis for SF3B4
############################################################

surv_data <- merged_data_test_5year

# Create survival object
surv_object <- Surv(
  time = surv_data$OS_TIME,
  event = surv_data$vital_status
)

# Fit Kaplan-Meier model
fit <- survfit(
  surv_object ~ risk_group,
  data = surv_data
)

# Plot Kaplan-Meier survival curve
ggsurvplot(
  fit,
  data = surv_data,
  pval = TRUE,
  conf.int = TRUE,
  risk.table = TRUE,
  legend.title = "Risk Group",
  legend.labs = c("High", "Low"),
  xlab = "Days",
  ylab = "Survival Probability",
  title = "SF3B4",
  palette = c("#E64B35", "#4DBBD5")
)