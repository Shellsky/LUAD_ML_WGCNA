############################################################
# Gene Overlap Analysis from Excel Lists
#
# Description:
# This script reads multiple gene lists from an Excel file
# and performs:
# 1. Intersection (genes shared across all lists)
# 2. Frequency analysis (genes appearing ≥ N times)
# 3. Optional Venn diagram visualization
############################################################


############################################################
# 1. Load required packages
############################################################

library(readxl)       # For reading Excel files
library(dplyr)        # Data manipulation
library(VennDiagram)  # For Venn diagram visualization
library(grid)         # For plotting Venn diagram


############################################################
# 2. Load gene list from Excel file
############################################################

# Replace with your actual file path
gene_df <- read_excel("Gene.xlsx")


############################################################
# 3. Convert each column into a gene list
############################################################

# Convert each column to character vector and remove NA values
gene_lists <- lapply(gene_df, function(col) {
  na.omit(as.character(col))
})


############################################################
# 4. Identify common genes across all lists
############################################################

# Find genes present in all columns (strict intersection)
common_all <- Reduce(intersect, gene_lists)

cat("Genes present in ALL gene lists:\n")
print(common_all)


############################################################
# 5. Identify genes appearing at least N times
############################################################

# Combine all gene lists into one vector
gene_union <- unlist(gene_lists)

# Count frequency of each gene
gene_freq_table <- table(gene_union)

# Extract genes appearing at least 2 times
dup_2plus <- names(gene_freq_table[gene_freq_table >= 2])

cat("\nGenes appearing in at least TWO lists:\n")
print(dup_2plus)


############################################################
# 6. (Optional) Visualize overlap using Venn Diagram
############################################################

venn.plot <- venn.diagram(
  x = gene_lists,
  category.names = names(gene_df),
  filename = NULL,          # Set to "venn_gene.png" to save file
  output = TRUE,
  imagetype = "png",
  euler.d = TRUE
)

grid.draw(venn.plot)