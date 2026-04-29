LUAD Biomarker Discovery Project

Overview

This study integrates weighted gene co-expression network analysis (WGCNA), survival analysis, and machine learning to identify prognostic biomarkers in lung adenocarcinoma (LUAD).

Workflow
1.Machine learning modeling (TabNet)
2.WGCNA module detection
3.Candidate gene selection (Survival analysis)
4.Venn
5.External validation using GEO datasets (GSE31210, GSE13213)
6.Functional analysis and protein interaction (STRING database)

Data Sources
TCGA-LUAD (Genomic Data Commons)
GSE31210 (GEO)
GSE13213 (GEO)


Repository Structure
/data/ – raw and processed datasets
/scripts/ – analysis scripts
/results/ – figures and tables


How to Run
Run scripts in the following order:
1.01_machine_learning_tabnet.py
2.02_WGCNA_analysis.R
3.03_gene_selection.R
4.04_Venn.R
5.05_external_validation.R


External Tools
STRING database for protein-protein interaction analysis: https://string-db.org/

Output
Results stored in /results/
Figures stored in /results/figures/

Reproducibility
All code and processed data required to reproduce the results are provided in this repository.
Users should install the required R and Python packages before running the scripts.

Recommended environment:
R ≥ 4.2
Python ≥ 3.9

License
This project is released under the CC-BY 4.0 license.


Large datasets are available at:
https://drive.google.com/drive/folders/12uFXK0SLzo4nl-CSISMZKrwTPNSKbZiO?usp=sharing
