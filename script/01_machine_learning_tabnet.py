# -*- coding: utf-8 -*-
"""
LUAD Machine Learning + Survival Analysis Pipeline

Author: Cheng-Hong Yang

Description:
This script integrates:
1. TabNet classification (baseline + weighted features)
2. Cox proportional hazards model
3. ROC curve comparison across models
4. Performance evaluation (Accuracy, Sensitivity, Specificity, F1, AUC)

Workflow:
- Load LUAD dataset
- Apply feature weighting (Additive / Multiplicative)
- Train TabNet models
- Fit CoxPH survival model
- Compare performance using ROC curves
"""

############################################################
# 1. Import required libraries
############################################################

import pandas as pd
import numpy as np

from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelEncoder
from sklearn.metrics import (
    accuracy_score, recall_score, f1_score,
    confusion_matrix, roc_curve, auc
)

from pytorch_tabnet.tab_model import TabNetClassifier
from lifelines import CoxPHFitter

import matplotlib.pyplot as plt
import matplotlib
matplotlib.use("TkAgg")

############################################################
# 2. File paths
############################################################

MODEL_PATH_LUAD = r"D:/HOSPITAL/.../tabnet_luad_model.zip"
DATA_PATH_LUAD  = r"D:/HOSPITAL/.../Combied_ML_LUAD.csv"
RISK_PATH_LUAD  = r"D:/HOSPITAL/.../Tab_risk_LUAD.csv"

############################################################
# 3. Load data
############################################################

df_LUAD   = pd.read_csv(DATA_PATH_LUAD)
risk_LUAD = pd.read_csv(RISK_PATH_LUAD)

############################################################
# 4. Feature matrix and label encoding
############################################################

# Feature matrix: remove first 2 columns and last 2 columns
X_base_LUAD = df_LUAD.iloc[:, 2:-2].values

# Target variable (binary survival status)
y_raw_LUAD  = df_LUAD["os"].values
y_enc_LUAD  = LabelEncoder().fit_transform(y_raw_LUAD)

############################################################
# 5. Risk weighting (Additive / Multiplicative)
############################################################

weights_add  = risk_LUAD["ADD_Quartile"].values
weights_mult = risk_LUAD["Mult_Quartile"].values

# Apply weighting (broadcasting)
X_add_LUAD  = X_base_LUAD * weights_add
X_mult_LUAD = X_base_LUAD * weights_mult

############################################################
# 6. Train-test split
############################################################

split_kws = dict(test_size=0.2, random_state=42, stratify=y_enc_LUAD)

Xtr, Xte, ytr, yte = train_test_split(X_base_LUAD, y_enc_LUAD, **split_kws)
Xtr_add, Xte_add, _, _ = train_test_split(X_add_LUAD, y_enc_LUAD, **split_kws)
Xtr_mult, Xte_mult, _, _ = train_test_split(X_mult_LUAD, y_enc_LUAD, **split_kws)

############################################################
# 7. Load pre-trained TabNet model
############################################################

tabnet_model = TabNetClassifier()
tabnet_model.load_model(MODEL_PATH_LUAD)

############################################################
# 8. Cox Proportional Hazards Model
############################################################

# Prepare dataset for Cox model
X_cox = df_LUAD.iloc[:, 2:].dropna(axis=1)

cox_model = CoxPHFitter()

train_data, test_data = train_test_split(X_cox, test_size=0.2, random_state=42)

# Fit Cox model
cox_model.fit(train_data, duration_col="os_time", event_col="os")

# Print summary
cox_model.print_summary()

############################################################
# 9. Cox model prediction and evaluation
############################################################

risk_scores_test = cox_model.predict_partial_hazard(
    test_data.drop(columns=["os_time", "os"])
)

# Binary classification using threshold = 1
pred_cox = (risk_scores_test > 1).astype(int)

y_true = test_data["os"].values

# Metrics
acc_cox  = accuracy_score(y_true, pred_cox)
sens_cox = recall_score(y_true, pred_cox)

tn, fp, fn, tp = confusion_matrix(y_true, pred_cox).ravel()
spec_cox = tn / (tn + fp)
f1_cox   = f1_score(y_true, pred_cox)

# ROC
cox_fpr, cox_tpr, _ = roc_curve(y_true, risk_scores_test)
roc_auc_cox = auc(cox_fpr, cox_tpr)

############################################################
# 10. Train TabNet models (Add / Mult weighting)
############################################################

# ADD weighting model
model_add = TabNetClassifier()
model_add.fit(
    Xtr_add, ytr,
    eval_set=[(Xte_add, yte)],
    max_epochs=10, patience=10
)

# MULT weighting model
model_mult = TabNetClassifier()
model_mult.fit(
    Xtr_mult, ytr,
    eval_set=[(Xte_mult, yte)],
    max_epochs=50, patience=10
)

############################################################
# 11. Prediction probabilities
############################################################

proba_tabnet = tabnet_model.predict_proba(Xte)[:, 1]
proba_add    = model_add.predict_proba(Xte_add)[:, 1]
proba_mult   = model_mult.predict_proba(Xte_mult)[:, 1]

############################################################
# 12. Model comparison (ROC + metrics)
############################################################

models = {
    "TabNet"      : proba_tabnet,
    "TabNet_Add"  : proba_add,
    "TabNet_Mult" : proba_mult
}

plt.figure(figsize=(5,5))

records = []

# Add CoxPH result
plt.plot(cox_fpr, cox_tpr, lw=1.5,
         label=f"CoxPH (AUC={roc_auc_cox:.2f})")

records.append({
    "Model": "CoxPH",
    "Accuracy": round(acc_cox, 3),
    "Sensitivity": round(sens_cox, 3),
    "Specificity": round(spec_cox, 3),
    "F1-score": round(f1_cox, 3),
    "AUC": round(roc_auc_cox, 3)
})

# Loop through ML models
for name, score in models.items():

    fpr, tpr, _ = roc_curve(y_true, score)
    auc_val = auc(fpr, tpr)

    plt.plot(fpr, tpr, label=f"{name} (AUC={auc_val:.2f})")

    preds = (score > 0.5).astype(int)

    acc  = accuracy_score(y_true, preds)
    sens = recall_score(y_true, preds)

    tn, fp, fn, tp = confusion_matrix(y_true, preds).ravel()
    spec = tn / (tn + fp)

    f1 = f1_score(y_true, preds)

    records.append({
        "Model": name,
        "Accuracy": round(acc, 3),
        "Sensitivity": round(sens, 3),
        "Specificity": round(spec, 3),
        "F1-score": round(f1, 3),
        "AUC": round(auc_val, 3)
    })

############################################################
# 13. Plot ROC curve
############################################################

plt.plot([0,1], [0,1], "--", lw=0.8)

plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
plt.title("ROC – LUAD: CoxPH vs TabNet")
plt.legend(loc="lower right")
plt.tight_layout()
plt.show()

############################################################
# 14. Output results
############################################################

results_df = pd.DataFrame(records).sort_values("AUC", ascending=False)

print(results_df)

# Optional: save results
# results_df.to_csv("combined_metrics.csv", index=False)