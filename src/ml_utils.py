import joblib
import json
import os
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import (
    precision_recall_curve, roc_auc_score, average_precision_score,
    accuracy_score, precision_score, recall_score, f1_score,
    confusion_matrix
)
import matplotlib.pyplot as plt
import numpy as np
import re


def normalize_chrom(c):
    """
    Turn any of 'chr1','1','CHR1' → '1', but leave 'GL000194.1','X','MT' alone.
    """
    s = str(c)
    # remove leading 'chr' (case-insensitive)
    s = re.sub(r'(?i)^chr', '', s)
    return s

def chrom_to_int(c):
    """
    Return the chromosome as an int if c is exactly '1'…'22'; else None.
    """
    s = str(c).lstrip('chr').upper()
    # fullmatch: entire string must be 1–2 digits
    if re.fullmatch(r'[1-9]|1[0-9]|2[0-2]', s):
        val = int(s)
        # only accept 1–15 for your use-case
        return val
    return None


def load_tmap_labels(tmap_path):
    df = pd.read_csv(tmap_path, sep='\t', comment='#', header=0)
    df = df[['qry_id', 'class_code']].rename(columns={'qry_id': 'transcript_id'})
    df['label'] = (df['class_code'] == '=').astype(int)
    return df[['transcript_id', 'label']]


def stratified_split(df, validation_chrom_file, label_col='label', ):
    # X_train, X_val, y_train, y_val =  train_test_split(df, df[label_col], test_size=test_size, stratify=df[label_col], random_state=seed)
    # split based on chromosome 
    df['chrom_number'] = df['chrom'].apply(chrom_to_int)
    
    train_mask = df['chrom_number'].between(1, 5)
    val_mask = ~train_mask

    df.drop(columns=['chrom_number'], inplace=True)

    X_train = df[train_mask].drop(columns=[label_col])
    X_val = df[val_mask].drop(columns=[label_col])
    y_train = df[train_mask][label_col]
    y_val = df[val_mask][label_col]

    with open(validation_chrom_file, 'w') as f:
        for chrom in X_val['chrom'].unique():
            f.write(f"{chrom}\n")

    print(f"Train size: {X_train.shape}, Validation size: {X_val.shape}")
    print(f"Train label distribution: {y_train.value_counts(normalize=True)}")
    print(f"Validation label distribution: {y_val.value_counts(normalize=True)}")
    return X_train, X_val, y_train, y_val
    

def evaluate_model(y_true, y_pred, y_prob, prdata_path, plot_path=None):
    precision, recall, _ = precision_recall_curve(y_true, y_prob)
    f1 = f1_score(y_true, y_pred)
    aupr = average_precision_score(y_true, y_prob)
    auc = roc_auc_score(y_true, y_prob)
    acc = accuracy_score(y_true, y_pred)
    prec = precision_score(y_true, y_pred, zero_division=0)
    rec = recall_score(y_true, y_pred)
    cm = confusion_matrix(y_true, y_pred)

    # PR curve dataframe
    pr_data = pd.DataFrame({"precision": precision, "recall": recall})

    pr_data.to_csv(prdata_path, index=False)
    if plot_path:
        plt.figure()
        plt.plot(recall, precision, label=f"PR AUC={aupr:.3f}")
        plt.xlabel("Recall")
        plt.ylabel("Precision")
        plt.title("Precision-Recall Curve")
        plt.grid()
        plt.savefig(plot_path)
        plt.close()

    return {
        "accuracy": acc,
        "precision": prec,
        "recall": rec,
        "f1": f1,
        "aupr": aupr,
        "auc": auc,
        "confusion_matrix": cm.tolist()
    }

def load_model(model_type, config):
    if model_type == "xgboost":
        import xgboost as xgb
        return xgb.XGBClassifier(
            n_estimators=config["n_estimators"],
            max_depth=config["max_depth"],
            learning_rate=config["learning_rate"],
            subsample=config["subsample"],
            colsample_bytree=config["colsample_bytree"],
            # use_label_encoder=False,
            objective="binary:logistic",
            eval_metric="aucpr"
        )

    elif model_type == "randomforest":
        return RandomForestClassifier(
            n_estimators=config["n_estimators"],
            max_depth=config["max_depth"],
            n_jobs=-1,
            random_state=42
        )
    else:
        raise ValueError(f"Unsupported model type: {model_type}")


def load_saved_model(model_type, model_path):
    if model_type == "xgboost":
        return json.load(open(model_path))  
    elif model_type == "randomforest":
        return joblib.load(model_path)
    else:
        raise ValueError(f"Unsupported model type: {model_type}")