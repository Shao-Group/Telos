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

def stratified_split(df, label_col='label', test_size=0.2, seed=42):
    # X_train, X_val, y_train, y_val =  train_test_split(df, df[label_col], test_size=test_size, stratify=df[label_col], random_state=seed)
    # split based on chromosome 
    chromosomes = df['chrom'].unique()
    # train_chromosomes, val_chromosomes = train_test_split(chromosomes, test_size=test_size, random_state=seed)
    train_chromosomes = ["chr" + str(i) for i in range(1,15)]
    val_chromosomes = [x for x in chromosomes if x not in train_chromosomes]
    print(f"Train chromosomes: {train_chromosomes}")
    print(f"Validation chromosomes: {val_chromosomes}")
    train_mask = df['chrom'].isin(train_chromosomes)
    val_mask = df['chrom'].isin(val_chromosomes)
    X_train = df[train_mask].drop(columns=[label_col])
    X_val = df[val_mask].drop(columns=[label_col])
    y_train = df[train_mask][label_col]
    y_val = df[val_mask][label_col]

    print(f"Train size: {X_train.shape}, Validation size: {X_val.shape}")
    print(f"Train label distribution: {y_train.value_counts(normalize=True)}")
    print(f"Validation label distribution: {y_val.value_counts(normalize=True)}")
    return X_train, X_val, y_train, y_val
    

def evaluate_model(y_true, y_pred, y_prob, plot_path=None):
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
    prefix = os.path.splitext(plot_path)[0]

    pr_data.to_csv(f"{prefix}_pr_data.csv", index=False)
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
            eval_metric="logloss"
        )
    elif model_type == "lightgbm":
        import lightgbm as lgb
        return lgb.LGBMClassifier(
            n_estimators=config["n_estimators"],
            max_depth=config["max_depth"],
            learning_rate=config["learning_rate"],
            subsample=config["subsample"],
            colsample_bytree=config["colsample_bytree"]
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
