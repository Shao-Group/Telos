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
    X_train, X_val, y_train, y_val =  train_test_split(df, df[label_col], test_size=test_size, stratify=df[label_col], random_state=seed)

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
