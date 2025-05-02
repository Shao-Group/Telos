import argparse
import yaml
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from joblib import dump

from utils.ml_utils import stratified_split, evaluate_model, load_model

def train_and_evaluate(df, model_type, config, site_type):
    drop = ["chrom", "position", "strand", "label"]
    feature_cols = [c for c in df.columns if c not in drop]
    # print(f"Feature columns: {feature_cols}")
    # print df column types
    # print(df.dtypes)
    X = df[feature_cols].select_dtypes(include=[int])
    y = df["label"]
    X_train, X_val, y_train, y_val = stratified_split(df)
    numeric_cols = df.select_dtypes(include=["number"]).columns.tolist()
    numeric_cols = [col for col in numeric_cols if col not in drop]
    X_train = X_train[numeric_cols]
    X_val = X_val[numeric_cols]

    model = load_model(model_type, config)
    model.fit(X_train, y_train)

    y_pred = model.predict(X_val)
    y_prob = model.predict_proba(X_val)[:, 1]

    model_dir = f"models/{site_type}"
    report_dir = f"reports/{site_type}"
    os.makedirs(model_dir, exist_ok=True)
    os.makedirs(report_dir, exist_ok=True)

    if model_type == "xgboost":
        model.save_model(f"{model_dir}/xgb_model.json")
    elif model_type == "lightgbm":
        model.booster_.save_model(f"{model_dir}/lgbm_model.txt")
    elif model_type == "randomforest":
        dump(model, os.path.join(model_dir, "rf_model.joblib"))

    metrics = evaluate_model(y_val, y_pred, y_prob, plot_path=f"{report_dir}/{model_type}_pr_curve.png")

    if model_type == "randomforest":
        importances = model.feature_importances_
        sorted_idx = np.argsort(importances)[::-1]
        sorted_feats = np.array(numeric_cols)[sorted_idx]
        sorted_vals = importances[sorted_idx]

        plt.figure(figsize=(8, 6))
        plt.barh(sorted_feats[:20][::-1], sorted_vals[:20][::-1])  # top 20
        plt.xlabel("Importance")
        plt.title(f"Top Feature Importances ({site_type.upper()} - RF)")
        plt.tight_layout()
        
        os.makedirs(report_dir, exist_ok=True)
        plt.savefig(os.path.join(report_dir, "feature_importance.png"))
        plt.close()

    with open(f"{report_dir}/{model_type}_metrics_summary.txt", "w") as f:    
        for k, v in metrics.items():
            if k == "confusion_matrix":
                f.write(f"{k}:\n")
                f.write(f"  TP: {v[1][1]}\n")
                f.write(f"  FP: {v[0][1]}\n")
                f.write(f"  FN: {v[1][0]}\n")
                f.write(f"  TN: {v[0][0]}\n")
            else:
                f.write(f"{k}: {v:.4f}\n")


    return metrics

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--config", required=True)
    parser.add_argument("--site_type", required=True, choices=["tss", "tes"])
    parser.add_argument("--model_type", required=True, choices=["xgboost", "lightgbm", "randomforest"])
    args = parser.parse_args()

    df = pd.read_csv(args.input)
    with open(args.config) as f:
        config = yaml.safe_load(f)

    metrics = train_and_evaluate(df, args.model_type, config, args.site_type)

    print(f"✅ {args.site_type.upper()} [{args.model_type}] - F1: {metrics['f1']:.4f}, AUPR: {metrics['aupr']:.4f}")

if __name__ == "__main__":
    main()
