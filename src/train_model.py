import argparse
import yaml
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from joblib import dump
from sklearn.preprocessing import StandardScaler
from config import config

from utils.ml_utils import stratified_split, evaluate_model, load_model

def train_and_evaluate(df, model_type, config, site_type, tool, out_dir, val_chrom_file):
    drop = ["chrom", "position", "strand", "label" , "soft_clip_entropy"]
    # drop += ['read_start_density', 'read_end_density', 'mean_mapq', 'std_mapq', 'nearest_splice_dist']
    features_to_normalize = ["total_reads", "read_start_density", "read_end_density", "soft_clip_mean", "soft_clip_max", "mean_mapq", "std_mapq", "strand_ratio", "coverage_before", "coverage_after", "delta_coverage", "nearest_splice_dist", "softclip_bias"]

    X_train, X_val, y_train, y_val = stratified_split(df, validation_chrom_file=val_chrom_file)
    numeric_cols = df.select_dtypes(include=["number"]).columns.tolist()
    numeric_cols = [col for col in numeric_cols if col not in drop]
    X_train = X_train[numeric_cols]
    X_val = X_val[numeric_cols]

    print(X_train.columns)
    if config['normalize']:
        print("Normalizing features...")
        scaler = StandardScaler()
        X_train[features_to_normalize] = scaler.fit_transform(X_train[features_to_normalize])
        X_val[features_to_normalize] = scaler.transform(X_val[features_to_normalize])

    model = load_model(model_type, config)
    model.fit(X_train, y_train)

    y_pred = model.predict(X_val)
    y_prob = model.predict_proba(X_val)[:, 1]

    model_dir = f"{out_dir}/models"
    report_dir = f"{out_dir}/reports/{site_type}"
    plot_dir = f"{out_dir}/plots/{site_type}"
    prediction_dir = f"{out_dir}/predictions/{site_type}"
    prediction_home = f"{out_dir}/predictions"
    os.makedirs(prediction_dir, exist_ok=True)
    os.makedirs(model_dir, exist_ok=True)
    os.makedirs(report_dir, exist_ok=True)
    os.makedirs(plot_dir, exist_ok=True)

    if model_type == "xgboost":
        model.save_model(f"{model_dir}/{site_type}_xgb_model.json")
    elif model_type == "lightgbm":
        model.booster_.save_model(f"{model_dir}/{site_type}_lgbm_model.txt")
    elif model_type == "randomforest":
        dump(model, os.path.join(model_dir, f"{site_type}_rf_model.joblib"))

    metrics = evaluate_model(y_val, y_pred, y_prob, plot_path=f"{plot_dir}/{tool}_{model_type}_pr_curve.png")

    if model_type == "randomforest":
        importances = model.feature_importances_
        sorted_idx = np.argsort(importances)[::-1]
        sorted_feats = np.array(numeric_cols)[sorted_idx]
        sorted_vals = importances[sorted_idx]

        # Save as CSV
        fi_df = pd.DataFrame({
            "feature": sorted_feats,
            "importance": sorted_vals
        })
        fi_df.to_csv(os.path.join(report_dir, f"{tool}_rf_feature_importance.csv"), index=False)


        plt.figure(figsize=(8, 6))
        plt.barh(sorted_feats[:20][::-1], sorted_vals[:20][::-1])  # top 20
        plt.xlabel("Importance")
        plt.title(f"Top Feature Importances ({site_type.upper()} - RF)")
        plt.tight_layout()
        
        os.makedirs(report_dir, exist_ok=True)
        plt.savefig(os.path.join(plot_dir, f"{tool}_feature_importance.png"))
        plt.close()

    with open(f"{report_dir}/{tool}_{model_type}_metrics_summary.txt", "w") as f:    
        for k, v in metrics.items():
            if k == "confusion_matrix":
                f.write(f"{k}:\n")
                f.write(f"  TP: {v[1][1]}\n")
                f.write(f"  FP: {v[0][1]}\n")
                f.write(f"  FN: {v[1][0]}\n")
                f.write(f"  TN: {v[0][0]}\n")
            else:
                f.write(f"{k}: {v:.4f}\n")

    # Save predictions as TSV
    pred_df = df.loc[y_val.index, ['chrom', 'position', 'label']].copy()
    pred_df["site_type"] = site_type.upper()
    pred_df["prediction"] = y_pred
    pred_df["probability"] = y_prob
    pred_df = pred_df[["site_type", "chrom", "position", "label", "prediction", "probability"]]
    pred_df["chrom"] = pred_df["chrom"].astype(str)
    pred_df["position"] = pred_df["position"].astype(int)
    pred_df["label"] = pred_df["label"].astype(int)
    pred_df["prediction"] = pred_df["prediction"].astype(int)
    pred_df["probability"] = pred_df["probability"].astype(float)

    pred_df.to_csv(os.path.join(prediction_dir, f"{tool}_{model_type}_predictions.tsv"), sep="\t", index=False)


    # save the predictions with features used for training and test set 
    X_df = df.loc[:, numeric_cols].copy()
    y_prob = model.predict_proba(X_df)[:, 1]
    X_df["site_type"] = site_type.upper()
    X_df["chrom"] = df["chrom"].astype(str)
    X_df["position"] = df["position"].astype(int)
    X_df["probability"] = y_prob
    X_df = X_df[ ["site_type", "chrom", "position"] + numeric_cols + ["probability"] ]

    print(f"X features df.shape: {X_df.shape}")

    X_df.to_csv(os.path.join(prediction_dir, f"{tool}_{model_type}_predictions_with_features.tsv"), sep="\t", index=False)
    return metrics, prediction_dir



def merge_tss_tes_predictions(pred_dir, tool, model_type):
    """
    Merge TSS and TES predictions into one TSV file if both exist.
    """
    tss_file = os.path.join(pred_dir, "tss", f"{tool}_{model_type}_predictions.tsv")
    tes_file = os.path.join(pred_dir, "tes", f"{tool}_{model_type}_predictions.tsv")
    merged_file = os.path.join(pred_dir, f"{tool}_{model_type}_merged.tsv")

    if os.path.exists(tss_file) and os.path.exists(tes_file):
        df_tss = pd.read_csv(tss_file, sep='\t')
        df_tes = pd.read_csv(tes_file, sep='\t')
        merged_df = pd.concat([df_tss, df_tes], ignore_index=True)
        merged_df.to_csv(merged_file, sep='\t', index=False)
        print(f"✅ [INFO] Merged TSS and TES predictions saved to {merged_file}")
    else:
        print(f"[INFO] TSS: {tss_file} or TES: {tes_file} prediction file missing. Merge skipped.")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--config", required=True)
    parser.add_argument("--site_type", required=True, choices=["tss", "tes"])
    parser.add_argument("--model_type", required=True, choices=["xgboost", "lightgbm", "randomforest"])
    parser.add_argument("--normalize", action='store_true', default=False, help="Normalize features")
    parser.add_argument("--out_dir", default="", help="Directory to save predictions")
    parser.add_argument("--val_chrom_file", default="", help="File containing validation chromosomes")
    args = parser.parse_args()

    # extract tool name
    tool = os.path.basename(args.input).split("_")[0]
    df = pd.read_csv(args.input)
    with open(args.config) as f:
        config = yaml.safe_load(f)

    config['normalize'] = args.normalize
    print("Configuration loaded:{}".format(config))

    metrics, prediction_dir = train_and_evaluate(df, args.model_type, config, args.site_type, tool, args.out_dir, args.val_chrom_file)

    print(f"✅ {args.site_type.upper()} [{args.model_type}] - F1: {metrics['f1']:.4f}, AUPR: {metrics['aupr']:.4f}")

    merge_tss_tes_predictions(pred_dir=f"{args.out_dir}/predictions", tool=tool, model_type=args.model_type)

if __name__ == "__main__":
    main()
    # Attempt to merge after each run
    

