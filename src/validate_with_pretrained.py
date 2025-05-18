import argparse
import yaml
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from joblib import dump
from sklearn.preprocessing import StandardScaler
from config import Config, load_config
from ml_utils import stratified_split, evaluate_model, load_model, load_tmap_labels, chrom_to_int, load_saved_model
from xgboost import XGBClassifier
from sklearn.metrics import (
    accuracy_score, roc_auc_score, average_precision_score,
    classification_report
)

def validate_stage2(df_tss, df_tes, project_config, model_type, pretrained_model):
    df_cov = pd.read_csv(project_config.cov_file, sep="\t")
    df_label = load_tmap_labels(project_config.tmap_file)

    df = df_cov.merge(df_tss, left_on=["tss_chrom", "tss_pos"], right_on=["chrom", "position"], how="inner")
    df = df.merge(df_tes, left_on=["tes_chrom", "tes_pos"], right_on=["chrom", "position"], how="inner", suffixes=("", "_tes"))
    df = df.merge(df_label, left_on=["transcript_id"], right_on=["transcript_id"], how="inner")

    df['chrom_num'] = df['tss_chrom'].apply(chrom_to_int)
    train_mask = df['chrom_num'].between(1, 5)
    X_train = df[train_mask]
    X_test  = df[~train_mask]
    y_train = X_train['label']
    y_test  = X_test['label']

    # select feature columns
    drop = [
        'chrom','position','chrom_tes','position_tes',
        'tss_chrom','tss_pos','tes_chrom','tes_pos',
        'site_type','site_type_tes',
        'ref_id','chrom_num','transcript_id','label'
    ]
    features = [c for c in df.columns
                if c not in drop and not c.startswith('tss_') and not c.startswith('tes_')]
    X_train = X_train[features]
    X_test  = X_test[features]

    clf = XGBClassifier(
            n_estimators=400,
            max_depth=6,
            learning_rate=0.1,
            subsample=0.8,
            colsample_bytree=0.8,
            objective='binary:logistic',
            random_state=42,
            eval_metric='aucpr'
        )
    clf.load_model(pretrained_model)

    y_pred = clf.predict(X_test)
    y_prob = clf.predict_proba(X_test)[:, 1]

    # save model
    clf.save_model(f"{project_config.models_output_dir}/xgboost_stage2_model.json")


    acc   = accuracy_score(y_test, y_pred)
    roc   = roc_auc_score(y_test, y_prob)
    aupr  = average_precision_score(y_test, y_prob)
    report_dict = classification_report(y_test, y_pred, digits=4, output_dict=True)

    print(f"\nAccuracy: {acc:.4f}")
    print(f"ROC AUC:  {roc:.4f}")
    print(f"AUPR:     {aupr:.4f}")
    print("\nClassification report:")
    print(classification_report(y_test, y_pred, digits=4))

    print("F1 score (macro):", report_dict['macro avg']['f1-score'])
    print("Precision (macro):", report_dict['macro avg']['precision'])
    print("Recall (macro):", report_dict['macro avg']['recall'])

    # save predictions
    pred_df = X_test.copy()
    pred_df['transcript_id'] = df.loc[~train_mask, 'transcript_id'].values
    pred_df['pred_prob']      = y_prob
    pred_df['pred_label']     = y_pred
    pred_df = pred_df[['transcript_id','pred_prob','pred_label']]

    pred_df.to_csv(os.path.join(project_config.predictions_output_dir, f"{model_type}_stage2_predictions.csv"),sep="\t" , index=False)


def validate_stage1(df, model_type, model_config, project_config, pretrained_model, site_type):
    drop = ["chrom", "position", "strand", "label" , "soft_clip_entropy"]

    X_train, X_val, y_train, y_val = stratified_split(df, validation_chrom_file=project_config.validation_chromosomes_file)
    numeric_cols = df.select_dtypes(include=["number"]).columns.tolist()
    numeric_cols = [col for col in numeric_cols if col not in drop]
    X_train = X_train[numeric_cols]
    X_val = X_val[numeric_cols]

    model = load_saved_model(model_type, pretrained_model, model_config)

    y_pred = model.predict(X_val)
    y_prob = model.predict_proba(X_val)[:, 1]

    
    if model_type == "xgboost":
        model.save_model(f"{project_config.models_output_dir}/{site_type}_xgb_model.json")
    elif model_type == "randomforest":
        dump(model, os.path.join(project_config.models_output_dir, f"{site_type}_randomforest_model.joblib"))

    pr_data_path = os.path.join(project_config.pr_data_dir, f"{site_type}_{model_type}_pr_data.csv")
    metrics = evaluate_model(y_val, y_pred, y_prob, pr_data_path )

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
        fi_df.to_csv(os.path.join(project_config.feature_importance_dir, f"{site_type}_rf_feature_importance.csv"), index=False)

    X_df = df.loc[:, numeric_cols].copy()
    y_prob = model.predict_proba(X_df)[:, 1]
    X_df["site_type"] = site_type.upper()
    X_df["chrom"] = df["chrom"].astype(str)
    X_df["position"] = df["position"].astype(int)
    X_df["probability"] = y_prob
    X_df = X_df[ ["site_type", "chrom", "position"] + numeric_cols + ["probability"] ]

    X_df.to_csv(os.path.join(project_config.predictions_output_dir, f"{site_type}_{model_type}_predictions.csv"), index=False)
    return metrics, model, X_df



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-config", required=True)
    parser.add_argument("--model-config-folder", required=True, help="Path to model config folder")
    parser.add_argument("--pretrained_tss_model", required=True, help="Path to pretrained TSS model")
    parser.add_argument("--pretrained_tes_model", required=True, help="Path to pretrained TES model")
    parser.add_argument("--pretrained_stage2_model", required=True, help="Path to pretrained stage 2 model")
    parser.add_argument("--model_type", required=True, help="Model type to use for Stage 1")
    args = parser.parse_args()

    project_config = load_config(args.project_config)
    project_config.validation_chromosomes_file = os.path.join(project_config.data_output_dir, "validation_chromosomes.txt")

    df_tss = pd.read_csv(project_config.tss_labeled_file)
    with open(os.path.join(args.model_config_folder, f"tss_{args.model_type}_config.yaml")) as f:
        model_config = yaml.safe_load(f)

    metrics, model, df_tss_pred = validate_stage1(df_tss, args.model_type, model_config, project_config, args.pretrained_tss_model, "tss")

    df_tes = pd.read_csv(project_config.tes_labeled_file)
    with open(os.path.join(args.model_config_folder, f"tes_{args.model_type}_config.yaml")) as f:
        model_config = yaml.safe_load(f)

    metrics, model, df_tes_pred = validate_stage1(df_tes, args.model_type, model_config, project_config, args.pretrained_tes_model, "tes")

    validate_stage2(df_tss_pred, df_tes_pred, project_config, args.model_type, args.pretrained_stage2_model)

    project_config.save_to_file(args.project_config)


if __name__ == "__main__":
    main()
    

