import argparse
import yaml
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import joblib
from sklearn.preprocessing import StandardScaler
from config import Config, load_config
from ml_utils import stratified_split, evaluate_model, load_model, load_tmap_labels, chrom_to_int, load_saved_model
from xgboost import XGBClassifier
from sklearn.metrics import (
    accuracy_score, roc_auc_score, average_precision_score,
    classification_report
)
from label_candidates import load_selected_features
from sklearn.ensemble import RandomForestClassifier

def validate_stage2(df_tss, df_tes, project_config, model_type, pretrained_model):
    df_cov = pd.read_csv(project_config.cov_file, sep="\t")
    df_label = load_tmap_labels(project_config.tmap_file)


    df = df_cov.merge(df_tss, left_on=["tss_chrom", "tss_pos"], right_on=["chrom", "position"], how="inner")
    df.drop(columns=["chrom", "position"], inplace=True)
    df = df.merge(df_tes, left_on=["tes_chrom", "tes_pos"], right_on=["chrom", "position"], how="inner", suffixes=("_tss", "_tes"))
    df = df.merge(df_label, left_on=["transcript_id"], right_on=["transcript_id"], how="inner")

    # MINIMAL MODIFICATION: Add transcript-level features before splitting
    df['transcript_length'] = np.abs(df['tes_pos'] - df['tss_pos'])
    df['log_transcript_length'] = np.log1p(df['transcript_length'])
    df['tss_confidence'] = df.get('probability_tss', 0.5)
    df['tes_confidence'] = df.get('probability_tes', 0.5) 
    df['min_confidence'] = np.minimum(df['tss_confidence'], df['tes_confidence'])
    df['confidence_product'] = df['tss_confidence'] * df['tes_confidence']
    
    # Add feature interactions for better transcript prediction
    # Exon-related interactions
    if 'exon_count' in df.columns:
        df['exon_density'] = df['exon_count'] / df['transcript_length'].clip(lower=1)  # Exons per kb
        df['confidence_exon_interaction'] = df['confidence_product'] * df['exon_count']
        df['coverage_per_exon'] = df['coverage'] / df['exon_count'].clip(lower=1)
    
    if 'total_exon_length' in df.columns and 'exon_count' in df.columns:
        df['avg_exon_length'] = df['total_exon_length'] / df['exon_count'].clip(lower=1)
    
    if 'max_exon_length' in df.columns and 'min_exon_length' in df.columns:
        df['exon_length_ratio'] = df['max_exon_length'] / df['min_exon_length'].clip(lower=1)
    
    # Coverage-transcript interactions
    df['coverage_length_ratio'] = df['coverage'] / df['transcript_length'].clip(lower=1)
    df['confidence_coverage_interaction'] = df['confidence_product'] * df['coverage']
    
    # TSS/TES confidence interactions
    df['confidence_sum'] = df['tss_confidence'] + df['tes_confidence']
    df['confidence_diff'] = abs(df['tss_confidence'] - df['tes_confidence'])
    
    # Transcript structure features
    if 'first_exon_length' in df.columns and 'last_exon_length' in df.columns:
        df['terminal_exon_ratio'] = df['first_exon_length'] / df['last_exon_length'].clip(lower=1)
    
    if 'mean_exon_length' in df.columns and 'transcript_length' in df.columns:
        df['exon_efficiency'] = df['mean_exon_length'] / df['transcript_length'].clip(lower=1)
    
    # Add log transformations for skewed features
    skewed_features = ['coverage', 'total_reads', 'transcript_length']
    for feature in skewed_features:
        if feature in df.columns:
            df[f'log_{feature}'] = np.log1p(df[feature])
    
    # Add exon length entropy (measure of exon length diversity)
    if all(col in df.columns for col in ['max_exon_length', 'min_exon_length', 'mean_exon_length']):
        # Calculate coefficient of variation as a proxy for entropy
        exon_length_cv = (df['max_exon_length'] - df['min_exon_length']) / df['mean_exon_length'].clip(lower=1)
        df['exon_length_entropy'] = np.log1p(exon_length_cv)
    
    # Add additional exon diversity features
    if 'exon_count' in df.columns and 'total_exon_length' in df.columns:
        df['exon_length_std'] = np.sqrt(df['exon_length_variance']) if 'exon_length_variance' in df.columns else 0
        df['exon_length_skewness'] = (df['max_exon_length'] - df['mean_exon_length']) / df['exon_length_std'].clip(lower=1)
    

    X_train, X_test, y_train, y_test, train_mask, test_mask = stratified_split(df, validation_chrom_file=project_config.validation_chromosomes_file, train_chrom_file=project_config.train_chromosomes_file, return_mask=True)
        
    # MINIMAL MODIFICATION: Better feature selection for stage 2
    # MINIMAL MODIFICATION: Better feature selection for stage 2
    drop = [
        'chrom','position','chrom_tes','position_tes',
        'tss_chrom','tss_pos','tes_chrom','tes_pos',
        'site_type_tss','site_type_tes','strand', 'strand_tes',
        'ref_id','chrom_num','transcript_id','label', 
        'label_tss', 'label_tes'
    ]
    # Select base features + key TSS/TES features + transcript features
    base_features = [c for c in df.columns if c not in drop]
    features = [f for f in base_features if f in df.columns]  # Only keep existing features

    X_train = X_train[features]
    X_test  = X_test[features]

    print(f"Training stage 2 with {len(features)} features")
    print(f"Key features: {features[:10]}")


    # MINIMAL MODIFICATION: Optimize hyperparameters for transcript prediction
    clf_xgb = XGBClassifier(
            n_estimators=500,  # More trees
            max_depth=8,       # Deeper for interactions
            learning_rate=0.05, # Lower learning rate
            subsample=0.8,
            colsample_bytree=0.8,
            objective='binary:logistic',
            random_state=42,
            eval_metric='aucpr'
            # scale_pos_weight=len(y_train[y_train==0]) / max(len(y_train[y_train==1]), 1)  # Handle imbalance
        )
    clf_rf = RandomForestClassifier(
        n_estimators=500,
        max_depth=8,
        random_state=42
    )
    
    # Load the pre-trained pipeline
    clf = joblib.load(pretrained_model)

    y_pred = clf.predict(X_test)
    y_prob = clf.predict_proba(X_test)[:, 1]

    # save model (optional - usually not needed in validation script)
    # joblib.dump(clf, f"{project_config.models_output_dir}/xgboost_stage2_model.pkl")


    acc   = accuracy_score(y_test, y_pred)
    roc   = roc_auc_score(y_test, y_prob, average='macro')
    aupr  = average_precision_score(y_test, y_prob, average='macro')
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
    pred_df['transcript_id'] = df.loc[test_mask, 'transcript_id'].values
    pred_df['pred_prob']      = y_prob
    pred_df['pred_label']     = y_pred
    pred_df = pred_df[['transcript_id','pred_prob','pred_label']]

    pred_df.to_csv(os.path.join(project_config.predictions_output_dir, f"{model_type}_stage2_predictions_val.csv"),sep="\t" , index=False)


    y_pred_train = clf.predict(X_train)
    y_prob_train = clf.predict_proba(X_train)[:, 1]
    pred_df_train = X_train.copy()
    pred_df_train['transcript_id'] = df.loc[train_mask, 'transcript_id'].values
    pred_df_train['pred_prob']      = y_prob_train
    pred_df_train['pred_label']     = y_pred_train
    pred_df_train = pred_df_train[['transcript_id','pred_prob','pred_label']]
    pred_df_train.to_csv(os.path.join(project_config.predictions_output_dir, f"{model_type}_stage2_predictions_train.csv"),sep="\t" , index=False)



def validate_stage1(df, model_type, model_config, project_config, pretrained_model, site_type):
    drop = ["chrom", "position", "strand", "label"]

    X_train, X_val, y_train, y_val = stratified_split(df, validation_chrom_file=project_config.validation_chromosomes_file, train_chrom_file=project_config.train_chromosomes_file)
    
    # Use selected features if available, otherwise use all numeric features
    # from label_candidates import load_selected_features
    selected_features = load_selected_features(project_config, site_type)
    
    if selected_features:
        # Use only selected features that exist in the DataFrame
        available_selected_features = [col for col in selected_features if col in df.columns]
        print(f"Using {len(available_selected_features)} selected features for {site_type}")
        numeric_cols = available_selected_features
    else:
        # Fallback to all numeric features
        raise ValueError(f"No selected features found for {site_type}, using all numeric features")
        # numeric_cols = df.select_dtypes(include=["number"]).columns.tolist()
        # numeric_cols = [col for col in numeric_cols if col not in drop]
    
    X_train = X_train[numeric_cols]
    X_val = X_val[numeric_cols]

    # model = load_saved_model(model_type, pretrained_model, model_config)
    model = joblib.load(pretrained_model)

    y_pred = model.predict(X_val)
    y_prob = model.predict_proba(X_val)[:, 1]

    
    # if model_type == "xgboost":
    #     model.save_model(f"{project_config.models_output_dir}/{site_type}_xgb_model.json")
    # elif model_type == "randomforest":
    #     dump(model, os.path.join(project_config.models_output_dir, f"{site_type}_randomforest_model.joblib"))

    pr_data_path = os.path.join(project_config.pr_data_dir, f"{site_type}_{model_type}_pr_data.csv")
    metrics = evaluate_model(y_val, y_pred, y_prob, pr_data_path )

    # save metrics to file
    metrics_list = []
    for metric_name, metric_value in metrics.items():
        if metric_name != 'confusion_matrix':  # Handle confusion matrix separately
            metrics_list.append({'metric': metric_name, 'value': metric_value})
    
    # Add confusion matrix as a formatted string
    cm_str = str(metrics['confusion_matrix'])
    metrics_list.append({'metric': 'confusion_matrix', 'value': cm_str})
    
    metrics_df = pd.DataFrame(metrics_list)
    metrics_file = os.path.join(project_config.metrics_output_dir, f"{site_type}_{model_type}_metrics.csv")
    metrics_df.to_csv(metrics_file, index=False)
    print(f"Metrics saved to {metrics_file}")

    if model_type == "randomforest":
        importances = model.named_steps['clf'].feature_importances_
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
    project_config.train_chromosomes_file = os.path.join(project_config.data_output_dir, "train_chromosomes.txt")

    # Check if selected feature files exist, otherwise use original labeled files
    
    print(f"Using original: {project_config.tss_labeled_file}")
    df_tss = pd.read_csv(project_config.tss_labeled_file, dtype={"chrom": str})
    print(f"df_tss.columns: {df_tss.columns}")
        
    with open(os.path.join(args.model_config_folder, f"{args.model_type}_config.yaml")) as f:
        model_config = yaml.safe_load(f)

    metrics, model, df_tss_pred = validate_stage1(df_tss, args.model_type, model_config, project_config, args.pretrained_tss_model, "tss")

    # Check if selected feature files exist for TES

    print(f"Selected features file not found, using original: {project_config.tes_labeled_file}")
    df_tes = pd.read_csv(project_config.tes_labeled_file, dtype={"chrom": str})
    print(f"df_tes.columns: {df_tes.columns}")
        
    with open(os.path.join(args.model_config_folder, f"{args.model_type}_config.yaml")) as f:
        model_config = yaml.safe_load(f)

    metrics, model, df_tes_pred = validate_stage1(df_tes, args.model_type, model_config, project_config, args.pretrained_tes_model, "tes")

    validate_stage2(df_tss_pred, df_tes_pred, project_config, args.model_type, args.pretrained_stage2_model)

    project_config.save_to_file(args.project_config)


if __name__ == "__main__":
    main()
    

