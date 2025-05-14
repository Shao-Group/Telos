#!/usr/bin/env python3
import argparse
import pandas as pd
import numpy as np
import re
import joblib
from xgboost import XGBClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report, roc_auc_score
from utils.ml_utils import chrom_to_int

def load_tmap_labels(tmap_path):
    # read in with header row (skip any comment lines)
    df = pd.read_csv(tmap_path, sep='\t', comment='#', header=0)

    # keep only the query‐ID and class code
    df = df[['qry_id', 'class_code']]

    # rename qry_id → transcript_id for consistency
    df = df.rename(columns={'qry_id': 'transcript_id'})

    # exact‐match '=' is positive (1), all others negative (0)
    df['label'] = (df['class_code'] == '=').astype(int)
    return df[['transcript_id', 'label']]


def train_model(tss_feat, tes_feat, cov_tsv, tmap_file, out_pred_file):
    # load
    df_tss   = pd.read_csv(tss_feat, sep='\t')
    df_tes   = pd.read_csv(tes_feat, sep='\t')
    df_cov   = pd.read_csv(cov_tsv, sep='\t')
    df_label = load_tmap_labels(tmap_file)

    print(f"Merging TSS, TES, and coverage data...")
    # print(f"TSS features columns: {df_tss.columns}")
    # print(f"DF coverage columns: {df_cov.columns}")

    # merge TSS
    df = df_cov.merge(
        df_tss,
        left_on=['tss_chrom','tss_pos'],
        right_on=['chrom',   'position'],
        how='inner'
    )

    print(f"After merging TSS: {df.shape}")
    # merge TES
    df = df.merge(
        df_tes,
        left_on=['tes_chrom','tes_pos'],
        right_on=['chrom',    'position'],
        how='inner',
        suffixes=('','_tes')
    )
    print(f"After merging TES: {df.shape}")

    # merge labels
    df = df.merge(df_label, on='transcript_id', how='inner')

    print(f"After merging labels: {df.shape}")

    # mark chrom number 
    df['chrom_num'] = df['tss_chrom'].apply(chrom_to_int)

    # split by chromosome
    train_mask = df['chrom_num'].between(1, 15)
    val_mask = ~train_mask
    
    train_df = df[train_mask]
    test_df  = df[val_mask]

    # prepare X/y
    drop_cols = ['chrom','position','chrom_tes','position_tes','tss_chrom', 'tss_pos', 'tes_chrom', 'tes_pos',
                 'site_type', 'site_type_tes',
                 'ref_id','class_code','chrom_num', 'transcript_id','label']
    features = [c for c in df.columns 
                if c not in drop_cols and not c.startswith('tss_') and not c.startswith('tes_')]
    X_train, y_train = train_df[features], train_df['label']
    X_test,  y_test  = test_df[features],  test_df['label']

    # train classifier
    clf = RandomForestClassifier(n_estimators=400, random_state=42)
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
    clf.fit(X_train, y_train)

    # eval on test set
    y_pred = clf.predict(X_test)
    y_prob = clf.predict_proba(X_test)[:,1]

    print("=== Test set classification report ===")
    print(classification_report(y_test, y_pred, digits=4))
    print(f"ROC AUC on test set: {roc_auc_score(y_test, y_prob):.4f}")
    print(f"AUPR on test set: {roc_auc_score(y_test, y_prob):.4f}")

    # save predicted transcript probabilities
    pred_df = test_df.copy()
    pred_df['pred_prob'] = y_prob
    pred_df['pred_label'] = y_pred
    pred_df = pred_df[['transcript_id', 'pred_prob', 'pred_label']]
    pred_df.to_csv(f"{out_pred_file}", sep='\t', index=False)
    
    # joblib.dump(clf, output_model)
    # print(f"Model saved to {output_model}")

def train_all_models(tools, models, data_home, pred_folder):
    """
    Train all models for all tools.
    """
    # tools.append("universe")
    for tool in tools:
        for model in models:
            pred_dir = f"{pred_folder}/transcripts/{tool}_{model}_transcript_predictions.tsv"
            tss_feat = f"{pred_folder}/tss/{tool}_{model}_predictions_with_features.tsv"
            tes_feat = f"{pred_folder}/tes/{tool}_{model}_predictions_with_features.tsv"
            cov_tsv = f"{data_home}/{tool}-cov.tsv"
            tmap_file = f"{data_home}/{tool}.{tool}.gtf.tmap"
            train_model(tss_feat, tes_feat, cov_tsv, tmap_file, pred_dir)
            
            pred_dir_universe = f"{pred_folder}/transcripts/{tool}_universe_{model}_transcript_predictions.tsv"
            tss_feat_universe = f"{pred_folder}/tss/universe_{model}_predictions_with_features.tsv"
            tes_feat_universe = f"{pred_folder}/tes/universe_{model}_predictions_with_features.tsv"
            train_model(tss_feat_universe, tes_feat_universe, cov_tsv, tmap_file, pred_dir_universe)

if __name__ == "__main__":
    p = argparse.ArgumentParser(
        description="Train classifier with chr1–15 as train, others as test"
    )
    p.add_argument("--tss-features", required=True)
    p.add_argument("--tes-features", required=True)
    p.add_argument("--cov-tsv",      required=True)
    p.add_argument("--tmap-file",    required=True)
    p.add_argument("--pred-folder", required=True)
    args = p.parse_args()
    train_model(args.tss_features, args.tes_features,
         args.cov_tsv, args.tmap_file,
         args.pred_folder)
