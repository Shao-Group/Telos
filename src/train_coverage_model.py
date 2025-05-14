#!/usr/bin/env python3
import argparse
import pandas as pd
import numpy as np
import re
import joblib
from xgboost import XGBClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import (
    classification_report,
    roc_auc_score,
    average_precision_score,
    accuracy_score
)
from utils.ml_utils import chrom_to_int

def load_tmap_labels(tmap_path):
    df = pd.read_csv(tmap_path, sep='\t', comment='#', header=0)
    df = df[['qry_id', 'class_code']].rename(columns={'qry_id': 'transcript_id'})
    df['label'] = (df['class_code'] == '=').astype(int)
    return df[['transcript_id', 'label']]

def train_model(tool, model,
                tss_feat, tes_feat,
                cov_tsv, tmap_file,
                out_pred_file):
    # load data
    df_tss   = pd.read_csv(tss_feat, sep='\t')
    df_tes   = pd.read_csv(tes_feat, sep='\t')
    df_cov   = pd.read_csv(cov_tsv,    sep='\t')
    df_label = load_tmap_labels(tmap_file)

    print(f"\n=== [{tool}/{model}] merging TSS features ===")
    m1 = df_cov.merge(
        df_tss,
        left_on=['tss_chrom','tss_pos'],
        right_on=['chrom','position'],
        how='outer',
        indicator=True
    )
    unmatched_tss = m1[m1['_merge']=='left_only']
    if not unmatched_tss.empty:
        print("Unmatched coverage → no TSS features for:")
        unmatched_tss_df = unmatched_tss[['tss_chrom','tss_pos', 'transcript_id']].drop_duplicates()
        print(f"Number of unmatched tss: {unmatched_tss_df.shape[0]}")
        # print(unmatched_tss_df.to_string(index=False))
    df = m1[m1['_merge']=='both'].drop(columns=['_merge','chrom','position'])

    print(f"\n=== [{tool}/{model}] merging TES features ===")
    m2 = df.merge(
        df_tes,
        left_on=['tes_chrom','tes_pos'],
        right_on=['chrom','position'],
        how='outer',
        indicator=True,
        suffixes=('','_tes')
    )
    unmatched_tes = m2[m2['_merge']=='left_only']
    if not unmatched_tes.empty:
        print("Unmatched after TSS→ no TES features for:")
        unmatched_tes_df = unmatched_tes[['tes_chrom','tes_pos', 'transcript_id']].drop_duplicates()
        print(f"Number of unmatched tes: {unmatched_tes_df.shape[0]}")

    df = m2[m2['_merge']=='both'].drop(columns=['_merge','chrom','position'])

    print(f"\n=== [{tool}/{model}] merging labels ===")
    m3 = df.merge(
        df_label,
        on='transcript_id',
        how='outer',
        indicator=True
    )
    unmatched_lbl = m3[m3['_merge']=='left_only']
    if not unmatched_lbl.empty:
        print("Unmatched after TES→ no label for transcript_id:")
        print(unmatched_lbl[['transcript_id']].drop_duplicates().to_string(index=False))
    df = m3[m3['_merge']=='both'].drop(columns=['_merge'])

    print(f"\nFinal merged shape: {df.shape}")

    # split train/test by chromosome
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

    # train
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

    # predict & evaluate
    y_pred = clf.predict(X_test)
    y_prob = clf.predict_proba(X_test)[:,1]

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
    pred_df.to_csv(out_pred_file, sep='\t', index=False)

    # return metrics
    return {
        'tool': tool,
        'model': model,
        'accuracy': acc,
        'roc_auc': roc,
        'aupr': aupr,
        'f1_macro': report_dict['macro avg']['f1-score'],
        'precision_macro': report_dict['macro avg']['precision'],
        'recall_macro': report_dict['macro avg']['recall']
    }

def train_all_models(tools, models, data_home, pred_folder):
    metrics_list = []
    for tool in tools:
        for model in models:
            pred_file = f"{pred_folder}/transcripts/{tool}_{model}_transcript_predictions.tsv"
            tss_feat  = f"{pred_folder}/tss/{tool}_{model}_predictions_with_features.tsv"
            tes_feat  = f"{pred_folder}/tes/{tool}_{model}_predictions_with_features.tsv"
            cov_tsv   = f"{data_home}/{tool}-cov.tsv"
            tmap_file = f"{data_home}/{tool}.{tool}.gtf.tmap"
            m = train_model(tool, model,
                            tss_feat, tes_feat,
                            cov_tsv, tmap_file,
                            pred_file)
            metrics_list.append(m)
            # universe variant
            pred_file_u = f"{pred_folder}/transcripts/{tool}_universe_{model}_transcript_predictions.tsv"
            tss_feat_u  = f"{pred_folder}/tss/universe_{model}_predictions_with_features.tsv"
            tes_feat_u  = f"{pred_folder}/tes/universe_{model}_predictions_with_features.tsv"
            m = train_model(tool + "_universe", model,
                            tss_feat_u, tes_feat_u,
                            cov_tsv, tmap_file,
                            pred_file_u)
            metrics_list.append(m)

    metrics_out = f"{pred_folder}/metrics_summary_tr_model.csv"
    # save metrics summary
    dfm = pd.DataFrame(metrics_list)
    dfm.to_csv(metrics_out, index=False)
    print(f"\nSaved metrics summary to {metrics_out}")

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--data-home",   required=True)
    p.add_argument("--pred-folder", required=True)
    p.add_argument("--tools",       nargs='+', required=True)
    p.add_argument("--models",      nargs='+', required=True)
    p.add_argument("--metrics-out", required=True,
                   help="Path to save CSV summary of accuracy/ROC/AUPR")
    args = p.parse_args()
    train_all_models(
        tools=args.tools,
        models=args.models,
        data_home=args.data_home,
        pred_folder=args.pred_folder
    )
