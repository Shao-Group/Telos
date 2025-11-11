import numpy as np
import argparse
import pandas as pd
import os
from tqdm import tqdm
from ml_utils import normalize_chrom
from config import load_config
import sys
from sklearn.feature_selection import SelectKBest, mutual_info_classif
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import RFECV
import glob

def load_candidates(candidate_path):
    """Load extracted candidate features."""
    return pd.read_csv(candidate_path, dtype={"chrom": str, "position": int, "strand": str})

def load_reference(reference_path):
    """Load reference annotation TSS/TES positions."""
    ref_df = pd.read_csv(reference_path, sep=" ", header=None,
                         names=["site_type", "chrom", "position", "pos_strand_cnt", "neg_strand_cnt"],
                         dtype={"chrom": str})
    ref_df["strand"] = ref_df.apply(
        lambda x: "+" if x["pos_strand_cnt"] > x["neg_strand_cnt"] else "-", axis=1
    )
    ref_df = ref_df[["site_type", "chrom", "position", "strand"]]
    # print(ref_df.head(5))
    return ref_df



def label_candidates(candidate_df, reference_df, site_type, max_distance=50):
    """Fast labeling using grouped and vectorized search."""
    
    reference_df = reference_df[
        reference_df['site_type'].str.lower() == site_type.lower()
    ].copy()
    reference_df['chrom'] = reference_df['chrom'].astype(str)
    reference_df['position'] = reference_df['position'].astype(int)

    # Group reference sites for fast lookup
    grouped_ref = reference_df.groupby(['chrom', 'strand'])

    labels = np.zeros(len(candidate_df), dtype=int)

    print(f"Labeling {site_type.upper()} candidates with {max_distance} nt tolerance ...")
    for i, row in candidate_df.iterrows():
        chrom, pos, strand = row['chrom'], row['position'], row['strand']
        if i % 1000 == 0:
            print(f"Labeled {i} candidates")
            sys.stdout.flush()
            
        try:
            group = grouped_ref.get_group((chrom, strand))
        except KeyError:
            continue  # No matching chrom+strand in reference
        
        # Use numpy for fast range check
        ref_positions = group['position'].values
        match_found = np.any(np.abs(ref_positions - pos) <= max_distance)
        labels[i] = 1 if match_found else 0

    candidate_df['label'] = labels

    print(f"Labeled {sum(labels)} {site_type.upper()} candidates")
    print(f"Positive {site_type.upper()} candidates: {sum(labels == 1)}")
    print(f"Negative {site_type.upper()} candidates: {sum(labels == 0)}")

    return candidate_df


    """
    Label individual soft-clipped sequences based on their candidate site positions.
    
    Args:
        softclip_df: DataFrame with soft-clipped sequences
        reference_df: DataFrame with reference annotations
        site_type: 'TSS' or 'TES'
        max_distance: Maximum distance for matching
        
    Returns:
        DataFrame with labels added
    """
    reference_df = reference_df[
        reference_df['site_type'].str.lower() == site_type.lower()
    ].copy()
    reference_df['chrom'] = reference_df['chrom'].astype(str)
    reference_df['position'] = reference_df['position'].astype(int)

    # Group reference sites for fast lookup
    grouped_ref = reference_df.groupby(['chrom', 'strand'])
    grouped_softclip_df = softclip_df.groupby(['chrom', 'strand', 'position'])

    labels = np.zeros(len(softclip_df), dtype=int)
    softclip_df['label'] = labels

    print(f"Labeling {site_type.upper()} soft-clipped sequences with {max_distance} nt tolerance ...")
    i = 0
    for key, row in grouped_softclip_df:
        if i % 1000 == 0:
            print(f"Labeled {i} soft-clipped sequences")
            sys.stdout.flush()
        i += 1

        chrom, strand, pos = key
        idx = grouped_softclip_df.indices[key]
        
        try:
            group = grouped_ref.get_group((chrom, strand))
        except KeyError:
            continue  # No matching chrom+strand in reference

        # Use numpy for fast range check
        ref_positions = group['position'].values
        match_found = np.any(np.abs(ref_positions - pos) <= max_distance)
        softclip_df.loc[idx, 'label'] = 1 if match_found else 0

    # softclip_df['label'] = labels
    return softclip_df

def perform_feature_selection(df, cfg, site_type):
    """
    Perform comprehensive feature selection for a given site type.
    
    Args:
        df: DataFrame with features and labels
        cfg: Configuration object
        site_type: 'TSS' or 'TES'
    
    Returns:
        selected_features: List of selected feature names
        feature_importance: Dictionary mapping features to importance scores
    """
    print(f"Performing feature selection for {site_type}...")
    
    # Separate features and labels
    exclude_cols = ['chrom', 'position', 'strand', 'label']
    feature_cols = [col for col in df.columns if col not in exclude_cols]
    
    df["chrom"] = df["chrom"].astype(str)
    df = df.copy()
    df["chrom"] = df["chrom"].apply(lambda x: x if x.startswith('chr') else f'chr{x}')
    df = df[df["chrom"].isin(["chr1"])]

    X = df[feature_cols]
    y = df['label']
    
    print(f"Total features before selection: {len(feature_cols)}")
    kept_features = feature_cols
    X_filtered = X

    # Step 1: Remove features with low variance
    # X_filtered, kept_features = _remove_low_variance_features(X, threshold=0.01)
    # print(f"Features after variance filtering: {len(kept_features)}")
    
    # # Step 2: Remove highly correlated features
    # X_filtered, kept_features = _remove_correlated_features(X_filtered, threshold=0.95)
    # print(f"Features after correlation filtering: {len(kept_features)}")
    
    # # Step 3: Univariate feature selection
    # X_filtered, kept_features, univariate_scores = _univariate_feature_selection(
    #     X_filtered, y, k_features=min(500, len(kept_features))
    # )
    # print(f"Features after univariate selection: {len(kept_features)}")
    
    # Step 4: Recursive feature elimination with Random Forest
    if len(kept_features) > 50000:  # Only do RFE if we have enough features
        X_filtered, kept_features, rf_importance = _recursive_feature_elimination(
            X_filtered, y, max_features=min(10, len(kept_features))
        )
        print(f"Features after recursive elimination: {len(kept_features)}")
    else:
        # Use Random Forest to get feature importance for remaining features
        rf = RandomForestClassifier(n_estimators=100, random_state=42, n_jobs=-1)
        rf.fit(X_filtered, y)
        rf_importance = dict(zip(kept_features, rf.feature_importances_))
    
    # Step 5: Final selection based on biological relevance and importance
    final_features = _final_biological_selection(kept_features, rf_importance, site_type)
    print(f"Final selected features: {len(final_features)}")
    
    # Save selected features to config directory
    _save_selected_features(cfg, final_features, rf_importance, site_type)
    
    return final_features, rf_importance



def _recursive_feature_elimination(X, y, max_features=200):
    """Perform recursive feature elimination with Random Forest."""
    rf = RandomForestClassifier(n_estimators=50, random_state=42, n_jobs=-1)
    
    # Use RFECV to find optimal number of features
    rfecv = RFECV(
        estimator=rf,
        step=1,
        cv=3,
        scoring='average_precision',  # Use AUPR for imbalanced data
        n_jobs=-1,
        min_features_to_select=min(50, max_features)
    )
    
    rfecv.fit(X, y)
    
    # Get selected features
    selected_mask = rfecv.support_
    kept_features = X.columns[selected_mask].tolist()
    
    # Get feature importance from the final estimator
    feature_importance = dict(zip(kept_features, rfecv.estimator_.feature_importances_))
    
    return X.loc[:, selected_mask], kept_features, feature_importance


def _final_biological_selection(features, importance_scores, site_type):
    """Apply biological knowledge to make final feature selection."""
    
    # Define features to drop
    features_to_drop = set()
    
    # 1. Drop redundant features
    redundant_features = [
        'gc_content_basic',  # keep gc_content
        'at_content_basic',  # can be calculated as 1 - gc_content
        # 'start_soft_clip_count', 'end_soft_clip_count',  # keep means
        # 'avg_clip_length',  # keep max_clip_length
        'total_reads',  # redundant with coverage features
        'total_clip_length'  # redundant with other clip features
    ]
    features_to_drop.update(redundant_features)
    
    # 2. Drop all 5prime, 3prime k-mer positional features
    positional_kmer_patterns = [
        '5prime', '3prime', 'middle_ratio'
    ]
    for feature in features:
        if any(pattern in feature for pattern in positional_kmer_patterns):
            features_to_drop.add(feature)
    
    # 3. Drop biological motif features
    biological_motif_patterns = [
        'cleavage_sites', 'downstream_elements', 'u_rich_regions',
        'pyrimidine_rich', 'polya_signals', 'degenerate_',
        '_total', '_terminus', '_terminal_ratio'
    ]
    for feature in features:
        if any(pattern in feature for pattern in biological_motif_patterns):
            features_to_drop.add(feature)
    
    # 4. Drop terminal enrichment features
    terminal_enrichment_patterns = ['_terminal_enrichment']
    for feature in features:
        if any(pattern in feature for pattern in terminal_enrichment_patterns):
            features_to_drop.add(feature)
    
    # 5. Keep fixed k-mers (kmer_XXX pattern) - don't drop these
    # 6. Keep homopolymer features (max_polyX) - don't drop these
    # 7. Keep soft-clip medians - don't drop these
    
    # Filter out features to drop
    filtered_features = [f for f in features if f not in features_to_drop]
    
    # sorted_features = sorted(filtered_features, key=lambda x: importance_scores.get(x, 0), reverse=True)
    # positive_features = [f for f in filtered_features if importance_scores.get(f, 0) > 0]
    # final_features = positive_features
    final_features = filtered_features

    print(f"Feature filtering for {site_type}:")
    print(f"  Original features: {len(features)}")
    print(f"  Features dropped: {len(features_to_drop)}")
    print(f"  Features kept: {len(filtered_features)}")
    
    # Sort by importance scores and return top 10
    # sorted_features = sorted(filtered_features, key=lambda x: importance_scores.get(x, 0), reverse=True)
    # final_features = filtered_features
    
    return final_features


def _save_selected_features(cfg, selected_features, importance_scores, site_type):
    """Save selected features and their importance scores."""
    # Use config paths for selected feature files
    if site_type.upper() == 'TSS':
        features_file = cfg.tss_selected_feature_file
    else:  # TES
        features_file = cfg.tes_selected_feature_file
    
    # Save selected features list
    with open(features_file, 'w') as f:
        for feat in selected_features:
            f.write(f"{feat}\n")
    
    # Save feature importance scores in the same directory
    features_dir = os.path.dirname(features_file)
    importance_file = os.path.join(features_dir, f'{cfg.data_name}_{site_type.lower()}_feature_importance.csv')
    importance_df = pd.DataFrame([
        {'feature': feat, 'importance': importance_scores.get(feat, 0)}
        for feat in selected_features
    ])
    importance_df = importance_df.sort_values('importance', ascending=False)
    importance_df.to_csv(importance_file, index=False)
    
    print(f"Saved selected features to {features_file}")
    print(f"Saved feature importance to {importance_file}")


def load_selected_features(cfg, site_type):
    """Load previously selected features."""
    # Use config paths for selected feature files
    if site_type.upper() == 'TSS':
        features_file = cfg.tss_selected_feature_file
    else:  # TES
        features_file = cfg.tes_selected_feature_file
    
    if os.path.exists(features_file):
        with open(features_file, 'r') as f:
            selected_features = [line.strip() for line in f.readlines()]
        print(f"Loaded {len(selected_features)} selected features for {site_type} from {features_file}")
        return selected_features
    else:
        print(f"No selected features file found for {site_type} at {features_file}")
        return None


def load_selected_features_from_pretrained(pretrained_model_dir, site_type):
    """
    Load selected features from a pretrained model directory (for cross-validation).
    
    Args:
        pretrained_model_dir: Path to the pretrained model directory
        site_type: 'TSS' or 'TES'
    
    Returns:
        List of selected feature names or None if not found
    """
    # Look for selected features file in the pretrained model's features directory
    # Navigate from models dir to features dir: models/../features/
    features_dir = os.path.join(os.path.dirname(pretrained_model_dir), 'features')
    
    # Try to find the selected features file with different naming patterns
    possible_patterns = [
        f"*_{site_type.lower()}_selected_features.txt",  # dataset_name_tss_selected_features.txt
        f"{site_type.lower()}_selected_features.txt"     # tss_selected_features.txt
    ]
    
    features_file = None
    
    for pattern in possible_patterns:
        search_path = os.path.join(features_dir, pattern)
        matches = glob.glob(search_path)
        if matches:
            features_file = matches[0]  # Take the first match
            break
    
    if features_file and os.path.exists(features_file):
        with open(features_file, 'r') as f:
            selected_features = [line.strip() for line in f.readlines()]
        print(f"Loaded {len(selected_features)} selected features for {site_type} from pretrained model: {features_file}")
        return selected_features
    else:
        print(f"No selected features file found for {site_type} in pretrained model directory: {features_dir}")
        return None
    

def main(args, cfg):
    # Load configuration
    reference_df = load_reference(cfg.ref_candidate_file)
    
    # Label feature candidates if they exist
    if os.path.exists(cfg.tss_feature_file) and os.path.exists(cfg.tes_feature_file):
        # Load candidates
        tss_candidate_df = load_candidates(cfg.tss_feature_file)
        tes_candidate_df = load_candidates(cfg.tes_feature_file)

        # Label candidates
        tss_labeled_df = label_candidates(tss_candidate_df, reference_df, "TSS", args.distance)
        tes_labeled_df = label_candidates(tes_candidate_df, reference_df, "TES", args.distance)

        # Save labeled candidates
        tss_labeled_df.to_csv(cfg.tss_labeled_file, index=False)
        tes_labeled_df.to_csv(cfg.tes_labeled_file, index=False)
        print(f"Labeled TSS candidates saved to: {cfg.tss_labeled_file}")
        print(f"Labeled TES candidates saved to: {cfg.tes_labeled_file}")

        # Perform feature selection on extracted features
        # Only do feature selection during training (not testing)
        if getattr(cfg, 'is_training', True):  # Default to True for backward compatibility
            
            # Load labeled data and perform feature selection for TSS
            if os.path.exists(cfg.tss_labeled_file):
                print("Performing feature selection for TSS...")
                tss_df = pd.read_csv(cfg.tss_labeled_file)
                if 'label' in tss_df.columns:
                    selected_tss_features, tss_importance = perform_feature_selection(tss_df, cfg, 'TSS')
                    
                    # Filter and save TSS features with selected features only
                    tss_filtered = tss_df[['chrom', 'position', 'strand'] + selected_tss_features + ['label']]
                    # tss_filtered.to_csv(cfg.tss_feature_file.replace('.csv', '_selected.csv'), index=False)
                    # print(f"Saved filtered TSS features to {cfg.tss_feature_file.replace('.csv', '_selected.csv')}")
                
            
            # Load labeled data and perform feature selection for TES
            if os.path.exists(cfg.tes_feature_file):
                print("Performing feature selection for TES...")
                tes_df = pd.read_csv(cfg.tes_labeled_file)
                if 'label' in tes_df.columns:
                    selected_tes_features, tes_importance = perform_feature_selection(tes_df, cfg, 'TES')
                    
                    # Filter and save TES features with selected features only
                    tes_filtered = tes_df[['chrom', 'position', 'strand'] + selected_tes_features + ['label']]
                    # tes_filtered.to_csv(cfg.tes_labeled_file.replace('.csv', '_selected.csv'), index=False)
                    # print(f"Saved filtered TES features to {cfg.tes_labeled_file.replace('.csv', '_selected.csv')}")
        
        # else:
        #     # For testing, load previously selected features and filter
        #     print("Loading previously selected features for testing...")
            
        #     # Filter TSS features using previously selected features
        #     if os.path.exists(cfg.tss_labeled_file):
        #         selected_tss_features = load_selected_features(cfg, 'TSS')
        #         # print(f"Selected TSS features: {selected_tss_features}")
        #         if selected_tss_features:
        #             tss_df = pd.read_csv(cfg.tss_labeled_file)
        #             available_features = [f for f in selected_tss_features if f in tss_df.columns]
        #             tss_filtered = tss_df[['chrom', 'position', 'strand'] + available_features]
        #             tss_filtered.to_csv(cfg.tss_labeled_file.replace('.csv', '_selected.csv'), index=False)
        #             print(f"Filtered TSS features using {len(available_features)} selected features")
            
        #     # Filter TES features using previously selected features
        #     if os.path.exists(cfg.tes_labeled_file):
        #         selected_tes_features = load_selected_features(cfg, 'TES')
        #         if selected_tes_features:
        #             tes_df = pd.read_csv(cfg.tes_labeled_file)
        #             available_features = [f for f in selected_tes_features if f in tes_df.columns]
        #             tes_filtered = tes_df[['chrom', 'position', 'strand'] + available_features]
        #             tes_filtered.to_csv(cfg.tes_labeled_file.replace('.csv', '_selected.csv'), index=False)
        #             print(f"Filtered TES features using {len(available_features)} selected features")
    
    # Label soft-clipped sequences if they exist
    # if hasattr(cfg, 'tss_softclip_file') and hasattr(cfg, 'tes_softclip_file'):
    #     if os.path.exists(cfg.tss_softclip_file) and os.path.exists(cfg.tes_softclip_file):
    #         # Load soft-clipped sequences
    #         tss_softclip_df = pd.read_csv(cfg.tss_softclip_file, dtype={"chrom": str})
    #         tes_softclip_df = pd.read_csv(cfg.tes_softclip_file, dtype={"chrom": str})
            
    #         # Label soft-clipped sequences
    #         tss_softclip_labeled_df = label_softclip_sequences(tss_softclip_df, reference_df, "TSS", args.distance)
    #         tes_softclip_labeled_df = label_softclip_sequences(tes_softclip_df, reference_df, "TES", args.distance)
            
    #         # Save labeled soft-clipped sequences
    #         tss_softclip_labeled_df.to_csv(cfg.tss_softclip_labeled_file, index=False)
    #         tes_softclip_labeled_df.to_csv(cfg.tes_softclip_labeled_file, index=False)
    #         print(f"Labeled TSS soft-clipped sequences saved to: {cfg.tss_softclip_labeled_file}")
    #         print(f"Labeled TES soft-clipped sequences saved to: {cfg.tes_softclip_labeled_file}")
            
    #         # Print summary statistics
    #         print(f"\nSoft-clipped sequence labeling summary:")
    #         print(f"TSS sequences: {len(tss_softclip_labeled_df)} total, {sum(tss_softclip_labeled_df['label'])} positive")
    #         print(f"TES sequences: {len(tes_softclip_labeled_df)} total, {sum(tes_softclip_labeled_df['label'])} positive")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Label candidate TSS/TES sites based on reference.")
    parser.add_argument('-d', '--distance', type=int, default=50, help="Maximum distance allowed for matching.")
    parser.add_argument('-c', '--config', type=str, required=True, help="Path to the configuration file.")
    args = parser.parse_args()
    cfg = load_config(args.config)
    main(args, cfg)
