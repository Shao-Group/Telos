#!/usr/bin/env python3
"""
Script to calculate true TSS/TES and false ones for each training dataset.
Generates both actual numbers and percentages in a tabular format.

This script follows the pattern from compile_metrics.py and uses the config files 
from project_config/ to locate the labeled data files.
"""

import os
import sys
import pandas as pd
from pathlib import Path
import argparse
import pickle
from typing import Dict, List, Tuple

# Add parent directory to path for imports
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from config import Config

pretty_model_names = {
    "baseline": "Baseline",
    "randomforest": "Telos-RF",
    "xgboost": "Telos-XGB"
}

dataset_to_run_accession = {
    "dRNA-ENCFF155CFF": "ENCFF155CFF",
    "dRNA-ENCFF771DIX": "ENCFF771DIX",
    "dRNA-NA12878": "dRNA NA12878",
    "dRNA-Hek293T": "ERR6053059",
    "cDNA-K562": "ERR6053079",
    "cDNA-NA12878": "cDNA NA12878",
    "cDNA-ENCFF263YFG": "ENCFF263YFG",
    "cDNA-ENCFF023EXJ": "ENCFF023EXJ",
    "pacbio_ENCFF450VAU": "ENCFF450VAU",
    "pacbio_ENCFF694DIE": "ENCFF694DIE",
    "pacbio_ENCFF563QZR": "ENCFF563QZR",
    "pacbio_ENCFF370NFS": "ENCFF370NFS",
    "SRR307903": "SRR307903",
    "SRR307911": "SRR307911",
    "SRR545695": "SRR545695",
    "SRR315334": "SRR315334",
    "SRR534307": "SRR534307",
    "SRR545723": "SRR545723",
    "SRR315323": "SRR315323",
    "SRR534319": "SRR534319",
    "SRR534291": "SRR534291",
    "SRR387661": "SRR387661"
}

def load_configs_from_directory(config_dir: str) -> List[Config]:
    """
    Load all training config files from the specified directory.
    
    Args:
        config_dir: Directory containing config pickle files
        
    Returns:
        List of training Config objects
    """
    train_configs = []
    
    config_files = [f for f in os.listdir(config_dir) if f.endswith('_config.pkl')]
    
    for config_file in config_files:
        try:
            config_path = os.path.join(config_dir, config_file)
            with open(config_path, 'rb') as f:
                config = pickle.load(f)
            
            # Only include training configs (exclude test configs)
            if 'train_output' in config.output_dir:
                train_configs.append(config)
            else:
                print(f"Skipping non-training config: {config_file}")
                
        except Exception as e:
            print(f"Error loading config {config_file}: {e}")
            continue
    
    return train_configs

def calculate_label_stats(labeled_file: str) -> Dict[str, int]:
    """
    Calculate true/false label statistics from a labeled TSV file.
    
    Args:
        labeled_file: Path to the labeled TSV file
        
    Returns:
        Dictionary with label statistics
    """
    try:
        df = pd.read_csv(labeled_file, dtype={"chrom": str})
        
        if 'label' not in df.columns:
            print(f"Warning: No 'label' column found in {labeled_file}")
            return {'total': 0, 'true_positive': 0, 'false_positive': 0}
        
        total_count = len(df)
        true_count = (df['label'] == 1).sum()
        false_count = (df['label'] == 0).sum()
        
        return {
            'total': total_count,
            'true_positive': true_count,
            'false_positive': false_count
        }
        
    except Exception as e:
        print(f"Error processing {labeled_file}: {e}")
        return {'total': 0, 'true_positive': 0, 'false_positive': 0}

def extract_dataset_info(config: Config) -> Dict[str, str]:
    """
    Extract dataset information from config object.
    
    Args:
        config: Config object
        
    Returns:
        Dictionary with dataset info
    """
    data_name = config.data_name
    
    # Split dataset and assembler from data_name
    if '_' in data_name:
        # Find the last underscore to split dataset from assembler
        last_underscore = data_name.rfind('_')
        dataset = data_name[:last_underscore]
        assembler = data_name[last_underscore + 1:]
    else:
        dataset = data_name
        assembler = 'unknown'
    
    return {
        'dataset': dataset_to_run_accession[dataset],
        'assembler': assembler,
        'data_name': data_name
    }

def compile_all_stats(configs: List[Config]) -> pd.DataFrame:
    """
    Compile label statistics for all training datasets.
    
    Args:
        configs: List of Config objects
        
    Returns:
        DataFrame with compiled statistics
    """
    all_data = []
    
    for config in configs:
        dataset_info = extract_dataset_info(config)
        
        # Process TSS data
        if os.path.exists(config.tss_labeled_file):
            tss_stats = calculate_label_stats(config.tss_labeled_file)
            tss_row = {
                **dataset_info,
                'site_type': 'TSS',
                'total_candidates': tss_stats['total'],
                'true_sites': tss_stats['true_positive'],
                'false_sites': tss_stats['false_positive'],
                'true_percentage': (tss_stats['true_positive'] / tss_stats['total'] * 100) if tss_stats['total'] > 0 else 0.0,
                'false_percentage': (tss_stats['false_positive'] / tss_stats['total'] * 100) if tss_stats['total'] > 0 else 0.0
            }
            all_data.append(tss_row)
            print(f"Processed TSS for {dataset_info['data_name']}: {tss_stats['true_positive']}/{tss_stats['total']} true sites")
        else:
            print(f"Warning: TSS labeled file not found for {dataset_info['data_name']}: {config.tss_labeled_file}")
        
        # Process TES data
        if os.path.exists(config.tes_labeled_file):
            tes_stats = calculate_label_stats(config.tes_labeled_file)
            tes_row = {
                **dataset_info,
                'site_type': 'TES',
                'total_candidates': tes_stats['total'],
                'true_sites': tes_stats['true_positive'],
                'false_sites': tes_stats['false_positive'],
                'true_percentage': (tes_stats['true_positive'] / tes_stats['total'] * 100) if tes_stats['total'] > 0 else 0.0,
                'false_percentage': (tes_stats['false_positive'] / tes_stats['total'] * 100) if tes_stats['total'] > 0 else 0.0
            }
            all_data.append(tes_row)
            print(f"Processed TES for {dataset_info['data_name']}: {tes_stats['true_positive']}/{tes_stats['total']} true sites")
        else:
            print(f"Warning: TES labeled file not found for {dataset_info['data_name']}: {config.tes_labeled_file}")
    
    if not all_data:
        print("No data found!")
        return pd.DataFrame()
    
    # Create DataFrame
    df = pd.DataFrame(all_data)
    
    # Reorder columns for better readability
    column_order = [
        'dataset', 'assembler', 'site_type', 'total_candidates', 
        'true_sites', 'false_sites', 'true_percentage', 'false_percentage', 'data_name'
    ]
    df = df[column_order]
    
    # Round percentages to 2 decimal places
    df['true_percentage'] = df['true_percentage'].round(2)
    df['false_percentage'] = df['false_percentage'].round(2)
    
    return df

def create_summary_table(df: pd.DataFrame) -> pd.DataFrame:
    """
    Create a summary table aggregated by dataset and site type.
    
    Args:
        df: DataFrame with detailed statistics
        
    Returns:
        DataFrame with summary statistics
    """
    summary_data = []
    
    # Group by dataset and site_type
    for (dataset, site_type), group in df.groupby(['dataset', 'site_type']):
        total_candidates = group['total_candidates'].sum()
        total_true = group['true_sites'].sum()
        total_false = group['false_sites'].sum()
        
        summary_row = {
            'dataset': dataset,
            'site_type': site_type,
            'assemblers': ', '.join(sorted(group['assembler'].unique())),
            'total_candidates': total_candidates,
            'true_sites': total_true,
            'false_sites': total_false,
            'true_percentage': (total_true / total_candidates * 100) if total_candidates > 0 else 0.0,
            'false_percentage': (total_false / total_candidates * 100) if total_candidates > 0 else 0.0
        }
        summary_data.append(summary_row)
    
    summary_df = pd.DataFrame(summary_data)
    
    # Round percentages
    summary_df['true_percentage'] = summary_df['true_percentage'].round(2)
    summary_df['false_percentage'] = summary_df['false_percentage'].round(2)
    
    return summary_df

def save_results(df: pd.DataFrame, summary_df: pd.DataFrame, output_prefix: str) -> None:
    """
    Save results to CSV files and create formatted tables including LaTeX.
    
    Args:
        df: Detailed DataFrame
        summary_df: Summary DataFrame
        output_prefix: Output file prefix
    """
    # Save detailed results
    detailed_file = f"{output_prefix}_detailed.csv"
    df.to_csv(detailed_file, index=False)
    print(f"Detailed results saved to: {detailed_file}")
    
    # Save summary results
    summary_file = f"{output_prefix}_summary.csv"
    summary_df.to_csv(summary_file, index=False)
    print(f"Summary results saved to: {summary_file}")
    
    # Create LaTeX table (compact two-assembler format)
    create_latex_table_compact(df, f"{output_prefix}_latex.txt")
    
    # Create formatted text tables
    detailed_table_file = f"{output_prefix}_detailed_table.txt"
    summary_table_file = f"{output_prefix}_summary_table.txt"
    
    with open(detailed_table_file, 'w') as f:
        f.write("="*120 + "\n")
        f.write("DETAILED TRUE/FALSE TSS/TES STATISTICS FOR TRAINING DATASETS\n")
        f.write("="*120 + "\n\n")
        f.write(df.to_string(index=False))
        f.write("\n\n")
        f.write("Legend:\n")
        f.write("- dataset: Original dataset name\n")
        f.write("- assembler: Assembly tool used (stringtie, isoquant, scallop2)\n")
        f.write("- site_type: TSS or TES\n")
        f.write("- total_candidates: Total number of candidate sites\n")
        f.write("- true_sites: Number of sites matching reference annotations (label=1)\n")
        f.write("- false_sites: Number of sites not matching reference annotations (label=0)\n")
        f.write("- true_percentage: Percentage of true sites\n")
        f.write("- false_percentage: Percentage of false sites\n")
    
    with open(summary_table_file, 'w') as f:
        f.write("="*100 + "\n")
        f.write("SUMMARY TRUE/FALSE TSS/TES STATISTICS BY DATASET\n")
        f.write("="*100 + "\n\n")
        f.write(summary_df.to_string(index=False))
        f.write("\n\n")
        f.write("Legend:\n")
        f.write("- dataset: Original dataset name\n")
        f.write("- site_type: TSS or TES\n")
        f.write("- assemblers: Assembly tools used for this dataset\n")
        f.write("- Statistics are aggregated across all assemblers for each dataset\n")
    
    print(f"Formatted detailed table saved to: {detailed_table_file}")
    print(f"Formatted summary table saved to: {summary_table_file}")

def create_latex_table(df: pd.DataFrame, output_file: str) -> None:
    """
    Create a LaTeX table in the preferred multirow format for true/false statistics.
    
    Args:
        df: DataFrame with true/false statistics
        output_file: Path to save the LaTeX file
    """
    # Sort data for proper grouping: by dataset, then assembler, then site_type (TSS first, then TES)
    df_sorted = df.copy()
    df_sorted['site_order'] = df_sorted['site_type'].map({'TSS': 0, 'TES': 1})
    df_sorted = df_sorted.sort_values(['dataset', 'assembler', 'site_order'])
    df_sorted = df_sorted.drop('site_order', axis=1)
    
    with open(output_file, 'w') as f:
        f.write("% LaTeX table snippet - True/False TSS/TES Statistics\n")
        f.write("% Use \\input{" + os.path.basename(output_file) + "} in your LaTeX document\n\n")
        
        # Create LaTeX table
        f.write("\\begin{table}[htbp]\n")
        f.write("\\centering\n")
        f.write("\\caption{True/False TSS/TES Statistics for Training Datasets}\n")
        f.write("\\label{tab:true_false_stats}\n")
        
        # Column alignment: left for dataset, left for assembler, left for site_type, center for numbers
        f.write("\\begin{tabular}{lllcccc}\n")
        f.write("\\toprule\n")
        
        # Header
        f.write("Dataset & Assembler & Site Type & Total & True Sites & False Sites & True \\% \\\\\n")
        f.write("\\midrule\n")
        
        # Generate data rows with multirow formatting
        current_dataset = None
        
        # First pass: count rows per dataset
        dataset_counts = df_sorted.groupby('dataset').size().to_dict()
        
        for idx, (_, row) in enumerate(df_sorted.iterrows()):
            dataset = row['dataset']
            assembler = row['assembler']
            site_type = row['site_type']
            
            # Format dataset name for display
            dataset_display = dataset.replace('_', '\\_')
            if dataset == 'pacbio_ENCFF450VAU':
                dataset_display = 'PacBio ENCFF450VAU'
            elif dataset == 'cDNA-NA12878':
                dataset_display = 'cDNA-NA12878'
            elif dataset == 'dRNA-NA12878':
                dataset_display = 'dRNA-NA12878'
            elif dataset == 'SRR307903':
                dataset_display = 'SRR307903'
            
            # Format assembler name
            assembler_display = assembler.replace('_', '\\_')
            if assembler == 'stringtie':
                assembler_display = 'StringTie'
            elif assembler == 'isoquant':
                assembler_display = 'IsoQuant'
            elif assembler == 'scallop2':
                assembler_display = 'Scallop2'
            
            row_parts = []
            
            # Handle dataset column with multirow
            if current_dataset != dataset:
                current_dataset = dataset
                dataset_row_count = dataset_counts[dataset]
                row_parts.append(f"\\multirow{{{dataset_row_count}}}{{*}}{{\\textbf{{{dataset_display}}}}}")
            else:
                row_parts.append("")
            
            # Add assembler and site type
            row_parts.extend([assembler_display, site_type])
            
            # Add statistics
            row_parts.append(f"{row['total_candidates']:,}")
            row_parts.append(f"{row['true_sites']:,}")
            row_parts.append(f"{row['false_sites']:,}")
            row_parts.append(f"{row['true_percentage']:.1f}\\%")
            
            f.write(" & ".join(row_parts) + " \\\\\n")
            
            # Add midrule after each dataset group (except the last one)
            if idx < len(df_sorted) - 1:
                next_dataset = df_sorted.iloc[idx + 1]['dataset']
                if next_dataset != current_dataset:
                    f.write("\\midrule\n")
        
        f.write("\\bottomrule\n")
        f.write("\\end{tabular}\n")
        f.write("\\end{table}\n")
    
    print(f"LaTeX table saved to: {output_file}")

def create_latex_table_compact(df: pd.DataFrame, output_file: str) -> None:
    """
    Create a LaTeX table matching the requested compact format:
    Dataset | Site Type | Assembler1 True | Assembler1 False | Assembler2 True | Assembler2 False
    
    Assumptions:
    - Assembler 1 = StringTie
    - Assembler 2 = IsoQuant
    - Percentages are per-assembler within each dataset+site type
    """
    # Keep only stringtie and isoquant
    # df2 = df[df['assembler'].isin(['stringtie', 'isoquant'])].copy()
    if df.empty:
        with open(output_file, 'w') as f:
            f.write('% No data for stringtie/isoquant to generate compact LaTeX table\n')
        print(f"LaTeX table saved to: {output_file}")
        return

    # Compute counts per (dataset, site_type, assembler)
    grp = df.groupby(['dataset', 'site_type', 'assembler'], as_index=False)[['true_sites', 'false_sites', 'total_candidates']].sum()

    # Pivot to wide with separate columns for each assembler
    wide_true = grp.pivot_table(index=['dataset', 'site_type'], columns='assembler', values='true_sites', fill_value=0)
    wide_false = grp.pivot_table(index=['dataset', 'site_type'], columns='assembler', values='false_sites', fill_value=0)
    wide_total = grp.pivot_table(index=['dataset', 'site_type'], columns='assembler', values='total_candidates', fill_value=0)

    # Helper to format "count (xx.x\%)"
    def fmt(count: int, total: int) -> str:
        pct = (count / total * 100.0) if total > 0 else 0.0
        return f"{count} ({pct:.1f}\\%)"

    # Build ordered rows (Dataset sorted, site type TSS then TES)
    idx = list(wide_true.index)
    idx_sorted = sorted(idx, key=lambda t: (t[0], 0 if t[1] == 'TSS' else 1))

    with open(output_file, 'w') as f:
        f.write("% LaTeX table snippet - Compact two-assembler format\n")
        f.write("% Use \\input{" + os.path.basename(output_file) + "} in your LaTeX document\n\n")
        f.write("\\begin{table}[htb]\n")
        f.write("\\centering\n")
        f.write("\\caption{Distribution of true and false candidate sites predicted by two assemblers in the four training datasets.}\n")
        f.write("\\label{tab:data-dist}\n\n")
        f.write("\\begin{tabular}{l l r r r r} \n")
        f.write("\\toprule % Top thick rule\n")
        f.write("\\textbf{Dataset} & \\textbf{\\begin{tabular}[c]{@{}l@{}}Site\\\\ Type\\end{tabular}} & \\multicolumn{2}{c}{\\textbf{Assembler 1}} & \\multicolumn{2}{c}{\\textbf{Assembler 2}} \\ \\\n")
        f.write("\\cmidrule(lr){3-4} \\cmidrule(lr){5-6} % Mid-rules for \"Assembler 1\" and \"Assembler 2\"\n")
        f.write("\\multicolumn{1}{c}{} & & \\textbf{True} & \\textbf{False} & \\textbf{True} & \\textbf{False} \\ \\\n")
        f.write("\\midrule % Medium rule to separate header from data\n")

        current_dataset = None
        for dataset, site_type in idx_sorted:
            true_s1 = int(wide_true.loc[(dataset, site_type)].get('stringtie', 0)) if ('stringtie' in wide_true.columns) else 0
            false_s1 = int(wide_false.loc[(dataset, site_type)].get('stringtie', 0)) if ('stringtie' in wide_false.columns) else 0
            tot_s1 = int(wide_total.loc[(dataset, site_type)].get('stringtie', 0)) if ('stringtie' in wide_total.columns) else 0

            # Assembler 2 depends on dataset: use scallop2 for SRR*, else isoquant
            assembler2_key = 'scallop2' if str(dataset).startswith('SRR') else 'isoquant'
            true_s2 = int(wide_true.loc[(dataset, site_type)].get(assembler2_key, 0)) if (assembler2_key in wide_true.columns) else 0
            false_s2 = int(wide_false.loc[(dataset, site_type)].get(assembler2_key, 0)) if (assembler2_key in wide_false.columns) else 0
            tot_s2 = int(wide_total.loc[(dataset, site_type)].get(assembler2_key, 0)) if (assembler2_key in wide_total.columns) else 0

            # Display-friendly dataset name
            dataset_display = dataset.replace('_', '\\_')

            # Row content
            dataset_cell = dataset_display if current_dataset != dataset else ""
            current_dataset = dataset

            s1_true_cell = fmt(true_s1, tot_s1)
            s1_false_cell = fmt(false_s1, tot_s1)
            s2_true_cell = fmt(true_s2, tot_s2)
            s2_false_cell = fmt(false_s2, tot_s2)

            f.write(f"{dataset_cell} & {site_type}  & {s1_true_cell} & {s1_false_cell} & {s2_true_cell} & {s2_false_cell}  \\\n")

            # Midrule between datasets (after TES or if next is a different dataset)
            # Determine next tuple
            next_idx = idx_sorted.index((dataset, site_type)) + 1
            if next_idx < len(idx_sorted):
                next_dataset, _ = idx_sorted[next_idx]
                if next_dataset != dataset:
                    f.write("\\midrule % Rule to separate datasets\n")

        f.write("\\bottomrule % Bottom thick rule\n")
        f.write("\\end{tabular}\n")
        f.write("\\end{table}\n")

    print(f"LaTeX table saved to: {output_file}")

def print_console_summary(df: pd.DataFrame, summary_df: pd.DataFrame) -> None:
    """
    Print a nicely formatted summary to console.
    
    Args:
        df: Detailed DataFrame
        summary_df: Summary DataFrame
    """
    print("\n" + "="*100)
    print("SUMMARY STATISTICS")
    print("="*100)
    print(summary_df.to_string(index=False))
    
    print("\n" + "="*120)
    print("DETAILED STATISTICS")
    print("="*120)
    print(df.to_string(index=False))
    
    # Print overall statistics
    print("\n" + "="*80)
    print("OVERALL STATISTICS")
    print("="*80)
    
    total_candidates = df['total_candidates'].sum()
    total_true = df['true_sites'].sum()
    total_false = df['false_sites'].sum()
    
    print(f"Total candidates across all datasets: {total_candidates:,}")
    print(f"Total true sites: {total_true:,} ({total_true/total_candidates*100:.2f}%)")
    print(f"Total false sites: {total_false:,} ({total_false/total_candidates*100:.2f}%)")
    
    # Statistics by site type
    print(f"\nBy site type:")
    for site_type in ['TSS', 'TES']:
        site_df = df[df['site_type'] == site_type]
        if not site_df.empty:
            site_total = site_df['total_candidates'].sum()
            site_true = site_df['true_sites'].sum()
            site_false = site_df['false_sites'].sum()
            print(f"  {site_type}: {site_total:,} total, {site_true:,} true ({site_true/site_total*100:.2f}%), {site_false:,} false ({site_false/site_total*100:.2f}%)")

def main():
    parser = argparse.ArgumentParser(
        description="Calculate true/false TSS/TES statistics for all training datasets"
    )
    parser.add_argument(
        '--output-prefix', '-o',
        default='true_false_tss_tes_stats',
        help='Output file prefix (default: true_false_tss_tes_stats)'
    )
    parser.add_argument(
        '--config-dir',
        default='project_config',
        help='Directory containing config pickle files (default: project_config)'
    )
    
    args = parser.parse_args()
    
    # Get the project root directory (assuming script is in src/plotters)
    script_dir = Path(__file__).parent
    project_root = script_dir.parent.parent
    
    # Convert relative paths to absolute paths
    if os.path.isabs(args.config_dir):
        config_dir = args.config_dir
    else:
        config_dir = str(project_root / args.config_dir)
    
    if not os.path.exists(config_dir):
        print(f"Error: Config directory {config_dir} does not exist!")
        return
    
    print(f"Loading configs from: {config_dir}")
    
    # Load training configs
    train_configs = load_configs_from_directory(config_dir)
    
    print(f"Found {len(train_configs)} training configs")
    
    if not train_configs:
        print("No valid training configs found!")
        return
    
    # Compile statistics
    print("\nCompiling label statistics...")
    df = compile_all_stats(train_configs)
    
    if df.empty:
        print("No statistics compiled!")
        return
    
    # Create summary
    summary_df = create_summary_table(df)
    
    # Determine output file paths
    if os.path.isabs(args.output_prefix):
        output_prefix = args.output_prefix
    else:
        output_prefix = str(project_root / args.output_prefix)
    
    # Save results
    save_results(df, summary_df, output_prefix)
    
    # Print console summary
    print_console_summary(df, summary_df)

if __name__ == "__main__":
    main()
