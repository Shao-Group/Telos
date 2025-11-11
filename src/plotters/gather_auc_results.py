import os
import sys
import csv
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from argparse import ArgumentParser

# Ensure we can import config loader
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from config import load_config, Config  # noqa: E402
from calculate_true_false_stats import pretty_model_names


def read_auc_file(auc_path: str) -> pd.DataFrame:
    """Read an AUC file into a DataFrame with columns: label, auc.

    Expected format (example):
        label,auc\n
        xgboost,574.724\n
        randomforest,578.232\n
        baseline_val,584.445\n
    Also supports simple CSV without header or generic 'tool,auc'.
    Adds a best-effort normalized column (auc_normalized) if values look scaled.
    """
    if not os.path.exists(auc_path):
        return pd.DataFrame(columns=["label", "auc"])

    try:
        df = pd.read_csv(auc_path)
        # Normalize column names if needed
        cols = {c.lower(): c for c in df.columns}
        label_col = cols.get("label") or ("label" if "label" in df.columns else None)
        auc_col = cols.get("auc") or ("auc" if "auc" in df.columns else None)
        if label_col and auc_col:
            df = df.rename(columns={label_col: "label", auc_col: "auc"})[["label", "auc"]]
            # Clean labels and coerce numeric
            df["label"] = df["label"].astype(str).str.strip().replace({"baseline_val": "baseline"})
            df["auc"] = pd.to_numeric(df["auc"])
            df = df.dropna(subset=["auc"]).reset_index(drop=True)
            # Add normalized column if values seem scaled
            max_auc = df["auc"].max() if not df.empty else 0.0
            norm = 10000
            
            df["auc_normalized"] = (df["auc"] / norm) if norm else pd.NA
            return df
    except Exception:
        pass

    # Fallback: manual parse
    rows = []
    with open(auc_path, "r") as f:
        reader = csv.reader(f)
        for row in reader:
            if not row:
                continue
            if row[0].strip().lower() == "label":
                continue
            if len(row) >= 2:
                try:
                    rows.append({"label": row[0], "auc": float(row[1])})
                except ValueError:
                    continue
    df = pd.DataFrame(rows, columns=["label", "auc"])
    if not df.empty:
        df["label"] = df["label"].astype(str).str.strip().replace({"baseline_val": "baseline"})
        df["auc"] = pd.to_numeric(df["auc"], errors="coerce")
        df = df.dropna(subset=["auc"]).reset_index(drop=True)
        max_auc = df["auc"].max()
        norm = None
        if max_auc > 1.0:
            norm = 1000.0 if max_auc <= 1000.0 else 10000.0
        df["auc_normalized"] = (df["auc"] / norm) if norm else pd.NA
    return df


def get_pretty_dataset_name(dataset_name: str) -> str:
    """Convert dataset name to pretty name using the mapping."""
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
    
    return dataset_to_run_accession.get(dataset_name, dataset_name)


def get_pretty_assembler_name(assembler: str) -> str:
    """Convert assembler name to pretty name."""
    assembler_mapping = {
        "isoquant": "IsoQuant",
        "stringtie": "StringTie",
        "scallop2": "Scallop2"
    }
    return assembler_mapping.get(assembler, assembler.title())


def format_auc_value(auc_value: float) -> str:
    """Format AUC value for LaTeX table."""
    if pd.isna(auc_value):
        return "---"
    # Convert to percentag
    # e and round to 2 decimal places
    return f"{auc_value:.4f}"


def generate_latex_table(df: pd.DataFrame, use_normalized: bool = True) -> str:
    """Generate LaTeX table code using booktabs package."""
    
    # Filter to only include datasets from data_names
    data_names = [
        "cDNA-K562",
        "cDNA-ENCFF263YFG", 
        "cDNA-NA12878", 
        "dRNA-Hek293T", 
        "dRNA-ENCFF771DIX",
        "dRNA-ENCFF155CFF",
        "pacbio_ENCFF694DIE",
        "pacbio_ENCFF563QZR",
        "pacbio_ENCFF370NFS",
        "SRR307911",
        "SRR545695",
        "SRR315334",
        "SRR534307",
        "SRR545723",
        "SRR315323",
        "SRR534319",
        "SRR534291",
        "SRR387661"
    ]
    
    # Filter dataframe to only include datasets in data_names
    df_filtered = df[df['dataset'].isin(data_names)].copy()
    
    # Choose which AUC column to use
    auc_column = 'auc_normalized' if use_normalized else 'auc'
    
    # Pivot the data to have columns for baseline, randomforest, xgboost
    pivot_df = df_filtered.pivot_table(
        index=['dataset', 'tool'], 
        columns='label', 
        values=auc_column,
        aggfunc='first'
    ).reset_index()
    
    # Rename columns for LaTeX
    pivot_df.columns.name = None
    
    # Ensure we have the expected columns
    expected_columns = ['baseline', 'randomforest', 'xgboost']
    for col in expected_columns:
        if col not in pivot_df.columns:
            pivot_df[col] = pd.NA
    
    # Sort by dataset and tool
    pivot_df = pivot_df.sort_values(['dataset', 'tool']).reset_index(drop=True)
    
    # Generate LaTeX code
    latex_lines = []
    latex_lines.append("\\begin{table}[htbp]")
    latex_lines.append("\\centering")
    latex_lines.append("\\caption{AUC Results for Transcript Prediction}")
    latex_lines.append("\\label{tab:auc_results}")
    latex_lines.append("\\begin{tabular}{l l c c c}")
    latex_lines.append("\\toprule")
    latex_lines.append("Dataset & Assembler & Baseline & Random Forest & XGBoost \\\\")
    latex_lines.append("\\midrule")
    
    current_dataset = None
    dataset_count = 0
    total_datasets = len(pivot_df['dataset'].unique())
    
    for _, row in pivot_df.iterrows():
        dataset = row['dataset']
        tool = row['tool']
        
        # Add dataset name only once per dataset group
        if dataset != current_dataset:
            # Add midrule after previous dataset (except for first dataset)
            if current_dataset is not None:
                latex_lines.append("\\midrule")
            
            pretty_dataset = get_pretty_dataset_name(dataset)
            pretty_tool = get_pretty_assembler_name(tool)
            latex_lines.append(f"{pretty_dataset} & {pretty_tool} & "
                             f"{format_auc_value(row['baseline'])} & "
                             f"{format_auc_value(row['randomforest'])} & "
                             f"{format_auc_value(row['xgboost'])} \\\\")
            current_dataset = dataset
            dataset_count += 1
        else:
            # Same dataset, different assembler
            pretty_tool = get_pretty_assembler_name(tool)
            latex_lines.append(f"& {pretty_tool} & "
                             f"{format_auc_value(row['baseline'])} & "
                             f"{format_auc_value(row['randomforest'])} & "
                             f"{format_auc_value(row['xgboost'])} \\\\")
    
    latex_lines.append("\\bottomrule")
    latex_lines.append("\\end{tabular}")
    latex_lines.append("\\end{table}")
    
    return "\n".join(latex_lines)


def classify_group(dataset: str) -> str:
    s = (dataset or "").lower()
    if s.startswith("pacbio") or "pacbio" in s:
        return "pacbio"
    if s.startswith("cdna") or "cdna" in s:
        return "cdna"
    if s.startswith("drna") or "drna" in s:
        return "drna"
    if s.startswith("srr"):
        return "srr"
    return "other"


def plot_auc_barplot(df: pd.DataFrame,
                     use_normalized: bool = True,
                     tool_filter: str | None = None,
                     output_path: str | None = None,
                     figsize=(10, 6)) -> str:
    """Plot grouped barplots of AUC values per dataset comparing labels.

    df: long-form dataframe with columns [dataset, tool, label, auc, auc_normalized]
    use_normalized: when True, use 'auc_normalized'; otherwise 'auc'
    agg: aggregation across tools per dataset ('mean' or 'max')
    tool_filter: if provided, restrict to a single assembler/tool
    output_path: path to save the figure; defaults under plots_individual
    """
    auc_col = "auc_normalized" if use_normalized else "auc"

    data = df.copy()

    data = (data.groupby(["dataset", "label"], as_index=False)[auc_col]
                    .max())

    # Pivot to columns per label (baseline, randomforest, xgboost)
    pivot = data.pivot_table(index=["dataset"], columns="label", values=auc_col, aggfunc="first").reset_index()
    pivot.columns.name = None

    # Prepare categories in order
    categories = [c for c in ["baseline", "randomforest", "xgboost"] if c in pivot.columns]
    labels = [
        pretty_model_names.get(c)
        for c in categories
    ]
    if not categories:
        raise ValueError("No AUC labels found to plot (expected baseline/randomforest/xgboost)")

    # Pretty names and grouping
    pivot["dataset_pretty"] = pivot["dataset"].apply(get_pretty_dataset_name)
    pivot["group"] = pivot["dataset"].apply(classify_group)
    group_order = ["pacbio", "cdna", "drna", "srr", "other"]
    pivot["group_sort_key"] = pivot["group"].apply(lambda g: group_order.index(g) if g in group_order else len(group_order))
    pivot = pivot.sort_values(["group_sort_key", "dataset_pretty"]).reset_index(drop=True)

    # X positions with spacing between datasets
    x_idx = np.arange(len(pivot))
    spacing_factor = 1.4
    x_pos = x_idx * spacing_factor
    width = 0.25 if len(categories) == 3 else min(0.8 / max(len(categories), 1), 0.3)

    fig, ax = plt.subplots(1, 1, figsize=figsize)

    for idx, col in enumerate(categories):
        offsets = [x_pos[i] + (idx - (len(categories) - 1) / 2) * width for i in range(len(pivot))]
        ax.bar(offsets, pivot[col].fillna(0.0), width=width, label=labels[idx])

    ax.set_xticks(list(x_pos))
    ax.set_xticklabels(pivot["dataset_pretty"], rotation=45, ha="right", fontsize=13)
    ax.set_ylabel("Area Under Precision-Recall Curve (AuPR)", fontsize=13)
    ax.legend(fontsize=16)
    ax.grid(axis="y", linestyle=":", alpha=0.4)

    # Headroom
    max_val = 0.0
    if len(pivot) > 0:
        max_val = (pivot[categories].max(axis=1).fillna(0.0)).max()
    upper = max_val + (0.03 * (1.0 if use_normalized else max(1.0, max_val)))
    ax.set_ylim(0.0, upper)

    # Group titles
    ranges = []
    if len(pivot) > 0:
        start_idx = 0
        current_group = pivot.loc[0, "group"]
        for i in range(1, len(pivot)):
            g = pivot.loc[i, "group"]
            if g != current_group:
                ranges.append((current_group, start_idx, i))
                current_group = g
                start_idx = i
        ranges.append((current_group, start_idx, len(pivot)))

    group_title_map = {
        "pacbio": "PacBio",
        "cdna": "cDNA",
        "drna": "dRNA",
        "srr": "Short Reads",
        "other": "Other",
    }

    y_text = ax.get_ylim()[1] - 0.06 * (ax.get_ylim()[1] if use_normalized else max(1.0, ax.get_ylim()[1]))
    for g, start, end in ranges:
        if end <= start:
            continue
        mid = (x_pos[start] + x_pos[end - 1]) / 2.0
        title = group_title_map.get(g, g.title())
        if g == "srr":
            mid = mid - 2.2 * spacing_factor
        ax.text(mid, y_text, title, ha="center", va="center", fontsize=18, fontweight="bold")

    # Group separators
    group_changes = []
    if len(pivot) > 0:
        prev_group = pivot.loc[0, "group"]
        for i in range(1, len(pivot)):
            g = pivot.loc[i, "group"]
            if g != prev_group:
                group_changes.append((x_pos[i - 1] + x_pos[i]) / 2.0)
                prev_group = g
    for xpos in group_changes:
        ax.axvline(x=xpos, color="#888888", linestyle="--", linewidth=2)

    plt.tight_layout()

    if output_path is None:
        output_dir = os.path.join("plots_individual", "auc_barplots")
        os.makedirs(output_dir, exist_ok=True)
        output_path = os.path.join(output_dir, "auc_barplot.pdf")

    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return output_path


def _filter_df_for_tool_mode(df: pd.DataFrame, tool_mode: str) -> pd.DataFrame:
    """Filter rows for a specific assembler mode.

    tool_mode: 'stringtie' or 'other' (scallop2 for SRR, isoquant otherwise)
    """
    data = df.copy()
    if tool_mode == "stringtie":
        return data[data["tool"].str.lower() == "stringtie"]
    if tool_mode == "other":
        # Keep rows where tool matches dataset-specific other assembler
        def keep_row(row) -> bool:
            dataset = str(row.get("dataset", ""))
            tool = str(row.get("tool", "")).lower()
            group = classify_group(dataset)
            expected = "scallop2" if group == "srr" else "isoquant"
            return tool == expected
        mask = data.apply(keep_row, axis=1)
        return data[mask]
    return data


def plot_auc_barplots_two_tools(df: pd.DataFrame,
                                use_normalized: bool = True,
                                agg: str = "mean",
                                output_dir: str | None = None,
                                split: str | None = None) -> tuple[str, str]:
    """Generate two barplots: one for StringTie, one for the dataset-appropriate other assembler.

    Returns: (path_stringtie, path_other)
    """
    data_str = _filter_df_for_tool_mode(df, "stringtie")
    data_oth = _filter_df_for_tool_mode(df, "other")

    if output_dir is None:
        output_dir = os.path.join("plots_individual", "auc_barplots")
    os.makedirs(output_dir, exist_ok=True)

    suffix = f"_{split}" if split else ""
    out_str = os.path.join(output_dir, f"auc_barplot_stringtie{suffix}.pdf")
    out_oth = os.path.join(output_dir, f"auc_barplot_other{suffix}.pdf")

    path1 = plot_auc_barplot(
        data_str,
        use_normalized=use_normalized,
        # agg=agg,
        tool_filter="stringtie",
        output_path=out_str,
    )
    path2 = plot_auc_barplot(
        data_oth,
        use_normalized=use_normalized,
        # agg=agg,
        tool_filter=None,  # already filtered per-dataset to correct assembler
        output_path=out_oth,
    )
    return path1, path2


def main():
    parser = ArgumentParser()
    parser.add_argument("--config_folder", default="project_config", help="Folder with *_config.pkl files")
    parser.add_argument("--split", choices=["train", "val"], default="val", help="Which AUC file to gather")
    parser.add_argument("--output", default=None, help="Output CSV path under plots_individual/")
    parser.add_argument("--latex", action="store_true", help="Generate LaTeX table output")
    parser.add_argument("--latex_output", default=None, help="Output LaTeX file path")
    parser.add_argument("--use_raw_auc", action="store_true", help="Use raw AUC values instead of normalized")
    parser.add_argument("--barplot", action="store_true", help="Generate grouped AUC barplot")
    parser.add_argument("--barplot_output", default=None, help="Output path for AUC barplot (PDF/PNG)")
    parser.add_argument("--barplot_agg", choices=["mean", "max"], default="mean", help="Aggregation across assemblers")
    parser.add_argument("--barplot_tool", default=None, help="Filter to a specific assembler (overrides aggregation)")
    parser.add_argument("--barplot_two_tools", action="store_true", help="Generate two barplots: StringTie and other assembler per dataset")
    args = parser.parse_args()

    # Resolve paths relative to project root
    script_dir = os.path.dirname(__file__)
    project_root = os.path.abspath(os.path.join(script_dir, "..", ".."))

    config_dir = args.config_folder
    if not os.path.isabs(config_dir):
        config_dir = os.path.join(project_root, config_dir)

    if not os.path.isdir(config_dir):
        raise FileNotFoundError(f"Config folder not found: {config_dir}")

    out_dir = os.path.join(project_root, "plots_individual")
    os.makedirs(out_dir, exist_ok=True)
    out_path = args.output or os.path.join(out_dir, f"all_auc_{args.split}.csv")

    records = []
    cfg_files = sorted([f for f in os.listdir(config_dir) if f.endswith("_config.pkl")])
    for cfg_name in cfg_files:
        cfg_path = os.path.join(config_dir, cfg_name)
        try:
            cfg: Config = load_config(cfg_path)
            if cfg.is_training:
                continue
        except Exception as e:
            print(f"Skipping {cfg_name}: {e}")
            continue

        auc_path = cfg.auc_file_train if args.split == "train" else cfg.auc_file_val
        df_auc = read_auc_file(auc_path)
        if df_auc.empty:
            print(f"No AUC data found for {cfg.data_name} at {auc_path}")
            continue

        # Derive dataset and tool from data_name (split at last underscore)
        data_name = getattr(cfg, "data_name", cfg_name.replace("_config.pkl", ""))
        if "_" in data_name:
            base, tool = data_name.rsplit("_", 1)
        else:
            base, tool = data_name, "unknown"

        for _, row in df_auc.iterrows():
            records.append({
                "dataset": base,
                "tool": tool,
                "label": row["label"],
                "auc": row["auc"],
                "auc_normalized": row.get("auc_normalized", pd.NA),
                "split": args.split,
                "auc_file": auc_path,
            })

    if not records:
        print("No AUC records collected.")
        # Still emit an empty file with headers for convenience
        pd.DataFrame(columns=["dataset", "tool", "label", "auc", "auc_normalized", "split", "auc_file"]).to_csv(out_path, index=False)
        print(f"Wrote empty AUC table to {out_path}")
        return

    out_df = pd.DataFrame(records)
    # Sort for readability
    out_df = out_df.sort_values(["dataset", "tool", "label"]).reset_index(drop=True)
    out_df.to_csv(out_path, index=False)
    print(f"Saved aggregated AUC table to {out_path} with {len(out_df)} rows")
    
    # Generate LaTeX table if requested
    if args.latex:
        latex_output_path = args.latex_output or os.path.join(out_dir, f"auc_results_table_{args.split}.tex")
        use_normalized = not args.use_raw_auc
        
        latex_code = generate_latex_table(out_df, use_normalized=use_normalized)
        
        # Write LaTeX to file
        with open(latex_output_path, 'w') as f:
            f.write(latex_code)
        
        print(f"Generated LaTeX table")
        print(f"Saved to: {latex_output_path}")
        
        # Also print the LaTeX code to stdout for easy copying
        print("\n" + "="*50)
        print("LaTeX CODE:")
        print("="*50)
        print(latex_code)

    # Generate AUC barplot if requested
    if args.barplot:
        use_normalized = not args.use_raw_auc
        barplot_out = args.barplot_output or os.path.join(out_dir, f"auc_barplot_{args.split}.pdf")
        out_path_plot = plot_auc_barplot(
            out_df,
            use_normalized=use_normalized,
            agg=args.barplot_agg,
            tool_filter=args.barplot_tool,
            output_path=barplot_out,
        )
        print(f"Saved AUPR barplot to {out_path_plot}")

    # Generate two AUC barplots (per assembler) if requested
    if args.barplot_two_tools:
        use_normalized = not args.use_raw_auc
        path_str, path_oth = plot_auc_barplots_two_tools(
            out_df,
            use_normalized=use_normalized,
            agg=args.barplot_agg,
            output_dir=os.path.join(out_dir, "auc_barplots"),
            split=args.split,
        )
        print(f"Saved AUPR barplots:\n  StringTie -> {path_str}\n  Other     -> {path_oth}")


if __name__ == "__main__":
    main()


