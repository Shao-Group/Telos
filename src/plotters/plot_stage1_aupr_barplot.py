import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from argparse import ArgumentParser

# Optional pretty names mapping reused from other scripts
try:
    from calculate_true_false_stats import dataset_to_run_accession, pretty_model_names
except Exception:
    dataset_to_run_accession = {}


def get_pretty_dataset_name(ds: str) -> str:
    if not ds:
        return ds
    return dataset_to_run_accession.get(ds, ds)


def classify_group(ds: str) -> str:
    s = (ds or "").lower()
    if s.startswith("pacbio") or "pacbio" in s:
        return "pacbio"
    if s.startswith("cdna") or "cdna" in s:
        return "cdna"
    if s.startswith("drna") or "drna" in s:
        return "drna"
    if s.startswith("srr"):
        return "srr"
    return "other"


def load_aupr_csv(csv_path: str) -> pd.DataFrame:
    if not os.path.exists(csv_path):
        raise FileNotFoundError(f"CSV not found: {csv_path}")
    df = pd.read_csv(csv_path)
    required_cols = {"dataset", "tool", "site", "model_type", "aupr"}
    missing = required_cols - set(df.columns)
    if missing:
        raise ValueError(f"CSV missing columns: {sorted(missing)}")
    return df


def prepare_barplot_dataframe(df: pd.DataFrame,
                              site: str,
                              tool_group: str | None = None) -> tuple[pd.DataFrame, list[str], list[str]]:
    data = df.copy()
    if site:
        data = data.loc[data["site"].str.lower() == site.lower()].copy()
    if tool_group:
        tg = tool_group.lower()
        if tg == "stringtie":
            data = data.loc[data["tool"].str.lower() == "stringtie"].copy()
        elif tg == "other":
            data = data.loc[data["tool"].str.lower() != "stringtie"].copy()
        else:
            raise ValueError(f"Unknown tool_group: {tool_group}")

    # Ensure uniqueness: exactly one row per (dataset, model_type)
    counts = (data.groupby(["dataset", "model_type"]).size().reset_index(name="n"))
    dupes = counts[counts["n"] != 1]
    if not dupes.empty:
        raise ValueError(
            f"Expected exactly one row per (dataset, model_type) after filtering; found duplicates: {dupes.to_dict(orient='records')}"
        )

    # Pivot to wide format without aggregation (unique keys)
    pivot = (data.pivot(index="dataset", columns="model_type", values="aupr")
                  .reset_index())
    pivot.columns.name = None

    # Categories and labels (keep only present columns)
    desired_order = ["baseline", "randomforest", "xgboost"]
    categories = [c for c in desired_order if c in pivot.columns]
    if not categories:
        raise ValueError("No model_type columns found to plot (expected baseline/randomforest/xgboost)")
    labels = [
        pretty_model_names.get(c)
        for c in categories
    ]


    # Pretty names and grouping
    pivot["dataset_pretty"] = pivot["dataset"].apply(get_pretty_dataset_name)
    pivot["group"] = pivot["dataset"].apply(classify_group)
    group_order = ["pacbio", "cdna", "drna", "srr"]
    pivot["group_sort_key"] = pivot["group"].apply(lambda g: group_order.index(g) if g in group_order else len(group_order))
    pivot = pivot.sort_values(["group_sort_key", "dataset_pretty"]).reset_index(drop=True)

    return pivot, categories, labels


def plot_stage1_aupr_barplot(df: pd.DataFrame,
                             site: str,
                             tool_group: str | None = None,
                             output_path: str | None = None,
                             figsize=(10, 6)) -> list:
    pivot, categories, labels = prepare_barplot_dataframe(df, site=site, tool_group=tool_group)
    
    # print(pivot.head())
    # Calculate improvements for this specific plot
    improvements = []
    
    for _, row in pivot.iterrows():
        baseline = row.get("baseline", 0)
        
        best_aupr = max([row.get(m, 0) for m in ["randomforest", "xgboost"] ])
        if best_aupr > 0:
            diff = best_aupr - baseline
            pct = (best_aupr / baseline - 1) * 100
            dataset_group = classify_group(row['dataset'])
            improvements.append({
                'site': site,
                'tool_group': tool_group,
                'dataset_group': dataset_group,
                'dataset': row['dataset'],
                'improvement_diff': diff,
                'improvement_pct': pct
            })
    
    avg_diff = np.mean([x['improvement_diff'] for x in improvements])
    avg_pct = np.mean([x['improvement_pct'] for x in improvements])
    print(f"{site}_{tool_group}: Avg improvement = {avg_diff:.4f} ({avg_pct:.1f}%) over {len(improvements)} datasets")

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
        max_val = (pivot[categories].max(axis=1)).max()
    upper = max(0.05, max_val + 0.28)
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
        "srr": "Short Reads"
    }

    y_text = ax.get_ylim()[1] - 0.06 * max(1.0, ax.get_ylim()[1])
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
        output_dir = os.path.join("plots_individual", "aupr_barplots")
        os.makedirs(output_dir, exist_ok=True)
        group_tag = tool_group if tool_group else "all"
        output_path = os.path.join(output_dir, f"stage1_aupr_barplot_{site}_{group_tag}.pdf")

    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return output_path, improvements


def main():
    parser = ArgumentParser()
    parser.add_argument("--csv", default=os.path.join("plots_individual", "stage1_aupr_values_val.csv"),
                        help="Path to AUPR CSV generated by generate_stage1_pr_curve.py")
    parser.add_argument("--site", choices=["tss", "tes"], help="Site to plot: tss or tes")
    parser.add_argument("--tool_group", choices=["stringtie", "other"], help="Tool group: stringtie or other")
    parser.add_argument("--output", default=None, help="Output file path for the barplot")
    parser.add_argument("--generate_all", action="store_true", help="Generate all four figures (tss/tes x stringtie/other)")
    args = parser.parse_args()

    df = load_aupr_csv(args.csv)

    all_improvements = []
    
    if args.generate_all:
        outputs = []
        for site in ["tss", "tes"]:
            for group in ["stringtie", "other"]:
                out, improvements = plot_stage1_aupr_barplot(
                    df=df,
                    site=site,
                    tool_group=group,
                    output_path=None
                )
                outputs.append(out)
                all_improvements.extend(improvements)
        print("Saved:")
        for p in outputs:
            print(p)
    else:
        if not args.site or not args.tool_group:
            raise SystemExit("Provide --site and --tool_group, or use --generate_all")
        out, improvements = plot_stage1_aupr_barplot(
            df=df,
            site=args.site,
            tool_group=args.tool_group,
            output_path=args.output
        )
        all_improvements.extend(improvements)
        print(f"Saved barplot to {out}")
    
    # Save improvements to CSV
    if all_improvements:
        improvements_df = pd.DataFrame(all_improvements)
        improvements_csv = os.path.join("plots_individual", "stage1_improvements_by_dataset_group.csv")
        improvements_df.to_csv(improvements_csv, index=False)
        print(f"Saved improvements to {improvements_csv}")
        
        # Print summary by dataset group
        summary = (improvements_df.groupby(['site', 'tool_group', 'dataset_group'])
                   .agg({'improvement_diff': ['mean', 'count'], 'improvement_pct': 'mean'})
                   .round(4))
        print("\n=== Summary by Dataset Group ===")
        print(f"{'Site':<4} {'Tool Group':<10} {'Dataset Group':<12} {'N':<3} {'Avg Diff':<10} {'Avg %':<8}")
        print("-" * 60)
        for (site, tool_group, dataset_group), row in summary.iterrows():
            n = int(row[('improvement_diff', 'count')])
            avg_diff = row[('improvement_diff', 'mean')]
            avg_pct = row[('improvement_pct', 'mean')]
            print(f"{site:<4} {tool_group:<10} {dataset_group:<12} {n:<3} {avg_diff:<10.4f} {avg_pct:<8.1f}%")


if __name__ == "__main__":
    main()


