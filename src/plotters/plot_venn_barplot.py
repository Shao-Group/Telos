import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from argparse import ArgumentParser

# Optional: reuse pretty names if available
try:
    # from generate_stage1_pr_curve import name_dict
    from calculate_true_false_stats import dataset_to_run_accession, pretty_model_names
except Exception:
    dataset_to_run_accession = {}


def load_summary(csv_path: str) -> pd.DataFrame:
    df = pd.read_csv(csv_path)
    # Normalize dataset names for display
    if dataset_to_run_accession:
        df["dataset_pretty"] = df["dataset"].apply(lambda x: dataset_to_run_accession.get(x, x))
    else:
        df["dataset_pretty"] = df["dataset"]
    return df


def plot_both_pct_barplot(baseline_dir: str,
                          predictions_root: str,
                          models: list,
                          output_path: str ,
                          site_type: str,
                          figsize=(12, 6)) -> str:
    baseline_csv = os.path.join(baseline_dir, f"venn_summary_{site_type}.csv")
    if not os.path.exists(baseline_csv):
        raise FileNotFoundError(f"Baseline summary not found: {baseline_csv}")
    df_base = load_summary(baseline_csv)[["dataset", "dataset_pretty", "both_pct"]].rename(
        columns={"both_pct": "baseline_both_pct"}
    )

    # Load models
    model_frames = []
    for model in models:
        model_csv = os.path.join(predictions_root, model, "venn_summary.csv")
        if not os.path.exists(model_csv):
            # Skip missing model silently to allow partial plots
            continue
        df_m = load_summary(model_csv)[["dataset", "dataset_pretty", "both_pct"]].rename(
            columns={"both_pct": f"{model}_both_pct"}
        )
        model_frames.append((model, df_m))

    if not model_frames:
        raise FileNotFoundError(f"No prediction summaries found under {predictions_root} for models: {models}")

    # Merge all on dataset
    df = df_base
    for model, df_m in model_frames:
        df = df.merge(df_m, on=["dataset", "dataset_pretty"], how="outer")

    # Group datasets: pacbio, cDNA, dRNA, SRR, then others
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

    df["group"] = df["dataset"].apply(classify_group)
    group_order = ["pacbio", "cdna", "drna", "srr", "other"]
    df["group_sort_key"] = df["group"].apply(lambda g: group_order.index(g) if g in group_order else len(group_order))
    df = df.sort_values(["group_sort_key", "dataset_pretty"]).reset_index(drop=True)

    # Prepare plot data
    categories = ["baseline_both_pct"] + [f"{m}_both_pct" for m, _ in model_frames]
    labels = [pretty_model_names.get(m) for m, _ in model_frames]
    labels = ["Baseline"] + labels
    print(labels)

    # Space datasets further apart by scaling x positions
    x_idx = np.arange(len(df))
    spacing_factor = 1.4  # >1 adds space between dataset groups
    x_pos = x_idx * spacing_factor
    width = 0.25 if len(categories) == 3 else min(0.8 / max(len(categories), 1), 0.3)

    fig, ax = plt.subplots(1, 1, figsize=figsize)

    for idx, col in enumerate(categories):
        offsets = [x_pos[i] + (idx - (len(categories) - 1) / 2) * width for i in range(len(df))]
        ax.bar(offsets, df[col].fillna(0.0) / 100.0, width=width, label=labels[idx])

    ax.set_xticks(list(x_pos))
    ax.set_xticklabels(df["dataset_pretty"], rotation=45, ha="right", fontsize=15)
    ax.set_ylabel("Jaccard Similarity", fontsize=16)
    # ax.set_title("Baseline vs Predictions (Jaccard Similarity)")
    ax.legend(fontsize=16, loc="upper right")
    ax.grid(axis="y", linestyle=":", alpha=0.4)

    # Add group titles centered above each contiguous group
    group_title_map = {
        "pacbio": "PacBio",
        "cdna": "cDNA",
        "drna": "dRNA",
        "srr": "Short Reads",
        "other": "Other",
    }

    # Ensure there is some headroom for titles
    max_frac = 0.0
    if len(df) > 0:
        max_frac = (df[categories].max(axis=1).fillna(0.0) / 100.0).max()
    upper = min(1.05, max_frac + 0.1)
    ax.set_ylim(0.0, upper)

    # Compute group ranges (start index inclusive, end index exclusive)
    ranges = []
    if len(df) > 0:
        start_idx = 0
        current_group = df.loc[0, "group"]
        for i in range(1, len(df)):
            g = df.loc[i, "group"]
            if g != current_group:
                ranges.append((current_group, start_idx, i))
                current_group = g
                start_idx = i
        ranges.append((current_group, start_idx, len(df)))

    # Place titles at the midpoint of each range
    y_text = ax.get_ylim()[1] - 0.06
    for g, start, end in ranges:
        if end <= start:
            continue
        mid = (x_pos[start] + x_pos[end - 1]) / 2.0
        title = group_title_map.get(g, g.title())
        if g == 'srr':
            mid = mid - 1.5*spacing_factor
        ax.text(mid, y_text, title, ha="center", va="center", fontsize=18, fontweight="bold")

    # Draw separators between groups
    group_changes = []
    if len(df) > 0:
        prev_group = df.loc[0, "group"]
        for i in range(1, len(df)):
            g = df.loc[i, "group"]
            if g != prev_group:
                # separator halfway between spaced positions
                group_changes.append((x_pos[i - 1] + x_pos[i]) / 2.0)
                prev_group = g
    for xpos in group_changes:
        ax.axvline(x=xpos, color="#888888", linestyle="--", linewidth=2)

    plt.tight_layout()

    # if output_path is None:
    #     output_dir = os.path.join("plots_individual", "venn_comparison")
    #     os.makedirs(output_dir, exist_ok=True)
    #     output_path = os.path.join(output_dir, "venn_both_pct_barplot.pdf")

    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return output_path


def main():
    parser = ArgumentParser()
    parser.add_argument("--site_type", required=True, help='Path to the configuration file')
    # parser.add_argument("--baseline_dir", default=os.path.join("plots_individual", "venn_baseline"),
    #                     help="Directory containing baseline venn_summary.csv")
    # parser.add_argument("--predictions_root", default=os.path.join("plots_individual", "venn_predictions"),
    #                     help="Root directory containing per-model folders with venn_summary.csv")
    # parser.add_argument("--models", nargs="*", default=["randomforest", "xgboost"],
    #                     help="Model folder names under predictions_root to include")
    # parser.add_argument("--output", default=None, help="Output file path for the barplot (PDF/PNG)")
    args = parser.parse_args()
    baseline_dir = os.path.join("plots_individual", "venn_baseline")
    predictions_root = os.path.join("plots_individual", "venn_predictions")
    models = ["randomforest", "xgboost"]
    output = os.path.join("plots_individual", "venn_comparison", f"venn_both_pct_barplot_{args.site_type}.pdf")

    out = plot_both_pct_barplot(
        baseline_dir=baseline_dir,
        predictions_root=predictions_root,
        models=models,
        output_path=output,
        site_type=args.site_type
    )
    print(f"Saved barplot to {out}")


if __name__ == "__main__":
    main()


