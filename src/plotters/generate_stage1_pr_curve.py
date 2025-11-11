import os, sys
sys.path.insert(
    0,
    os.path.abspath(
        os.path.join(os.path.dirname(__file__), "..")
    )
)

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from argparse import ArgumentParser
from config import load_config, Config
from sklearn.metrics import auc


# —— USER SETTINGS ——
data_names = [
    # "cDNA-NA12878","dRNA-ENCFF155CFF", "pacbio_ENCFF450VAU", "SRR307903"
        "cDNA-K562",
        "cDNA-ENCFF263YFG", 
        "cDNA-NA12878", 
        "dRNA-Hek293T", 
        "dRNA-ENCFF771DIX",
        # "dRNA-NA12878",
        "dRNA-ENCFF155CFF",
        "pacbio_ENCFF694DIE",
        "pacbio_ENCFF563QZR",
        "pacbio_ENCFF370NFS",
        "SRR307911",
        "SRR545695",
        "SRR315334",
        "SRR534307",
        "SRR545723",
        # "SRR307911",
        "SRR315323",
        "SRR534319",
        "SRR534291",
        "SRR387661"
]
name_dict = {
    "cDNA-K562" : "K562 cDNA",
    "cDNA-ENCFF263YFG" : "ENCFF263YFG cDNA",
    "cDNA-NA12878" : "NA12878 cDNA",
    "dRNA-Hek293T" : "Hek293T dRNA",
    "dRNA-ENCFF771DIX" : "ENCFF771DIX dRNA",
    "dRNA-ENCFF155CFF" : "ENCFF155CFF dRNA",
    "pacbio_ENCFF694DIE" : "ENCFF694DIE pacbio",
    "pacbio_ENCFF563QZR" : "ENCFF563QZR pacbio",
    "pacbio_ENCFF370NFS" : "ENCFF370NFS pacbio",
    # "pacbio_ENCFF212HLP" : "ENCFF212HLP pacbio",
    "SRR307911" : "SRR307911",
    "SRR545695" : "SRR545695",
    "SRR315334" : "SRR315334",
    "SRR534307" : "SRR534307",
    "SRR545723" : "SRR545723",
    "SRR307911" : "SRR307911",
    "SRR315323" : "SRR315323",
    "SRR534319" : "SRR534319",
    "SRR534291" : "SRR534291",
    "SRR387661" : "SRR387661"
}

# Define colors for different curve types
# Use a single base color per assembler and vary linestyles for baseline/models
assembler_base_colors = {
    "assembler1": "#1f77b4",  # blue
    "assembler2": "#ff7f0e"   # orange
}

# Color shade factors per model (same hue, different lightness)
model_shade_factors = {
    "xgboost": 0.75,       # darker
    "randomforest": 1.00,  # base
    "logistic": 1.15       # lighter
}
BASELINE_SHADE_FACTOR = 0.50  # darkest for baseline
BASELINE_LINESTYLE = "dashed"  # dotted

# Utility to lighten/darken a hex color while preserving hue and saturation
from matplotlib.colors import to_rgb, to_hex
import colorsys

def shade_hex_color(hex_color: str, lightness_factor: float) -> str:
    r, g, b = to_rgb(hex_color)
    h, l, s = colorsys.rgb_to_hls(r, g, b)
    # scale lightness around 0.5 for smoother changes
    new_l = max(0.0, min(1.0, l * lightness_factor))
    nr, ng, nb = colorsys.hls_to_rgb(h, new_l, s)
    return to_hex((nr, ng, nb))

STAGE1_PR_PLOT_FOLDER = "plots_individual/stage1_pr_curves/"

# ———————

def plot_pr(config_folder, is_train):
    site_types = ['tss', 'tes']
    tss_folder = os.path.join(STAGE1_PR_PLOT_FOLDER, "tss")
    tes_folder = os.path.join(STAGE1_PR_PLOT_FOLDER, "tes")
    os.makedirs(tss_folder, exist_ok=True)
    os.makedirs(tes_folder, exist_ok=True)
    suffix = "train" if is_train else "val"
    # Collect AUPR rows across all datasets/tools/sites
    aupr_records = []

    for i, name in enumerate(data_names):
        # if "ENCFF563QZR" not in name:
        #     continue
        tools =  ["stringtie", "scallop2"] if name.startswith("SRR") else ["stringtie", "isoquant"]
        tool_map = {
            "stringtie" : "assembler1",
            "scallop2" : "assembler2",
            "isoquant" : "assembler2"
        }

        configs  = {
            tools[0] : load_config(os.path.join(config_folder, f"{name}_{tools[0]}_config.pkl")),
            tools[1] : load_config(os.path.join(config_folder, f"{name}_{tools[1]}_config.pkl"))
        }
        for j, site in enumerate(site_types):
            fig, ax = plt.subplots(1, 1, figsize=(5, 5))
            print("Generating PR curve for {} {}".format(name, site))
            pr_files  = {
                tools[0] : os.listdir(configs[tools[0]].pr_data_dir), 
                tools[1] : os.listdir(configs[tools[1]].pr_data_dir)
            }
            
            # print(pr_files)
            end_precisions = []  # collect final precision values for ylim adjustment
            for tool, files in pr_files.items():
                for fname in files:
                    # if not fname.endswith(f"_{suffix}_pr_data.csv"):
                    #     continue
                    # assert fname.endswith(f"_{suffix}_pr_data.csv")
                    # print(tool + " " + fname)
                    # if "isoquant" not in tool:
                    #     continue
                    if site not in fname:
                        continue

                    # site_model = fname.replace(f"_{suffix}_pr_data.csv", "")
                    # if  "randomforest" in site_model :
                    #     continue
                    pr_df = pd.read_csv(os.path.join(configs[tool].pr_data_dir, fname))

                    model_type = fname.split("_")[1]  # e.g., "logistic"
                    label = f"{tool.title()} {model_type.title()}"
                    assembler_key = tool_map[tool]
                    base_color = assembler_base_colors.get(assembler_key, "#333333")
                    color = shade_hex_color(base_color, model_shade_factors.get(model_type, 1.0))

                    # print(f"Plotting {label} with color {color}")
                    # print(f"PR data: {pr_df.head(3)}")
                    
                    # full PR curve
                    line = ax.plot(
                        pr_df["recall"],
                        pr_df["precision"],
                        label=label,
                        color=color,
                        linewidth=2
                    )[0]

                    # Compute AUPR for this model curve using sklearn.metrics.auc (trapezoidal rule)
                    pr_sorted = pr_df.sort_values("recall")
                    aupr_val = float(auc(pr_sorted["recall"].values, pr_sorted["precision"].values))
                    aupr_records.append({
                        "dataset": name,
                        "tool": tool,
                        "site": site,
                        "model_type": model_type,
                        "aupr": aupr_val
                    })

                # baseline: coverage-based PR curve
                labeled_data_file = configs[tool].tss_labeled_file if site == "tss" else configs[tool].tes_labeled_file
                baseline_df = pd.read_csv(labeled_data_file, dtype={"chrom": str, "position": int})
                # print(f"Number of original points: {len(baseline_df)} for {tool} {site}")

                # print([c for c in baseline_df.columns])
                df_cov = pd.read_csv(configs[tool].cov_file, sep="\t", dtype={"tss_chrom": str, "tes_chrom": str})
                # print(f"Coverage file shape: {df_cov.shape}")

                df_cov = df_cov[[f"{site}_chrom", f"{site}_pos", "coverage"]]
                
                # Check for duplicates in coverage file
                n_unique_positions = df_cov[[f"{site}_chrom", f"{site}_pos"]].drop_duplicates().shape[0]
                # print(f"Unique positions in coverage file: {n_unique_positions}, Total rows: {len(df_cov)} for {tool} {site}")
                
                if len(df_cov) > n_unique_positions:
                    # print(f"WARNING: Coverage file has duplicates! Averaging coverage values...")
                    # Average coverage for duplicate positions
                    df_cov = df_cov.groupby([f"{site}_chrom", f"{site}_pos"])['coverage'].mean().reset_index()
                
                baseline_df = baseline_df.merge(df_cov, left_on=["chrom", "position"], right_on=[f"{site}_chrom", f"{site}_pos"], how="inner")
                # print(f"Number of filtered points: {len(baseline_df)} for {tool} {site}")
                # # print sample missing points
                # if baseline_df[baseline_df["coverage"].isna()].shape[0] > 0:
                #     print(baseline_df[baseline_df["coverage"].isna()].sample(min(10, baseline_df[baseline_df["coverage"].isna()].shape[0])))

                with open (configs[tool].validation_chromosomes_file, "r") as fp:
                    val_chrom = [line.rstrip("\n") for line in fp]

                # Filter to validation chromosomes and sort by coverage
                baseline_filtered = baseline_df.loc[baseline_df["chrom"].isin(val_chrom)].copy()
                baseline_filtered = baseline_filtered.sort_values('coverage', ascending=False)
                
                # Generate PR curve by gradually removing low coverage points
                precisions = []
                recalls = []
                n_total = len(baseline_filtered)
                n_positive = baseline_filtered['label'].sum()

                # print(f"Baseline columns: {baseline_filtered.columns}")
                # print(baseline_filtered[['chrom', 'position', f'{site}_chrom', f'{site}_pos', 'coverage', 'label']].head(50))
                # print(baseline_filtered[['chrom', 'position', f'{site}_chrom', f'{site}_pos', 'coverage', 'label']].tail(50))
                
                for k in range(1, n_total + 1):
                    # Keep top i highest coverage points
                    subset = baseline_filtered.iloc[:k]
                    n_kept = len(subset)
                    n_tp = subset['label'].sum()  # true positives in kept subset
                    n_fp = n_kept - n_tp  # false positives in kept subset
                    
                    # Precision = TP / (TP + FP)
                    precision = n_tp / n_kept 
                    # Recall = TP / (total positive)  
                    recall = n_tp / n_positive 
                    
                    precisions.append(precision)
                    recalls.append(recall)
                
                # Plot coverage-based baseline PR curve
                assembler_key = tool_map[tool]
                baseline_label = f"{assembler_key.title()} Baseline"
                baseline_color = shade_hex_color(assembler_base_colors.get(assembler_key, "#333333"), BASELINE_SHADE_FACTOR)
                
                baseline_line = ax.plot(
                    recalls,
                    precisions,
                    label=baseline_label,
                    linestyle=BASELINE_LINESTYLE,
                    linewidth=2,
                    color=baseline_color,
                    alpha=0.8
                )[0]

                # Compute AUPR for baseline curve
                if len(recalls) > 1 and len(precisions) == len(recalls):
                    # recalls are generated in increasing order
                    aupr_baseline = float(auc(recalls, precisions))
                else:
                    aupr_baseline = float("nan")
                aupr_records.append({
                    "dataset": name,
                    "tool": tool,
                    "site": site,
                    "model_type": "baseline",
                    "aupr": aupr_baseline
                })

                # add a point at the end of the curve
                ax.plot(
                    [recalls[-1]],
                    [precisions[-1]],
                    marker="o",
                    color=baseline_color,
                    markersize=5,
                    alpha=0.5
                )
                
                # record final precision of this tool's baseline for y-axis lower bound
                if len(precisions) > 0:
                    end_precisions.append(precisions[-1])
                
            # formatting per‐subplot
            # ax.set_title(f"{name_dict.get(name, name)} - {site.upper()}", fontsize=12)
            ax.set_xlabel("Recall", fontsize=14)
            ax.set_ylabel("Precision", fontsize=14)
            ax.tick_params(axis='both', labelsize=12)
            # Adjust y-axis lower bound to slightly below the minimum of final precisions
            if len(end_precisions) > 0:
                y_min = max(0.0, min(end_precisions) * 0.95)
                ax.set_ylim(y_min, 1.05)

            # ax.legend(loc='lower left', fontsize=9, frameon=True)

            plt.tight_layout()
            out_name = f"stage1_pr_{name}_{suffix}.pdf"
            out_path = os.path.join(tss_folder if site == "tss" else tes_folder, out_name)
            plt.savefig(out_path, dpi=300, format="pdf", bbox_inches="tight")
            plt.close(fig)
            print(f"Saved: {out_path}")

    # Save AUPR records to CSV inside plots_individual
    try:
        os.makedirs("plots_individual", exist_ok=True)
        aupr_df = pd.DataFrame(aupr_records, columns=["dataset", "tool", "site", "model_type", "aupr"])
        aupr_csv_path = os.path.join("plots_individual", f"stage1_aupr_values_{suffix}.csv")
        aupr_df.to_csv(aupr_csv_path, index=False)
        print(f"Saved AUPR table: {aupr_csv_path}")
    except Exception as e:
        print(f"Failed to save AUPR CSV: {e}")


def main():
    os.makedirs(STAGE1_PR_PLOT_FOLDER, exist_ok=True)
    parser = ArgumentParser()
    parser.add_argument('--config_folder', required=True, help='Path to the configuration file')
    parser.add_argument('--is_train', action='store_true', help='Is training chromosomes')
    args = parser.parse_args()
    plot_pr(args.config_folder, args.is_train)

if __name__ == "__main__":
    main()