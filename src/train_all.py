import subprocess
import os
from datetime import datetime
import pandas as pd 

def main(tools, log_dir, train_dir, out_dir):
    # tools = ["stringtie", "isoquant", "universe"]
    tools.append("universe")
    site_types = ["tss", "tes"]
    models = ["xgboost", "randomforest"]
    # models = ["randomforest"]
    
    if not os.path.exists(train_dir):
        print(f"Missing training data directory: {train_dir}")
        return
    
    os.makedirs(log_dir, exist_ok=True)


    log_path = os.path.join(log_dir, "train_benchmark.log")

    with open(log_path, "w") as log:
        log.write(f" Training started: {datetime.now()}\n")

        for tool in tools:
            if tool == "universe":
                # concatenate all the tss and tes files from other tool types and drop duplicates
                # this is to create a universe dataset that contains all the tss and tes from all tools
            
                tss_files = [os.path.join(train_dir, f"{t}_tss_labeled.csv") for t in tools if t != "universe"]
                tes_files = [os.path.join(train_dir, f"{t}_tes_labeled.csv") for t in tools if t != "universe"]
                all_tss = pd.concat([pd.read_csv(f) for f in tss_files if os.path.exists(f) or print(f"file {f} not found")], ignore_index=True).reset_index(drop=True)
                all_tes = pd.concat([pd.read_csv(f) for f in tes_files if os.path.exists(f) or print(f"file {f} not found")], ignore_index=True).reset_index(drop=True)
                # drop duplicates
                all_tss = all_tss.drop_duplicates(subset=["chrom", "position", "strand"], keep="first")
                all_tes = all_tes.drop_duplicates(subset=["chrom", "position", "strand"], keep="first")
                # save to csv
                all_tss.to_csv(os.path.join(train_dir, "universe_tss_labeled.csv"), index=False)
                all_tes.to_csv(os.path.join(train_dir, "universe_tes_labeled.csv"), index=False)
                print("✅ Universe files created successfully.")
                

            for site in site_types:
                for model_type in models:
                    input_file = os.path.join(train_dir, f"{tool}_{site}_labeled.csv")
                    config_file = f"configs/{site}_config.yaml"
                    site_tag = f"{tool.upper()} - {site.upper()} - {model_type.upper()}"

                    if not os.path.exists(input_file):
                        print(f" Skipping missing input: {input_file}")
                        continue
                    
                    print(f"▶️  [START] {site_tag}")
                    cmd = [
                        "python", "src/train_model.py",
                        "--input", input_file,
                        "--config", config_file,
                        "--site_type", site,
                        "--model_type", model_type,
                        "--out_dir", out_dir,
                        "--val_chrom_file", os.path.join(train_dir, "validation_chromosomes.txt")
                    ]

                    print(f"🔁 Running: {site_tag}")
                    result = subprocess.run(cmd, capture_output=True, text=True)

                    log.write(f"\n===== {site_tag} =====\n")
                    log.write(result.stdout)
                    log.flush()
                    if result.returncode == 0:
                        print(f"✅ [DONE]  {site_tag}")
                    else:
                        print(f"❌ [FAIL]  {site_tag}")

                    if result.stderr.strip():
                        log.write(f"\n[stderr]\n{result.stderr}")


if __name__ == "__main__":
    main()