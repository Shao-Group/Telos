import subprocess
import os
from datetime import datetime

tools = ["stringtie", "isoquant"]
site_types = ["tss", "tes"]
models = ["xgboost", "randomforest"]
# models = ["randomforest"]
log_dir = "logs"
os.makedirs(log_dir, exist_ok=True)

log_path = os.path.join(log_dir, "train_benchmark.log")

with open(log_path, "w") as log:
    log.write(f" Training started: {datetime.now()}\n")

    for tool in tools:
        for site in site_types:
            for model_type in models:
                input_file = f"data_train/{tool}_{site}_labeled.csv"
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
                    "--model_type", model_type
                ]

                print(f"🔁 Running: {site_tag}")
                result = subprocess.run(cmd, capture_output=True, text=True)

                log.write(f"\n===== {site_tag} =====\n")
                log.write(result.stdout)
                if result.returncode == 0:
                    print(f"✅ [DONE]  {site_tag}")
                else:
                    print(f"❌ [FAIL]  {site_tag}")

                if result.stderr.strip():
                    log.write(f"\n[stderr]\n{result.stderr}")
