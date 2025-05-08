import os
import re
import pandas as pd
import argparse
import numpy as np

def main(root_dir, output_dir, hardness="hard"):
    # Ensure the root directory exists
    if not os.path.exists(root_dir):
        raise FileNotFoundError(f"The directory {root_dir} does not exist.")

    # Regex for extracting metrics
    metric_pattern = re.compile(r'(\w*\s\w+ level):\s+([0-9.]+)\s+\|\s+([0-9.]+)')

    # Where we store all rows
    rows = []

    # Map directory to filter name
    filter_dirs = {
        'val-baseline': 'baseline',
        'filtered-randomforest': 'randomforest',
        'filtered-xgboost': 'xgboost'
    }

    for subdir, filter_name in filter_dirs.items():
        full_path = os.path.join(root_dir, subdir)
        for file in os.listdir(full_path):
            if file.endswith('.stats'):
                tool = 'isoquant' if 'isoquant' in file else 'stringtie'
                stats_path = os.path.join(full_path, file)
                with open(stats_path, 'r') as f:
                    content = f.read()

                for match in metric_pattern.finditer(content):
                    level = match.group(1).replace(' level', '').lower()
                    level.replace(' ', '_')
                    # print(level)
                    if (level.lower().strip() in ['base', 'exon', 'locus']):
                        # print(f"Skipping {level} level for {tool} in {filter_name}")
                        continue
                    sens = float(match.group(2))
                    prec = float(match.group(3))
                    
                    rows.append({
                        'tool': tool,
                        'filter': filter_name,
                        'level': level,
                        'sensitivity': sens,
                        'precision': prec
                    })

    # Create DataFrame
    df = pd.DataFrame(rows)

    # Pivot to wide format for better Markdown output
    summary = df.pivot_table(index=['tool', 'filter'], 
                             columns='level', 
                             values=['sensitivity', 'precision'])

    # Flatten MultiIndex columns
    summary.columns = [f"{stat}_{level}" for stat, level in summary.columns]
    summary.reset_index(inplace=True)

    # Save to Markdown
    with open(os.path.join(output_dir, f'summary-{hardness}.md'), 'w') as f:
        f.write(summary.to_markdown(index=False))

    # Optional: also save to CSV
    # summary.to_csv(os.path.join(output_dir, 'summary.tsv'), index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compile GFFCompare results.")
    parser.add_argument('-r',"--root_dir", dest = "root", type=str, required=True,  help="Root directory containing the results.")
    parser.add_argument('-o', "--output_dir", dest="output", type=str, default="out/", help="Directory to save the summary files.")
    parser.add_argument('-d', "--hardness", dest="hardness", type=str, default="hard", help="Hardness level of the results (default: hard).")
    args = parser.parse_args()
    
    main(args.root, args.output, args.hardness)
