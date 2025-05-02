import pandas as pd
import os

def sanity_check(df, name):
    print(f"\n🧪 Sanity Check: {name}")
    
    # 1. Check column presence
    required_cols = ["chrom", "position", "strand", "label"]
    for col in required_cols:
        if col not in df.columns:
            print(f"❌ Missing column: {col}")
            return

    # 2. Check chromosome formatting
    chr_prefix = df['chrom'].astype(str).str.startswith("chr").sum()
    print(f"✅ Chromosome format: {chr_prefix}/{len(df)} start with 'chr'")

    # 3. Label distribution
    label_counts = df['label'].value_counts().to_dict()
    print(f"✅ Label counts: {label_counts}")

    # 4. Strand values
    strand_vals = df['strand'].unique().tolist()
    if all(s in ["+", "-"] for s in strand_vals):
        print(f"✅ Valid strands: {strand_vals}")
    else:
        print(f"❌ Invalid strand values: {strand_vals}")

    # 5. NaN/Inf check
    nan_cols = df.isna().sum()[df.isna().sum() > 0]
    inf_cols = df.columns[df.apply(lambda s: s.isin([float("inf"), float("-inf")]).any())]
    if nan_cols.empty and inf_cols.empty:
        print("✅ No NaNs or Infs in features.")
    else:
        print(f"❌ NaNs in columns: {nan_cols.to_dict()}")
        print(f"❌ Infs in columns: {inf_cols.tolist()}")

    # 6. Feature value ranges
    feature_cols = [c for c in df.columns if c not in ["chrom", "position", "strand", "site_type", "label"]]
    stats = df[feature_cols].describe()
    print(f"📊 Feature summary:\n{stats.T[['min', 'max', 'mean', 'std']].round(2)}")

def main():
    input_dir = "data_train"
    for fname in os.listdir(input_dir):
        if fname.endswith(".csv"):
            path = os.path.join(input_dir, fname)
            df = pd.read_csv(path)
            sanity_check(df, name=fname)

if __name__ == "__main__":
    main()
