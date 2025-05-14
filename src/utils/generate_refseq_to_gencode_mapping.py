#!/usr/bin/env python3
import pandas as pd
import argparse
from ml_utils import normalize_chrom

p = argparse.ArgumentParser(
    description="Generate RefSeq to Gencode mapping file."
)
p.add_argument( '--ref2ucsc', type=str, required=True,
    help="Path to RefSeq to UCSC mapping file."
)
p.add_argument( '--gencode2ucsc', type=str, required=True,
    help="Path to Gencode to UCSC mapping file."
)
p.add_argument( '--output', type=str, required=True,
    help="Path to output file."
)
args = p.parse_args()
# 1) Load your mapping files.
#    Adjust sep, header, and column names if your files have headers.

ref2ucsc = pd.read_csv(
    args.ref2ucsc, sep="\t", header=None,
    names=["refseq", "ucsc"]
)
gencode2ucsc = pd.read_csv(
    args.gencode2ucsc, sep="\t", header=None,
    names=["gencode", "ucsc"]
)

# 2) Normalize the chromosome names.
ref2ucsc["ucsc"] = ref2ucsc["ucsc"].apply(normalize_chrom)
print(ref2ucsc["ucsc"])

# 2) Merge on the UCSC column.
#    This will keep only entries present in both maps.
merged = pd.merge(
    ref2ucsc,
    gencode2ucsc,
    on="ucsc",
    how="inner"
)

# 3) Extract the RefSeq→Gencode pairs.
ref2gencode = merged[["refseq", "gencode"]]

# 4) Save to a new file.
ref2gencode.to_csv(
    args.output,
    sep="\t",
    index=False,
    header=False
)

print(f"Mapped {len(ref2gencode)} RefSeq IDs to Gencode IDs.")
