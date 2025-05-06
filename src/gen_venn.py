import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2

def load_site_ids(filepath):
    df = pd.read_csv(filepath)
    return set(df['chrom'].astype(str) + ":" + df['position'].astype(str) + ":" + df['strand'].astype(str))

# find intersection of sites if they are withing 100bp
def find_intersection(site1, site2):
    # Split the site strings into components
    chrom1, pos1, strand1 = site1.split(":")
    chrom2, pos2, strand2 = site2.split(":")

    # Check if they are on the same chromosome, same strand, and within 50bp
    if chrom1 == chrom2 and strand1 == strand2 and abs(int(pos1) - int(pos2)) <= 50:
        return True
    return False

# find intersection of two sets based on the above function
def find_intersection_sets(set1, set2):
    intersection = set()
    for site1 in set1:
        for site2 in set2:
            if find_intersection(site1, site2):
                intersection.add(site1)
                break
    return intersection

paths = {
    "TSS": {
        "IsoQuant": "data_train/isoquant_tss_labeled.csv",
        "StringTie": "data_train/stringtie_tss_labeled.csv"
    },
    "TES": {
        "IsoQuant": "data_train/isoquant_tes_labeled.csv",
        "StringTie": "data_train/stringtie_tes_labeled.csv"
    }
}

# Create subplot figure
fig, axes = plt.subplots(1, 2, figsize=(12, 6))
for i, (site_type, tool_paths) in enumerate(paths.items()):
    isoquant_sites = load_site_ids(tool_paths["IsoQuant"])
    stringtie_sites = load_site_ids(tool_paths["StringTie"])

    intersection = find_intersection_sets(isoquant_sites, stringtie_sites)
    only_isoquant = len(isoquant_sites - intersection)
    only_stringtie = len(stringtie_sites - intersection)
    both = len(intersection)

    venn = venn2(subsets=(only_isoquant, only_stringtie, both), 
          set_labels=('IsoQuant', 'StringTie'), 
          ax=axes[i])
    # Increase font size
    for text in venn.set_labels:
        if text:
            text.set_fontsize(16)
    for text in venn.subset_labels:
        if text:
            text.set_fontsize(16)
    axes[i].set_title(f'{site_type} Site Overlap')

# Final adjustments
plt.tight_layout()
plt.savefig("out/plots/merged_venn_diagram.png")
plt.show()

# print("Venn diagrams saved as out/plots/venn_tss.png and out/plots/venn_tes.png.")
