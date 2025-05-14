import os
import glob
import re
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

def plot_pr_curves(folder_path, out_folder, tool, auc_file):
    """
    Reads all ROC-like files in the specified folder, parses sensitivity (recall)
    and precision values, and plots Precision-Recall curves for each file.
    """
    plt.figure(figsize=(6, 4))

    # get auc map where key is the tool name and value is the auc in the two column csv
    auc_map = {}
    with open(auc_file, 'r') as f:
        for line in f:
            if line.startswith('tool'):
                continue
            tool_name, auc = line.strip().split(',')
            auc_map[tool_name] = float(auc)


    # Enable minor ticks and set their spacing
    # plt.minorticks_on()
    # plt.gca().xaxis.set_minor_locator(ticker.AutoMinorLocator(1))  # 4 minor intervals per major
    # plt.gca().yaxis.set_minor_locator(ticker.AutoMinorLocator(2))

    
    # Find all files in the folder
    file_paths = glob.glob(os.path.join(folder_path, '*'))
    # out_folder = f'out/{data_name}/plots/'

    for file_path in file_paths:
        recalls = []
        precisions = []
        
        with open(file_path, 'r') as f:
            for line in f:
                if line.startswith('ROC:'):
                    # Extract sensitivity and precision values
                    rec_match = re.search(r'sensitivity\s*=\s*([0-9.]+)', line)
                    prec_match = re.search(r'precision\s*=\s*([0-9.]+)', line)
                    if rec_match and prec_match:
                        recalls.append(float(rec_match.group(1)))
                        precisions.append(float(prec_match.group(1)))
        
        # Plot the curve for this file
        label = os.path.basename(file_path)
        label = label + " \n[ AuPR:" + f"{auc_map[label.split('.')[0].split('-updated-cov')[0]]/10000.0:.4f}" + " ]"
        # print(label)
        label = label.replace('-', ' ').replace('.roc', '').title()
        label = label.replace('Updated Cov', '' ) # '\n(Updated Coverage)')
        label = label.replace('Aupr', 'AuPR')
        if tool.lower() in label.lower() and 'universe' not in label.lower():
            plt.plot(recalls, precisions, label=label)

    # Draw major & minor grid
    # plt.grid(which='major', color='gray', linestyle='-', linewidth=0.5)
    # plt.grid(which='minor', color='gray', linestyle=':', linewidth=0.3)

    plt.xlabel('Recall (%)')
    plt.ylabel('Precision (%)')
    # plt.title('Precision-Recall Curves')
    plt.legend(loc='lower left', fontsize=11)
    plt.grid(False)
    # plt.grid(True, linewidth=0.1)   
    plt.tight_layout()
    # plt.subplots_adjust(right=0.6, top=0.8)
    plt.savefig(f'{out_folder}/pr_curves_compare_{tool}.pdf', format='pdf', dpi=300, bbox_inches='tight')
    # plt.show()

# Replace with the path to your folder containing the ROC files
plot_pr_curves('/datadisk1/ixk5174/long_reads_compare/out/gffcomp-results/pacbio_ENCFF694DIE/updated-cov/roc', 'out/pacbio_ENCFF694DIE/refSeq/plots', 'isoquant', 'out/pacbio_ENCFF694DIE/refSeq/predictions/transcripts/auc.csv')
# plot_pr_curves('/datadisk1/ixk5174/long_reads_compare/out/gffcomp-results/pacbio_ENCFF694DIE/updated-cov/roc', 'out/pacbio_ENCFF694DIE/plots', 'stringtie')