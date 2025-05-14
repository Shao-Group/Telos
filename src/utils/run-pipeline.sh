#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

# ─── CONFIG ─────────────────────────────────────────────────────────────────────
home="/datadisk1/ixk5174/tools/tss-tes-project"
LOG_DIR="$home/logs"
cd $home

# Load conda
eval "$(conda shell.bash hook)"
# Load conda environment
conda activate irtesam-berth

# ─── MAIN LOOP ───────────────────────────────────────────────────────────────────
echo "✅ Running pipeline for Long read data with gencode..."


#  ─── SKIP PREVIOUS RUNS ────────────────────────────────────────────────────────────────
echo "✨ Running pipeline for direct RNA data..."

mkdir -p "$LOG_DIR/nanopore_dRNA_NA12878/gencode"
python src/pipeline.py --method long \
        --bam_file data/nanopore_dRNA_NA12878/NA12878-DirectRNA.sorted.bam \
            --data_name nanopore_dRNA_NA12878 \
                --reference gencode > "$LOG_DIR/nanopore_dRNA_NA12878/gencode/run.log" 2>&1


mkdir -p "$LOG_DIR/pacbio_ENCFF450VAU/gencode"
echo "✨ Running pipeline for cDNA data..."
python src/pipeline.py --method long \
        --bam_file nanopore_cDNA_NA12878/NA12878-cDNA.sorted.bam \
            --data_name nanopore_cDNA_NA12878 \
                --reference gencode > "$LOG_DIR/nanopore_cDNA_NA12878/gencode/run.log" 2>&1

mkdir -p "$LOG_DIR/pacbio_ENCFF450VAU/gencode"
echo "✨ Running pipeline for PacBio data ENCFF450VAU..."
python src/pipeline.py --method long \
        --bam_file pacbio_ENCFF450VAU/ENCFF450VAU.sorted.bam \
            --data_name pacbio_ENCFF450VAU \
                --reference gencode > "$LOG_DIR/pacbio_ENCFF450VAU/gencode/run.log" 2>&1

mkdir -p "$LOG_DIR/pacbio_ENCFF694DIE/gencode"
echo "✨ Running pipeline for PacBio data ENCFF694DIE..."
python src/pipeline.py --method long \
        --bam_file pacbio_ENCFF694DIE/ENCFF694DIE.sorted.bam \
            --data_name pacbio_ENCFF694DIE \
                --reference gencode > "$LOG_DIR/pacbio_ENCFF694DIE/gencode/run.log" 2>&1    


mkdir -p "$LOG_DIR/SRR307903_hisat/ensembl"
echo "✨ Running pipeline for Short read data SRR307903..."
python src/pipeline.py --method short \
        --bam_file SRR307903_hisat/hisat.sorted.bam \
            --data_name SRR307903_hisat \
                --reference ensembl > "$LOG_DIR/SRR307903_hisat/ensembl/run.log" 2>&1

mkdir -p "$LOG_DIR/SRR307911_hisat/ensembl"
echo "✨ Running pipeline for Short read data SRR307911..."
python src/pipeline.py --method short \
        --bam_file SRR307911_hisat/hisat.sorted.bam \
            --data_name SRR307911_hisat \
                --reference ensembl > "$LOG_DIR/SRR307911_hisat/ensembl/run.log" 2>&1


# --------- SKIP PREVIOUS RUNS ---------

mkdir -p "$LOG_DIR/nanopore_dRNA_NA12878/refSeq"
echo "✅ Running pipeline for Long read data with refseq..."
python src/pipeline.py --method long \
        --bam_file data/nanopore_dRNA_NA12878/NA12878-DirectRNA.sorted.bam \
            --data_name nanopore_dRNA_NA12878 \
                --reference refSeq > "$LOG_DIR/nanopore_dRNA_NA12878/refSeq/run.log" 2>&1

mkdir -p "$LOG_DIR/nanopore_cDNA_NA12878/refSeq"
echo "✨ Running pipeline for cDNA data..."
python src/pipeline.py --method long \
        --bam_file nanopore_cDNA_NA12878/NA12878-cDNA.sorted.bam \
            --data_name nanopore_cDNA_NA12878 \
                --reference refSeq > "$LOG_DIR/nanopore_cDNA_NA12878/refSeq/run.log" 2>&1

mkdir -p "$LOG_DIR/pacbio_ENCFF450VAU/refSeq"
echo "✨ Running pipeline for PacBio data ENCFF450VAU..."
python src/pipeline.py --method long \
        --bam_file pacbio_ENCFF450VAU/ENCFF450VAU.sorted.bam \
            --data_name pacbio_ENCFF450VAU \
                --reference refSeq > "$LOG_DIR/pacbio_ENCFF450VAU/refSeq/run.log" 2>&1

mkdir -p "$LOG_DIR/pacbio_ENCFF694DIE/refSeq"
echo "✨ Running pipeline for PacBio data ENCFF694DIE..."
python src/pipeline.py --method long \
        --bam_file pacbio_ENCFF694DIE/ENCFF694DIE.sorted.bam \
            --data_name pacbio_ENCFF694DIE \
                --reference refSeq > "$LOG_DIR/pacbio_ENCFF694DIE/refSeq/run.log" 2>&1