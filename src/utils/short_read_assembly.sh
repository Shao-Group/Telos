#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

eval "$(conda shell.bash hook)"

# ─── CONFIG ─────────────────────────────────────────────────────────────────────
INPUT_DIR="/datadisk1/shared/data/encode10"
OUT_DIR="/datadisk1/ixk5174/long_reads_compare/out/sr-out-refseq"
SAMTOOLS_ENV="samtools-irtesam"
SCALLOP2_ENV="scallop2"
STRINGTIE_BIN="/datadisk1/ixk5174/tools/stringtie-3.0.0.Linux_x86_64/stringtie"
THREADS=4

# ─── MAIN LOOP ───────────────────────────────────────────────────────────────────
mkdir -p "$OUT_DIR"

for sample_dir in "$INPUT_DIR"/SRR*; do
  sample=$(basename "$sample_dir")
  for aligner in hisat star; do
    bam="$sample_dir/${aligner}.sort.bam"
    if [[ ! -f "$bam" ]]; then
      echo "⚠️  Skipping: $bam not found"
      continue
    fi

    out_subdir="$OUT_DIR/${sample}_${aligner}"
    bam_sorted="$out_subdir/${aligner}.sorted.bam"

    mkdir -p "$out_subdir"
    
    if [[ -f "$bam_sorted" ]]; then
      echo "⚠️  Skipping: $bam_sorted already exists"
      continue
    fi
    # 1) Sort & index with samtools
    echo "🔄 Sorting & indexing $sample/$aligner..."
    cp "$bam" "$bam_sorted"
    conda activate "$SAMTOOLS_ENV"
    samtools index      "$bam_sorted"
    conda deactivate

    echo "✨ Launching StringTie & Scallop2 in parallel for $sample/$aligner …"
    # continue
    # 2) Run both tools in the background
    (
    # StringTie run
    "$STRINGTIE_BIN" \
        -o "$out_subdir/stringtie.gtf" \
        "$bam_sorted" \
    && echo "✅ StringTie done for $sample/$aligner"
    ) &

    (
    # Scallop2 run (needs conda)
    conda activate "$SCALLOP2_ENV"
    scallop2 \
        -i "$bam_sorted" \
        -o "$out_subdir/scallop2.gtf" > "$out_subdir/scallop2.log" \
    && echo "✅ Scallop2 done for $sample/$aligner"
    conda deactivate
    ) &

    # Wait for both to finish before continuing
    wait
    
    echo "🏁 Completed both assemblers for $sample/$aligner"
    echo
  done
done

echo "🏁 All done!"
