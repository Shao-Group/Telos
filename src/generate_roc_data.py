#!/usr/bin/env python3
import logging
import subprocess
import os
from pathlib import Path
from dataclasses import dataclass, field
from typing import List
import re
from collections import defaultdict
import logging
from config import Config
import argparse

class ROCPipeline:
    def __init__(self, cfg: Config, gffcompare_env: str):
        self.cfg = cfg
        self.gffcompare_env = gffcompare_env
        # ensure output dirs exist
        self.multi_exon_count = self._count_multi_exon_transcripts()

        self.gtfformat_home = Path(self.cfg.rnaseqtools_dir) / "gtfformat"
        self.gtfcuff_home   = Path(self.cfg.rnaseqtools_dir) / "gtfcuff"
        print("Initializing directories...")
        print("Multi exon count:", self.multi_exon_count)
        self.baseline_gtf = None


    def _run(self, cmd: List[str], cwd: Path = None, **kwargs):
        """Run a subprocess with logging and error checking."""
        logging.info(f"→ {' '.join(cmd)}  (cwd={cwd})")
        p = subprocess.run(cmd, cwd=cwd, check=True, **kwargs)
        return p
 
    def _gffcompare_cmd(self, *args) -> List[str]:
        """Wrap gffcompare in a conda-run command."""
        return ["conda", "run", "-n", self.gffcompare_env, "gffcompare", *args]

    def build_gtfformat(self):
        """Clean + make the gtfformat tool."""
        self._run(["make", "clean"], cwd=self.cfg.home)
        self._run(["make"],       cwd=self.cfg.home)

    def update_coverage(self, model: str):
        """Invoke gtfformat update-cov for one (tool,model)."""
        
        pred_file = os.path.join(self.cfg.predictions_output_dir, f"{model}_stage2_predictions.csv")
        pred_file = Path(pred_file).absolute()
        out_file = Path(self.cfg.updated_cov_dir) / f"{model}-updated-cov.gtf"
        out_file = out_file.absolute()

        logging.info(f"[{model}] update-cov")
        self._run(
            ["./gtfformat", "update-transcript-cov", str(self.baseline_gtf), str(pred_file), str(out_file)],
            cwd=self.gtfformat_home
        )
        return out_file

    def run_gffcompare(self, gtf: Path, label: str, cwd: Path):
        """Run gffcompare inside the conda env."""
        logging.info(f"[{label}] gffcompare")
        ref_anno = Path(self.cfg.ref_anno).absolute()
        cmd = self._gffcompare_cmd("-r", str(ref_anno), "-o", label, str(gtf))
        self._run(cmd, cwd=cwd)

    def generate_roc(self, label: str):
        """Run gtfcuff roc on a .tmap → .roc file."""
        tmap = Path(self.cfg.updated_cov_dir) / f"{label}.{label}-updated-cov.gtf.tmap"
        roc  = Path(self.cfg.transcript_pr_data) / f"{label}-updated-cov.roc"
        logging.info(f"[{label}] gtfcuff pr → {roc.name}")
        with roc.open("w") as out_fh:
            self._run(
                ["./gtfcuff", "roc", str(tmap.absolute()), str(self.multi_exon_count), "cov"],
                cwd=self.gtfcuff_home,
                stdout=out_fh
            )
    
    def get_aupr(self, label: str, tmap: Path):
        """Run gtfcuff roc on a .tmap → auc value"""
        logging.info(f"[{label}] gtfcuff auc → {self.cfg.auc_file}")
        p = self._run(
            ["./gtfcuff", "auc", str(tmap), str(self.multi_exon_count)],
            cwd=self.gtfcuff_home,
            capture_output=True,
            text=True
        )
        print(p.stdout)
        match = re.search(r"auc\s*=\s*(\d+\.?\d*)", p.stdout)
        assert match, f"Failed to parse AUC from output: {p.stdout}"    
        auc = float(match.group(1))
        logging.info(f"[{label}] auc = {auc}")
        # append to the auc file
        with open(self.cfg.auc_file, "a") as f:
            f.write(f"{label},{auc}\n")


    def process_model(self, model: str):
        """Full pipeline for one tool/model pair."""
        
        updated_gtf = self.update_coverage(model)
        self.run_gffcompare(updated_gtf, model, cwd=self.cfg.gffcompare_dir)
        self.generate_roc(model)

        tmap = Path(self.cfg.updated_cov_dir) / f"{model}.{model}-updated-cov.gtf.tmap"
        self.get_aupr(model, tmap.absolute())

    def process_all(self):
        """Run build + all tool/model combinations + baseline."""
        
        with open(self.cfg.auc_file, "w") as f:
            f.write("label,auc\n")

        self.generate_baseline()

        for model_type in ["xgboost", "randomforest"]:
            self.process_model(model_type)
        
        self.process_baseline()

    def process_baseline(self):
        """Run gffcompare & ROC for the original baseline gtfs."""
        logging.info("[baseline] processing original GTFs")
        label = "baseline"
        self.run_gffcompare(self.baseline_gtf, label, cwd=self.cfg.gffcompare_dir)
       
        tmap = Path(self.cfg.data_output_dir) / f"{label}.{label}-chrom-filtered.gtf.tmap"
        roc_out_file = Path(self.cfg.transcript_pr_data) / f"{label}.roc"
        
        tmap = tmap.absolute()
        roc_out_file = roc_out_file.absolute()

        with roc_out_file.open("w") as out_fh:
            self._run(
                ["./gtfcuff", "roc", str(tmap), str(self.multi_exon_count), "cov"],
                cwd=self.gtfcuff_home,
                stdout=out_fh
            )
        self.get_aupr(label, tmap)           
    
    def generate_baseline(self):
        """Filter the baseline GTFs with validation chromosomes."""
        logging.info("[baseline] filtering original GTFs")
        
        input_gtf = Path(self.cfg.gtf_file).absolute()
        self.baseline_gtf = Path(self.cfg.data_output_dir) / f"baseline-chrom-filtered.gtf"
        self.baseline_gtf = self.baseline_gtf.absolute()
        val_chromosome_file = Path(self.cfg.validation_chromosomes_file).absolute()
        self._run(
            ["./gtfformat", "filter-chrom", str(input_gtf), str(val_chromosome_file), str(self.baseline_gtf)],
            cwd=self.gtfformat_home
        )
        
    def _count_multi_exon_transcripts(self) -> int:
        """
        Parse the reference GTF and return the number of transcripts
        that have more than one 'exon' entry.
        """
        exon_counts = defaultdict(int)
        with open(self.cfg.ref_anno, 'r') as gtf:
            for line in gtf:
                if line.startswith('#'):
                    continue
                cols = line.rstrip('\n').split('\t')
                if cols[2] != 'exon':
                    continue
                attrs = cols[8]
                m = re.search(r'transcript_id\s+"([^"]+)"', attrs)
                if m:
                    exon_counts[m.group(1)] += 1

        multi = sum(1 for cnt in exon_counts.values() if cnt > 1)
        logging.info(f"Reference has {multi} multi-exon transcripts")
        return multi




def main():
    parser = argparse.ArgumentParser(description="Generate ROC data for RNAseq tools.")
    parser.add_argument("--project-config", required=True, help="Path to the project config file.")
    parser.add_argument("--gffcompare-env", required=True, help="Conda environment for gffcompare.")

    p = parser.parse_args()
    # Parse command line arguments
    cfg = Config.load_from_file(p.project_config)
    gffcompare_env = p.gffcompare_env

    # Create and run the pipeline
    pipeline = ROCPipeline(cfg, gffcompare_env)
    pipeline.process_all()

    

if __name__ == "__main__":
    # logging.basicConfig(level=logging.INFO)
    main()