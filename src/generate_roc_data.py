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
        self.multi_exon_count_train = 0
        self.multi_exon_count_val = 0

        self.gtfformat_home = Path(self.cfg.rnaseqtools_dir) / "gtfformat"
        self.gtfcuff_home   = Path(self.cfg.rnaseqtools_dir) / "gtfcuff"
        print("Initializing directories...")
        # print("Multi exon count train:", self.multi_exon_count_train)
        # print("Multi exon count val:", self.multi_exon_count_val)
        self.ref_anno_train = None
        self.ref_anno_val = None


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
        
        pred_file_train = os.path.join(self.cfg.predictions_output_dir, f"{model}_stage2_predictions_train.csv")
        pred_file_val = os.path.join(self.cfg.predictions_output_dir, f"{model}_stage2_predictions_val.csv")
        pred_file_train = Path(pred_file_train).absolute()
        pred_file_val = Path(pred_file_val).absolute()
        out_file_train = Path(self.cfg.updated_cov_dir) / f"{model}-updated-cov-train.gtf"
        out_file_val = Path(self.cfg.updated_cov_dir) / f"{model}-updated-cov-val.gtf"
        out_file_train = out_file_train.absolute()
        out_file_val = out_file_val.absolute()

        logging.info(f"[{model}] update-cov")
        self._run(
            ["./gtfformat", "update-transcript-cov", str(self.baseline_gtf_train), str(pred_file_train), str(out_file_train)],
            cwd=self.gtfformat_home
        )
        self._run(
            ["./gtfformat", "update-transcript-cov", str(self.baseline_gtf_val), str(pred_file_val), str(out_file_val)],
            cwd=self.gtfformat_home
        )
        return out_file_train, out_file_val

    def run_gffcompare(self, gtf: Path, label: str, cwd: Path, is_train: bool):
        """Run gffcompare inside the conda env."""
        logging.info(f"[{label}] gffcompare")
        ref_anno = Path(self.ref_anno_train if is_train else self.ref_anno_val).absolute()
        cmd = self._gffcompare_cmd("-r", str(ref_anno), "-o", label, str(gtf))
        self._run(cmd, cwd=cwd)

    def generate_roc(self, label: str):
        """Run gtfcuff roc on a .tmap → .roc file."""
        tmap_train = Path(self.cfg.updated_cov_dir) / f"{label}.{label}-updated-cov-train.gtf.tmap"
        tmap_val = Path(self.cfg.updated_cov_dir) / f"{label}.{label}-updated-cov-val.gtf.tmap"
        roc_train = Path(self.cfg.transcript_pr_data) / f"{label}-updated-cov-train.roc"
        roc_val = Path(self.cfg.transcript_pr_data) / f"{label}-updated-cov-val.roc"
        logging.info(f"[{label}] gtfcuff pr → {roc_train.name} {roc_val.name}")
        with roc_train.open("w") as out_fh:
            self._run(
                ["./gtfcuff", "roc", str(tmap_train.absolute()), str(self.multi_exon_count_train), "cov"],
                cwd=self.gtfcuff_home,
                stdout=out_fh
            )
        with roc_val.open("w") as out_fh:
            self._run(
                ["./gtfcuff", "roc", str(tmap_val.absolute()), str(self.multi_exon_count_val), "cov"],
                cwd=self.gtfcuff_home,
                stdout=out_fh
            )
    
    def get_aupr(self, label: str, tmap: Path, is_train: bool):
        """Run gtfcuff roc on a .tmap → auc value"""
        logging.info(f"[{label}] gtfcuff auc → {self.cfg.auc_file_train if is_train else self.cfg.auc_file_val}")
        p = self._run(
            ["./gtfcuff", "auc", str(tmap), str(self.multi_exon_count_train if is_train else self.multi_exon_count_val)],
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
        with open(self.cfg.auc_file_train if is_train else self.cfg.auc_file_val, "a") as f:
            f.write(f"{label},{auc}\n")


    def process_model(self, model: str):
        """Full pipeline for one tool/model pair."""
        
        updated_gtf_train, updated_gtf_val = self.update_coverage(model)
        self.run_gffcompare(updated_gtf_train, model, cwd=self.cfg.gffcompare_dir_train, is_train=True)
        self.run_gffcompare(updated_gtf_val, model, cwd=self.cfg.gffcompare_dir_val, is_train=False)
        self.generate_roc(model)

        tmap_train = Path(self.cfg.updated_cov_dir) / f"{model}.{model}-updated-cov-train.gtf.tmap"
        tmap_val = Path(self.cfg.updated_cov_dir) / f"{model}.{model}-updated-cov-val.gtf.tmap"
        self.get_aupr(model, tmap_train.absolute(), is_train=True)
        self.get_aupr(model, tmap_val.absolute(), is_train=False)

    def process_all(self):
        """Run build + all tool/model combinations + baseline."""

        with open(self.cfg.auc_file_train, "w") as f:
            f.write("label,auc\n")

        with open(self.cfg.auc_file_val, "w") as f:
            f.write("label,auc\n")

        self.generate_baseline()

        for model_type in ["xgboost", "randomforest"]:
            self.process_model(model_type)
        
        self.process_baseline()

    def process_baseline(self):
        """Run gffcompare & ROC for the original baseline gtfs."""
        logging.info("[baseline] processing original GTFs")
        label_train = "baseline_train"
        label_val = "baseline_val"
        self.run_gffcompare(self.baseline_gtf_train, label_train, cwd=self.cfg.gffcompare_dir_train, is_train=True)
        self.run_gffcompare(self.baseline_gtf_val, label_val, cwd=self.cfg.gffcompare_dir_val, is_train=False)
       
        tmap_train = Path(self.cfg.data_output_dir) / f"{label_train}.{label_train}-chrom.gtf.tmap"
        tmap_val = Path(self.cfg.data_output_dir) / f"{label_val}.{label_val}-chrom.gtf.tmap"
        roc_out_file_train = Path(self.cfg.transcript_pr_data) / f"{label_train}.roc"
        roc_out_file_val = Path(self.cfg.transcript_pr_data) / f"{label_val}.roc"
        
        tmap_train = tmap_train.absolute()
        tmap_val = tmap_val.absolute()
        roc_out_file_train = roc_out_file_train.absolute()
        roc_out_file_val = roc_out_file_val.absolute()

        with roc_out_file_train.open("w") as out_fh:
            self._run(
                ["./gtfcuff", "roc", str(tmap_train), str(self.multi_exon_count_train), "cov"],
                cwd=self.gtfcuff_home,
                stdout=out_fh
            )
        with roc_out_file_val.open("w") as out_fh:
            self._run(
                ["./gtfcuff", "roc", str(tmap_val), str(self.multi_exon_count_val), "cov"],
                cwd=self.gtfcuff_home,
                stdout=out_fh
            )
        self.get_aupr(label_train, tmap_train, is_train=True)           
        self.get_aupr(label_val, tmap_val, is_train=False)           
    
    def generate_baseline(self):
        """Filter the baseline GTFs with validation chromosomes."""
        logging.info("[baseline] filtering original GTFs")
        
        input_gtf = Path(self.cfg.gtf_file).absolute()
        self.baseline_gtf_train = Path(self.cfg.data_output_dir) / f"baseline_train-chrom.gtf"
        self.baseline_gtf_val = Path(self.cfg.data_output_dir) / f"baseline_val-chrom.gtf"
        self.baseline_gtf_train = self.baseline_gtf_train.absolute()
        self.baseline_gtf_val = self.baseline_gtf_val.absolute()
        val_chromosome_file = Path(self.cfg.validation_chromosomes_file).absolute()
        train_chromosome_file = Path(self.cfg.train_chromosomes_file).absolute()

        self.ref_anno_train = Path(self.cfg.data_output_dir) / f"ref_anno-chrom-train.gtf"
        self.ref_anno_val = Path(self.cfg.data_output_dir) / f"ref_anno-chrom-val.gtf"
        self.ref_anno_train = self.ref_anno_train.absolute()
        self.ref_anno_val = self.ref_anno_val.absolute()

        self._run(
            ["./gtfformat", "filter-chrom", str(input_gtf), str(train_chromosome_file), str(self.baseline_gtf_train)],
            cwd=self.gtfformat_home
        )
        self._run(
            ["./gtfformat", "filter-chrom", str(input_gtf), str(val_chromosome_file), str(self.baseline_gtf_val)],
            cwd=self.gtfformat_home
        )
        self._run(
            ["./gtfformat", "filter-chrom", str(Path(self.cfg.ref_anno).absolute()), str(train_chromosome_file), str(self.ref_anno_train)],
            cwd=self.gtfformat_home
        )
        self._run(
            ["./gtfformat", "filter-chrom", str(Path(self.cfg.ref_anno).absolute()), str(val_chromosome_file), str(self.ref_anno_val)],
            cwd=self.gtfformat_home
        )
        self.multi_exon_count_train = self._count_multi_exon_transcripts(is_train=True)
        self.multi_exon_count_val = self._count_multi_exon_transcripts(is_train=False)
        
    def _count_multi_exon_transcripts(self, is_train: bool) -> int:
        """
        Parse the reference GTF and return the number of transcripts
        that have more than one 'exon' entry.
        """
        exon_counts = defaultdict(int)
        with open(self.ref_anno_train if is_train else self.ref_anno_val, 'r') as gtf:
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
    logging.basicConfig(level=logging.INFO)
    main()