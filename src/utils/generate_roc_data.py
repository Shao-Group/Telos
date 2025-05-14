#!/usr/bin/env python3
import logging
import subprocess
import sys
from pathlib import Path
from dataclasses import dataclass, field
from typing import List
import re
from collections import defaultdict
import logging

@dataclass
class PipelineConfig:
    data_name: str
    ref_anno: Path
    val_chrom: Path
    anno_name: str 
    pred_dir: Path = field(init = False) 
    auc_file: Path = field(init = False)
    # Repositories / tools
    home: Path = Path("/datadisk1/ixk5174/tools/rnaseqtools/gtfformat")
    home_cuff: Path = Path("/datadisk1/ixk5174/tools/rnaseqtools/gtfcuff")
    gffcompare_env: str = "gffcompare"

    # Data directories
    data_home: Path = Path("/datadisk1/ixk5174/long_reads_compare/out/gffcomp-results/")
    data_dir: Path = field(init=False)
    project_out_dir: Path = Path("/datadisk1/ixk5174/tools/tss-tes-project/out/")
    project_data_dir: Path = Path("/datadisk1/ixk5174/tools/tss-tes-project/data/")
    out_dir: Path = field(init=False)
    roc_out_dir: Path = field(init=False)
    # ref_anno: Path = Path("/datadisk1/ixk5174/long_reads_compare/anno/refSeq_anno.gtf")
    # val_chrom: Path = Path(f"/datadisk1/ixk5174/tools/tss-tes-project/data_train/")
    # Parameters
    tools: List[str]  = field(default_factory=lambda: ["stringtie", "isoquant"])
    models: List[str] = field(default_factory=lambda: ["xgboost", "randomforest"])

    def __post_init__(self):
        self.data_home = self.data_home / self.data_name
        self.data_dir    = self.data_home / "val-baseline"
        self.out_dir     = self.data_home / "updated-cov"
        self.roc_out_dir = self.out_dir  / "roc"
        self.pred_dir = self.project_out_dir / self.data_name / self.anno_name / "predictions/transcripts"
        self.project_data_dir = self.project_data_dir / self.data_name
        self.ref_anno = self.ref_anno.absolute()
        self.val_chrom = self.val_chrom.absolute()
        self.auc_file = self.pred_dir / f"auc.csv"

class TSSPipeline:
    def __init__(self, cfg: PipelineConfig):
        self.cfg = cfg
        # ensure output dirs exist
        self.cfg.out_dir.mkdir(parents=True, exist_ok=True)
        self.cfg.roc_out_dir.mkdir(parents=True, exist_ok=True)
        self.multi_exon_count = self._count_multi_exon_transcripts()


    def _run(self, cmd: List[str], cwd: Path = None, **kwargs):
        """Run a subprocess with logging and error checking."""
        logging.info(f"→ {' '.join(cmd)}  (cwd={cwd})")
        proc = subprocess.run(cmd, cwd=cwd, check=True, **kwargs)
        # print(proc.stdout)
        # print(proc.stderr)

    def _gffcompare_cmd(self, *args) -> List[str]:
        """Wrap gffcompare in a conda-run command."""
        return ["conda", "run", "-n", self.cfg.gffcompare_env, "gffcompare", *args]

    def build_gtfformat(self):
        """Clean + make the gtfformat tool."""
        self._run(["make", "clean"], cwd=self.cfg.home)
        self._run(["make"],       cwd=self.cfg.home)

    def update_coverage(self, tool: str, model: str, universe: bool = False):
        """Invoke gtfformat update-cov for one (tool,model)."""
        data_file = self.cfg.data_dir / f"{tool}-chrom-filtered.gtf"
        pred_name = f"{f'{tool}_universe_' if universe else f'{tool}_'}{model}_transcript_predictions.tsv"
        pred_file = self.cfg.pred_dir / pred_name

        suffix = f"universe-{model}" if universe else model
        out_file = self.cfg.out_dir / f"{tool}-{suffix}-updated-cov.gtf"

        logging.info(f"[{tool}/{suffix}] update-cov")
        self._run(
            ["./gtfformat", "update-transcript-cov", str(data_file), str(pred_file), str(out_file)],
            cwd=self.cfg.home
        )
        return out_file, f"{tool}-{suffix}"

    def run_gffcompare(self, gtf: Path, label: str, cwd: Path):
        """Run gffcompare inside the conda env."""
        logging.info(f"[{label}] gffcompare")
        cmd = self._gffcompare_cmd("-r", str(self.cfg.ref_anno), "-o", label, str(gtf))
        self._run(cmd, cwd=cwd)

    def generate_roc(self, label: str):
        """Run gtfcuff roc on a .tmap → .roc file."""
        tmap = self.cfg.out_dir / f"{label}.{label}-updated-cov.gtf.tmap" 
        roc  = self.cfg.roc_out_dir / f"{label}-updated-cov.roc"
        logging.info(f"[{label}] gtfcuff roc → {roc.name}")
        with roc.open("w") as out_fh:
            self._run(
                ["./gtfcuff", "roc", str(tmap), str(self.multi_exon_count), "cov"],
                cwd=self.cfg.home_cuff,
                stdout=out_fh
            )
    
    def get_aupr(self, label: str):
        """Run gtfcuff roc on a .tmap → auc value"""
        tmap = self.cfg.out_dir / f"{label}.{label}-updated-cov.gtf.tmap" 
        tmp_file = self.cfg.out_dir / f"{label}-updated-cov.auc"
        logging.info(f"[{label}] gtfcuff auc → {tmp_file.name}")
        with tmp_file.open("w") as out_fh:
            self._run(
                ["./gtfcuff", "auc", str(tmap), str(self.multi_exon_count)],
                cwd=self.cfg.home_cuff,
                stdout=out_fh
            )
        
        # parse the auc value
        with open(tmp_file, "r") as f:
            match = re.search(r"auc\s*=\s*(\d+\.?\d*)", f.read())
            if match:
                auc = float(match.group(1))
                logging.info(f"[{label}] auc = {auc}")
            else:
                auc = None
                exit(55)
        
        # append to the auc file
        with self.cfg.auc_file.open("a") as f:
            f.write(f"{label},{auc}\n")


    def process_model(self, tool: str, model: str):
        """Full pipeline for one tool/model pair."""
        # 1) update-cov for standard & universe
        for is_universe in (False, True):
            gtf, label = self.update_coverage(tool, model, universe=is_universe)
            # 2) gffcompare
            self.run_gffcompare(gtf, label, cwd=self.cfg.out_dir)
            # 3) ROC
            self.generate_roc(label)
            self.get_aupr(label)

    def process_all(self):
        """Run build + all tool/model combinations + baseline."""
        # self.build_gtfformat()
        self.generate_baseline()

        for tool in self.cfg.tools:
            for model in self.cfg.models:
                self.process_model(tool, model)

        self.process_baseline()

    def process_baseline(self):
        """Run gffcompare & ROC for the original baseline gtfs."""
        logging.info("[baseline] processing original GTFs")
        for tool in self.cfg.tools:
            gtf   = self.cfg.data_dir / f"{tool}-chrom-filtered.gtf"
            label = f"{tool}-baseline"

            # gffcompare
            self.run_gffcompare(gtf, label, cwd=self.cfg.data_dir)
            # ROC
            # self.generate_roc(label)
            tmap = self.cfg.data_dir / f"{label}.{tool}-chrom-filtered.gtf.tmap"
            roc_out_file = self.cfg.roc_out_dir / f"{label}.roc"
            # print(f"roc: {out_fh}")
            with roc_out_file.open("w") as out_fh:
                self._run(
                    ["./gtfcuff", "roc", str(tmap), str(self.multi_exon_count), "cov"],
                    cwd=self.cfg.home_cuff,
                    stdout=out_fh
                )
            
            # AUC
            tmp_file = self.cfg.out_dir / f"{label}.auc"    
            logging.info(f"[{label}] gtfcuff auc → {tmp_file.name}")
            with tmp_file.open("w") as out_fh:
                self._run(
                    ["./gtfcuff", "auc", str(tmap), str(self.multi_exon_count)],
                    cwd=self.cfg.home_cuff,
                    stdout=out_fh
                )
            # parse the auc value
            with open(tmp_file, "r") as f:
                match = re.search(r"auc\s*=\s*(\d+\.?\d*)", f.read())
                if match:
                    auc = float(match.group(1))
                    logging.info(f"[{label}] auc = {auc}")
                else:
                    auc = None
                    exit(55)
            # append to the auc file
            with self.cfg.auc_file.open("a") as f:
                f.write(f"{label},{auc}\n")
            
    
    def generate_baseline(self):
        """Filter the baseline GTFs with validation chromosomes."""
        logging.info("[baseline] filtering original GTFs")
        self.cfg.data_dir.mkdir(parents=True, exist_ok=True)
        for tool in self.cfg.tools:
            gtf = self.cfg.project_data_dir / f"{tool}.gtf"
            out_gtf = self.cfg.data_dir / f"{tool}-chrom-filtered.gtf"
            self._run(
                ["./gtfformat", "filter-chrom", str(gtf), str(self.cfg.val_chrom), str(out_gtf)],
                cwd=self.cfg.home
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




def main(data_name: str, tools: List[str], annotation: str, val_chrom_file: str, anno_name: str):

    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    cfg = PipelineConfig(data_name = data_name, tools=tools, ref_anno=Path(annotation), val_chrom=Path(val_chrom_file), anno_name= anno_name)
    pipeline = TSSPipeline(cfg)
    try:
        pipeline.process_all()
        return cfg.roc_out_dir
    except subprocess.CalledProcessError as e:
        logging.error(f"Step failed (exit code {e.returncode}).")
        sys.exit(e.returncode)

