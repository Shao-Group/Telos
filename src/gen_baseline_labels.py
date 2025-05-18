#!/usr/bin/env python3
import logging
import subprocess
import sys
from pathlib import Path
from dataclasses import dataclass, field
from typing import List
import argparse
from config import Config

@dataclass
class TranscriptLabelingConfig:
    rnaseq_dir: Path
    ref_anno: Path    
    data_file: Path
    out_dir: Path 
    prefix: str 
    # Repositories / tools
    home_gtfformat: Path = field(init=False)

    def __post_init__(self):        
        self.home_gtfformat = self.rnaseq_dir / "gtfformat"
        self.ref_anno = self.ref_anno.absolute()
        self.data_file = self.data_file.absolute()
        self.home_gtfformat = self.home_gtfformat.absolute()
        

class TranscriptLabelingPipeline:
    def __init__(self, cfg: TranscriptLabelingConfig):
        self.cfg = cfg
        # ensure output dirs exist
        self.cfg.out_dir.mkdir(parents=True, exist_ok=True)
        

    def _run(self, cmd: List[str], cwd: Path = None, **kwargs):
        """Run a subprocess with logging and error checking."""
        logging.info(f"→ {' '.join(cmd)}  (cwd={cwd})")
        subprocess.run(cmd, cwd=cwd, check=True, **kwargs)
    
    def get_coverage(self):
        """Get coverage for the tool."""

        cov_tsv = self.cfg.out_dir / f"{self.cfg.prefix}_cov.tsv"
        self._run(
            ["./gtfformat", "get-cov", str(self.cfg.data_file), str(cov_tsv.absolute())],
            cwd=self.cfg.home_gtfformat
        )
        return str(cov_tsv)
    
    def generate_candidate_sites(self):
        """Generate candidate sites for the tool."""
        candidate_file = self.cfg.out_dir / f"{self.cfg.prefix}_candidates.tsv"
        with open (candidate_file, "w") as cf:
            self._run(
                ["./gtfformat", "TSSTES", str(self.cfg.data_file)], stdout=cf,
                cwd=self.cfg.home_gtfformat
            )
        ref_candidate_file = self.cfg.out_dir / f"_ref_candidates.tsv"
        with open(ref_candidate_file, "w") as rcf:
            self._run(
                ["./gtfformat", "TSSTES", str(self.cfg.ref_anno)], stdout=rcf,
                cwd=self.cfg.home_gtfformat
            )
        
        return str(candidate_file), str(ref_candidate_file)

    def process_model(self):
        """Full pipeline for one tool/model pair."""
        cov_file = self.get_coverage()
        candidate_file, ref_candidate_file = self.generate_candidate_sites()
        return cov_file, candidate_file, ref_candidate_file




def main(project_config):
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    cfg = TranscriptLabelingConfig(rnaseq_dir=Path(project_config.rnaseqtools_dir), data_file=Path(project_config.gtf_file), 
                                   ref_anno=Path(project_config.ref_anno), out_dir=Path(project_config.data_output_dir),
                                   prefix=project_config.data_name)
    pipeline = TranscriptLabelingPipeline(cfg)
    try:
        return pipeline.process_model()
        # return cfg.roc_out_dir
    except subprocess.CalledProcessError as e:
        logging.error(f"Step failed (exit code {e.returncode}).")
        sys.exit(e.returncode)

