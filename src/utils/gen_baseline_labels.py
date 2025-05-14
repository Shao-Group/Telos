#!/usr/bin/env python3
import logging
import subprocess
import sys
from pathlib import Path
from dataclasses import dataclass, field
from typing import List


@dataclass
class TranscriptLabelingConfig:
    data_name: str
    ref_anno: Path 
    pred_dir: Path = field(init = False) 
    # Repositories / tools
    home: Path = Path("/datadisk1/ixk5174/tools/rnaseqtools/gtfformat")
    home_cuff: Path = Path("/datadisk1/ixk5174/tools/rnaseqtools/gtfcuff")
    gffcompare_env: str = "gffcompare"

    # Data directories
    
    data_home: Path = Path("/datadisk1/ixk5174/tools/tss-tes-project/data/")
    out_dir: Path = field(init=False)
    
    # Parameters
    tools: List[str]  = field(default_factory=lambda: ["stringtie", "isoquant"])
    models: List[str] = field(default_factory=lambda: ["xgboost", "randomforest"])

    def __post_init__(self):
        self.data_home = self.data_home / self.data_name
        self.out_dir     = self.data_home / "gffcmp"
        self.ref_anno = self.ref_anno.absolute()
        # self.val_chrom = self.val_chrom / f"{self.anno_name}_validation_chromosomes.txt"    

class TranscriptLabelingPipeline:
    def __init__(self, cfg: TranscriptLabelingConfig):
        self.cfg = cfg
        # ensure output dirs exist
        self.cfg.out_dir.mkdir(parents=True, exist_ok=True)
        # self.cfg.roc_out_dir.mkdir(parents=True, exist_ok=True)

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

    def get_coverage(self, tool: str):
        """Get coverage for the tool."""
        out_gtf = self.cfg.data_home / f"{tool}.gtf"
        cov_tsv = self.cfg.data_home / f"{tool}-cov.tsv"
        self._run(
            ["./gtfformat", "get-cov", str(out_gtf), str(cov_tsv)],
            cwd=self.cfg.home
        )
    
    def generate_candidate_sites(self, tool: str):
        """Generate candidate sites for the tool."""
        input_gtf = self.cfg.data_home / f"{tool}.gtf"
        candidate_file = self.cfg.data_home / f"{self.cfg.data_name}_{tool}_candidates.tsv"
        with open (candidate_file, "w") as cf:
            self._run(
                ["./gtfformat", "TSSTES", str(input_gtf)], stdout=cf,
                cwd=self.cfg.home
            )

    def run_gffcompare(self, gtf: Path, label: str, cwd: Path):
        """Run gffcompare inside the conda env."""
        logging.info(f"[{label}] gffcompare")
        cmd = self._gffcompare_cmd("-r", str(self.cfg.ref_anno), "-o", label, str(gtf))
        self._run(cmd, cwd=cwd)

    def process_model(self, tool: str):
        """Full pipeline for one tool/model pair."""
        out_gtf = self.cfg.data_home / f"{tool}.gtf"
        # 1) update-cov for standard & universe
        if tool == "isoquant":
            # isoquant has a different format for the GTF
            # we need to convert it to the standard GTF format
            
            tpm_file = self.cfg.data_home / f"{tool}.tpm"
            self._run(
                ["./gtfformat", "update-tpm", str(out_gtf), str(tpm_file), str(out_gtf)],
                cwd=self.cfg.home
            )

        self.get_coverage(tool)
        self.generate_candidate_sites(tool)
        # 2) gffcompare
        self.run_gffcompare(out_gtf, f"{tool}", cwd=self.cfg.out_dir)
           

    def process_all(self):
        """Run build + all tool/model combinations + baseline."""
        for tool in self.cfg.tools:
            self.process_model(tool)



def main(data_name: str, tools: List[str], reference: str):
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    cfg = TranscriptLabelingConfig(data_name = data_name, tools=tools, ref_anno=Path(reference))
    pipeline = TranscriptLabelingPipeline(cfg)
    try:
        pipeline.process_all()
        # return cfg.roc_out_dir
    except subprocess.CalledProcessError as e:
        logging.error(f"Step failed (exit code {e.returncode}).")
        sys.exit(e.returncode)



if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python transcript_labeling.py <data_name>")
        sys.exit(1)
    data_name = sys.argv[1]
    main(data_name)