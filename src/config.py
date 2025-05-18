import os
import pickle

class Config:
    _instance = None

    def __init__(self, bam_file, gtf_file, prefix, output_dir, rnaseqtools_dir, ref_anno_gtf, tmap_file):
        # your existing defaults
        self.window_size = 200
        self.min_mapq = 10
        self.soft_clip_window = 10
        self.splice_site_window = 100
        self.coverage_window = 100
        self.density_window = 100
        self.normalize = False

        # parameters
        self.data_name       = prefix
        self.output_dir      = output_dir
        self.rnaseqtools_dir = rnaseqtools_dir
        self.bam_file        = bam_file
        self.gtf_file        = gtf_file
        self.ref_anno        = ref_anno_gtf
        self.tmap_file        = tmap_file

        # placeholders
        self.cov_file           = None
        self.candidate_file     = None
        self.ref_candidate_file = None
        self.validation_chromosomes_file = None

        self._make_dirs()

    def _make_dirs(self):
        # create structured output directories
        base = self.output_dir
        self.data_output_dir        = os.path.join(base, "data")
        self.features_output_dir    = os.path.join(base, "features")
        self.reports_output_dir     = os.path.join(base, "reports")
        self.predictions_output_dir = os.path.join(base, "predictions")
        self.models_output_dir      = os.path.join(base, "models")
        self.updated_cov_dir        = os.path.join(base, "updated_cov")

        self.transcript_pr_data     = os.path.join(self.reports_output_dir, "transcript_pr_data")
        self.pr_data_dir            = os.path.join(self.reports_output_dir, "pr_data")
        self.feature_importance_dir = os.path.join(self.reports_output_dir, "feature_importance")
        self.gffcompare_dir        = os.path.join(self.reports_output_dir, "gffcompare")

        for d in (self.data_output_dir, self.features_output_dir,
                  self.reports_output_dir, self.predictions_output_dir,
                  self.models_output_dir, self.pr_data_dir,
                  self.feature_importance_dir, self.gffcompare_dir, 
                  self.transcript_pr_data, self.updated_cov_dir):
            os.makedirs(d, exist_ok=True)


        # feature / label files
        p = self.data_name
        self.tss_feature_file   = os.path.join(self.features_output_dir, f"{p}_tss_features.tsv")
        self.tes_feature_file   = os.path.join(self.features_output_dir, f"{p}_tes_features.tsv")
        self.tss_labeled_file   = os.path.join(self.features_output_dir, f"{p}_tss_labeled.tsv")
        self.tes_labeled_file   = os.path.join(self.features_output_dir, f"{p}_tes_labeled.tsv")
        self.auc_file           = os.path.join(self.transcript_pr_data, f"{p}_auc.csv")

    @classmethod
    def create(cls, bam_file, gtf_file, prefix, output_dir, rnaseqtools_dir, ref_anno_gtf, tmap_file):
        """Initialize singleton (only once)."""
        if cls._instance is None:
            cls._instance = cls(bam_file, gtf_file, prefix, output_dir, rnaseqtools_dir, ref_anno_gtf, tmap_file)
        return cls._instance

    @classmethod
    def get(cls):
        """Retrieve the existing singleton."""
        if cls._instance is None:
            raise ValueError("Config not created yet. Call Config.create(...) first.")
        return cls._instance

    @classmethod
    def reset(cls):
        cls._instance = None

    @classmethod
    def save_to_file(cls, path: str):
        """Pickle the singleton out to disk."""
        cfg = cls.get()
        with open(path, 'wb') as f:
            pickle.dump(cfg, f, protocol=pickle.HIGHEST_PROTOCOL)
    
    @classmethod
    def load_from_file(cls, path: str):
        """Unpickle from disk and set as the singleton."""
        with open(path, 'rb') as f:
            cfg = pickle.load(f)
        cls._instance = cfg
        return cfg

    # your existing setters…
    def set_cov_file(self, cov_file):           self.cov_file = cov_file
    def set_candidate_file(self, candidate):    self.candidate_file = candidate
    def set_ref_candidate_file(self, ref_cand): self.ref_candidate_file = ref_cand

# Convenience aliases
create_config            = Config.create
get_config               = Config.get
reset_config             = Config.reset
save_config              = Config.save_to_file
load_config              = Config.load_from_file
set_cov_file             = Config.set_cov_file
set_candidate_file       = Config.set_candidate_file
set_ref_candidate_file   = Config.set_ref_candidate_file
