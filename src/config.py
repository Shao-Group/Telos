class config:
    def __init__(self, method):
        self.window_size = 200
        self.min_mapq = 10
        self.candidate_method = method
        self.candidate_sites_file = f"data/{self.candidate_method}.tsstes"
        self.bam_file = "data/NA12878-cDNA.sorted.bam"
        self.tss_output_file, self.tes_output_file = self.get_output_file()
        self.soft_clip_window = 10
        self.splice_site_window = 100
        self.coverage_window = 100
        self.density_window = 100
        self.normalize = False
    
    def get_output_file(self):
        return f"features/{self.candidate_method}_tss.csv", f"features/{self.candidate_method}_tes.csv"
        