import os 
class config:
    def __init__(self, method, bam_file, candidate_sites_file,data_name):
        self.window_size = 200
        self.min_mapq = 10
        self.candidate_method = method
        self.data_name = data_name
        self.candidate_sites_file = candidate_sites_file
        self.bam_file = bam_file#"data/NA12878-cDNA.sorted.bam"
        self.tss_output_file, self.tes_output_file = self.get_output_file()
        self.soft_clip_window = 10
        self.splice_site_window = 100
        self.coverage_window = 100
        self.density_window = 100
        self.normalize = False
        self.create_output_files()
    
    def get_output_file(self):
        return f"features/{self.data_name}/{self.candidate_method}_tss.csv", f"features/{self.data_name}/{self.candidate_method}_tes.csv"
    
    def create_output_files(self):
        os.makedirs(os.path.dirname(self.tss_output_file), exist_ok=True)
        os.makedirs(os.path.dirname(self.tes_output_file), exist_ok=True)