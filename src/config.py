class config:
    def __init__(self):
        self.window_size = 200
        self.min_mapq = 10
        self.candidate_sites_file = "data/isoquant.tsstes"
        self.bam_file = "data/NA12878-cDNA.sorted.bam"
        self.tss_output_file, self.tes_output_file = self.get_output_file()
    
    def get_output_file(self):
        candidate_method = self.candidate_sites_file.split('/')[-1].split('.')[0]
        return f"features/{candidate_method}_tss.csv", f"features/{candidate_method}_tes.csv"
        