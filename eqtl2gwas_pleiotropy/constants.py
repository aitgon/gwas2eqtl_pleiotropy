import os
from eqtl2gwas_pleiotropy.PathManager import PathManager

# Raw colocalization data
coloc_all_path = os.path.join(PathManager.get_project_path(), "out", "coloc_all")
coloc_raw_tsv_path = os.path.join(coloc_all_path, "coloc_all_20220329.tsv")

public_coloc_all_tsv_path = os.path.join(coloc_all_path, "coloc_all_20220329.tsv")
coloc_h4_tsv_path = os.path.join(PathManager.get_project_path(), "out", "filter_h4.py/coloc_h4.tsv")

# OpenGWAS info
# opengwas_metadata_url = 'http://gwas-api.mrcieu.ac.uk/gwasinfo'

# EBI eQTL public metadate
eqtl_metadata_url = 'https://raw.githubusercontent.com/eQTL-Catalogue/eQTL-Catalogue-resources/master/tabix/tabix_ftp_paths.tsv'

h4_cutoff = 0.8  # coloc cutoff

region_bin = 100000  # region bin to define as pleiotropic

# Pleiotropic regions with 5 categories
pleiotropic_regions_5 = [
    [2, 27508073, 27519736],
    [3, 141375367, 141428572],
    [11, 61783884, 61855668],
    [16, 28743363, 28908020],
]

# Pleiotropic regions with 4 categories
pleiotropic_regions_4 = [
    [5, 132254564, 132485305],
    [6, 31629923, 32460012],
    [16, 28478025, 28620700],
    [19, 48700572, 48714803],
]

# Pyplot constants
label_fontsize = 20
tick_fontsize = 10
scatter_dot_size = 50
