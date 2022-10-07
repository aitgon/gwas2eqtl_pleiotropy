import os
from gwas2eqtl_pleiotropy.PathManager import PathManager

public_data_dir = os.path.join(os.environ['HOME'], "Software", "public")

# Raw colocalization data
coloc_raw_tsv_path = os.path.join(PathManager.get_outdir_path(), "coloc_all/genome/5e-08/1000000/coloc.tsv")
coloc_all_path = os.path.join(PathManager.get_project_path(), "out", "coloc_all")

public_coloc_all_tsv_path = os.path.join(coloc_all_path, "coloc_all_20220329.tsv")
coloc_h4_tsv_path = os.path.join(PathManager.get_project_path(), "out", "filter_h4.py/coloc_h4.tsv")

# EBI eQTL public metadate
eqtl_metadata_url = 'https://raw.githubusercontent.com/eQTL-Catalogue/eQTL-Catalogue-resources/master/tabix/tabix_ftp_paths.tsv'

h4_cutoff = 0.8  # coloc cutoff
snp_h4_cutoff = 0.5  # pp cutoff of snp
region_bin = 100000  # region bin to define as pleiotropic

# seaborn theme
seaborn_theme_dic = {'style': "darkgrid"}

# Pyplot constants
label_fontsize = 30
tick_fontsize = 24
scatter_dot_size = 50
dpi = 150  # for publication 600, for ongoing 150

# Boxplot constants
alpha = 0.5
boxplot_kwargs = {'linewidth': 3, 'notch': True, 'palette': "rocket_r", 'showfliers': False, 'boxprops': dict(alpha=alpha)}

# statannotations
annotator_config_dic = {'fontsize': 16}