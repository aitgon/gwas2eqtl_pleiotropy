import os
from gwas2eqtl_pleiotropy.PathManager import PathManager

public_data_dir = os.path.join(os.environ['HOME'], "Software", "public")

# EBI eQTL public metadate
eqtl_metadata_url = 'https://raw.githubusercontent.com/eQTL-Catalogue/eQTL-Catalogue-resources/master/tabix/tabix_ftp_paths.tsv'

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

# mhc locus, excluded in gwas2eqtl
mhc_chrom = 6
mhc_start38 = 25000000
mhc_end38 = 35000000
