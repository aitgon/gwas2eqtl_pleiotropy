from eqtl2gwas_pleiotropy.PathManager import PathManager
from eqtl2gwas_pleiotropy.URL import URL
from pybedtools import BedTool

import matplotlib.pyplot as plt
import os
import pathlib


#%%
from eqtl2gwas_pleiotropy.constants import label_fontsize

plt.rcParams["figure.figsize"] = (8, 6)
ylabel = "Prob. Density"

#%% Input
variant_pleio_flank_0_bed_path = os.path.join(PathManager.get_outdir_path(), "cmpt_count_per_rsid.py", "variant_pleio_{}_flank_{}.bed")

#%% Output
if not '__file__' in locals():
    outdir_path = os.path.join(PathManager.get_outdir_path(), "plt_hist_remap_crm.py")
else:
    outdir_path = os.path.join(PathManager.get_outdir_path(), os.path.basename(__file__))
pathlib.Path(outdir_path).mkdir(parents=True, exist_ok=True)

url_str = "http://remap.univ-amu.fr/storage/remap2022/hg38/MACS2/remap2022_crm_macs2_hg38_v1_0.bed.gz"
remap_crm_path = URL(url_str, data_public_dir="/home/gonzalez/Software/public").download()

#%%
ylim = [0, 0.07]
xlim = [0, 1000]
bins = 50
# flank = 50  # 0 or 50
remap_crm_bed_obj = BedTool(remap_crm_path)

for flank in [0, 50]:
    pleioprev_count_lst = \
    remap_crm_bed_obj.intersect(BedTool(variant_pleio_flank_0_bed_path.format(1, flank)), wb=True).to_dataframe(
        names=range(0, 16))[3].str.split(',').apply(len)
    for pleio_i in range(2, 6):
        pleioi_count_lst = remap_crm_bed_obj.intersect(BedTool(variant_pleio_flank_0_bed_path.format(pleio_i, flank)), wb=True).to_dataframe(names=range(0, 16))[3].str.split(',').apply(len)
        plt.ylim(ylim)
        plt.xlim(xlim)
        plt.hist(pleioprev_count_lst, alpha=0.5, label='pleio {}'.format(1), density=True, bins=bins)
        plt.hist(pleioi_count_lst, alpha=0.5, label='pleio {}'.format(pleio_i), density=True, bins=bins)
        plt.grid(axis='y')
        plt.legend(loc='upper right', fontsize=label_fontsize)
        png_path = os.path.join(outdir_path, "remap_crm_pleio{}_flank_{}_hist.png".format(pleio_i, flank))
        plt.savefig(png_path)
        plt.clf()
        plt.close()
