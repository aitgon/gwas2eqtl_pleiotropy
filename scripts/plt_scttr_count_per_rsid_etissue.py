from eqtl2gwas_pleiotropy.PathManager import PathManager
from matplotlib import pyplot as plt

import matplotlib
import os
import pandas
import pathlib


#%% Parameters
from eqtl2gwas_pleiotropy.constants import pleiotropic_regions_5, label_fontsize, tick_fontsize, scatter_dot_size

if not '__file__' in locals():
    __file__ = "plt_scttr_count_per_rsid_etissue.py"
indir_path = os.path.join(PathManager.get_project_path(), "out", "cmpt_count_per_rsid.py")
outdir_path = os.path.join(PathManager.get_project_path(), "out", os.path.basename(__file__))
tsv_path = os.path.join(indir_path, "count_per_rsid_etissue.tsv")
size = 10
title = "Count Per SNP: Tissues"
count_col_name = "etissue_subcategory_count"
ylim = [0, 80]

#%% Paths
pathlib.Path(outdir_path).mkdir(parents=True, exist_ok=True)

#%% Data

count_per_rsid_gwas_df = pandas.read_csv(tsv_path, sep="\t")
count_per_rsid_gwas_df['pos'] = count_per_rsid_gwas_df['pos'].astype('int')

#%% Loop over regions
for chrom, start, end in pleiotropic_regions_5:
    #%% scatter 16:28 528 527-28 904 206
    # chrom = 16
    # start = 28500000
    # end = 29000000
    count_per_rsid_gwas_region_df = count_per_rsid_gwas_df.copy()
    count_per_rsid_gwas_region_df = count_per_rsid_gwas_region_df.loc[count_per_rsid_gwas_region_df['chrom'] == chrom, ]
    count_per_rsid_gwas_region_df = count_per_rsid_gwas_region_df.loc[count_per_rsid_gwas_region_df['pos'] >= start, ]
    count_per_rsid_gwas_region_df = count_per_rsid_gwas_region_df.loc[count_per_rsid_gwas_region_df['pos'] <= end, ]

    #%%
    plt.scatter(count_per_rsid_gwas_region_df['pos'] / 1000000, count_per_rsid_gwas_region_df[count_col_name], c='blue', s=scatter_dot_size)
    plt.grid(True)
    plt.xticks(fontsize=tick_fontsize)
    plt.yticks(fontsize=tick_fontsize)
    plt.title(title, fontsize=label_fontsize)
    plt.xlabel("Position Chr{} [Mbp]".format(chrom), fontsize=label_fontsize)
    plt.ylabel("GWAS Disease Category Count", fontsize=label_fontsize)
    plt.ylim(ylim)
    png_path = os.path.join(outdir_path, "count_per_rsid_etissue_chr{}_start{}_end{}_categories5.png".format(chrom, start, end))
    plt.savefig(png_path)
    plt.close()
