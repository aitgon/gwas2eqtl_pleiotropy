from eqtl2gwas_pleiotropy.PathManager import PathManager
from matplotlib import pyplot as plt
from eqtl2gwas_pleiotropy.constants import pleiotropic_regions_5, tick_fontsize, label_fontsize, pleiotropic_regions_4, \
    scatter_dot_size

import matplotlib
import os
import pandas
import pathlib

#%% Output

if not '__file__' in locals():
    __file__ = "plt_scttr_count_per_rsid_gwas.py"
outdir_path = os.path.join(PathManager.get_project_path(), "out", os.path.basename(__file__))
pathlib.Path(outdir_path).mkdir(parents=True, exist_ok=True)

#%% Parameters
indir_path = os.path.join(PathManager.get_project_path(), "out", "cmpt_count_per_rsid.py")
tsv_path = os.path.join(indir_path, "count_per_rsid_gwas.tsv")
title = "Count Per SNP: Disease Category"
count_col_name = "gwas_subcategory_count"
ylim = [0, 5.1]
c = 'blue'

#%% Paths
pathlib.Path(outdir_path).mkdir(parents=True, exist_ok=True)

count_per_rsid_gwas_df = pandas.read_csv(tsv_path, sep="\t")
count_per_rsid_gwas_df['pos'] = count_per_rsid_gwas_df['pos'].astype('int')

# font = {'size': 16}
# matplotlib.rc('font', **font)

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
    plt.grid(True)
    plt.xticks(fontsize=tick_fontsize)
    plt.yticks(fontsize=tick_fontsize)
    plt.scatter(count_per_rsid_gwas_region_df['pos'] / 1000000, count_per_rsid_gwas_region_df[count_col_name], c='blue', s=scatter_dot_size)
    plt.title(title, fontsize=label_fontsize)
    plt.xlabel("Position Chr{} [Mbp]".format(chrom), fontsize=label_fontsize)
    plt.ylabel("GWAS Disease Category Count", fontsize=label_fontsize)
    plt.ylim(ylim)
    png_path = os.path.join(outdir_path, "count_per_rsid_gwas_chr{}_start{}_end{}_categories5.png".format(chrom, start, end))
    plt.savefig(png_path)
    plt.close()

#%% Loop over regions
for chrom, start, end in pleiotropic_regions_4:
    #%% scatter 16:28 528 527-28 904 206
    # chrom = 16
    # start = 28500000
    # end = 29000000
    count_per_rsid_gwas_region_df = count_per_rsid_gwas_df.copy()
    count_per_rsid_gwas_region_df = count_per_rsid_gwas_region_df.loc[count_per_rsid_gwas_region_df['chrom'] == chrom, ]
    count_per_rsid_gwas_region_df = count_per_rsid_gwas_region_df.loc[count_per_rsid_gwas_region_df['pos'] >= start, ]
    count_per_rsid_gwas_region_df = count_per_rsid_gwas_region_df.loc[count_per_rsid_gwas_region_df['pos'] <= end, ]


    #%%
    plt.grid(True)
    plt.xticks(fontsize=tick_fontsize)
    plt.yticks(fontsize=tick_fontsize)
    plt.scatter(count_per_rsid_gwas_region_df['pos'] / 1000000, count_per_rsid_gwas_region_df[count_col_name], c='blue', s=scatter_dot_size)
    plt.title(title, fontsize=label_fontsize)
    plt.xlabel("Position Chr{} [Mbp]".format(chrom), fontsize=label_fontsize)
    plt.ylabel("GWAS Disease Category Count", fontsize=label_fontsize)
    plt.ylim(ylim)
    png_path = os.path.join(outdir_path, "count_per_rsid_gwas_chr{}_start{}_end{}_categories4.png".format(chrom, start, end))
    plt.savefig(png_path)
    plt.close()
