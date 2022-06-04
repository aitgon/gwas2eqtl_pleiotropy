from eqtl2gwas_pleiotropy.PathManager import PathManager
from matplotlib import pyplot as plt
from eqtl2gwas_pleiotropy.constants import tick_fontsize, label_fontsize, scatter_dot_size, dpi

import os
import pandas
import pathlib

plt.rcParams["figure.figsize"] = (8, 6)

#%% input dir cmpt_count_per_rsid
basename_str = "count_per_rsid_gwas.tsv"
indir_path = os.path.join(PathManager.get_project_path(), "out", "cmpt_count_per_rsid.py")
tsv_path = os.path.join(indir_path, basename_str)
count_per_rsid_df = pandas.read_csv(tsv_path, sep="\t")

#%% input regions
region_window_100000_tsv_path = os.path.join(PathManager.get_outdir_path(), "cmpt_pleiotropic_regions.py", "region_window_100000.tsv")
region_window_100000_df = pandas.read_csv(region_window_100000_tsv_path, sep="\t")

#%% Output
if not '__file__' in locals():
    __file__ = "plt_scttr_count_per_rsid_gwas.py"
outdir_path = os.path.join(PathManager.get_project_path(), "out", os.path.basename(__file__))
pathlib.Path(outdir_path).mkdir(parents=True, exist_ok=True)

#%% Plot parameters
title = "Count Per SNP: Disease Category"
count_col_name = "gwas_subcategory_count"
ylim = [0, 5.1]
c = 'blue'

count_per_rsid_df['pos'] = count_per_rsid_df['pos'].astype('int')

#%% Loop over regions
pleiotropic_regions_df = region_window_100000_df.loc[region_window_100000_df['gwas_category_count'] >= 5, ['chrom', 'start', 'end', 'gwas_category_count']]
for rowi, row in pleiotropic_regions_df.iterrows():
    chrom = row['chrom']
    start = row['start']
    end = row['end']
    gwas_category_count = row['gwas_category_count']
    #%% scatter 16:28 528 527-28 904 206
    # chrom = 16
    # start = 28500000
    # end = 29000000
    count_per_rsid_gwas_region_df = count_per_rsid_df.copy()
    count_per_rsid_gwas_region_df = count_per_rsid_gwas_region_df.loc[count_per_rsid_gwas_region_df['chrom'] == chrom, ]
    count_per_rsid_gwas_region_df = count_per_rsid_gwas_region_df.loc[count_per_rsid_gwas_region_df['pos'] >= start, ]
    count_per_rsid_gwas_region_df = count_per_rsid_gwas_region_df.loc[count_per_rsid_gwas_region_df['pos'] <= end, ]


    #%%
    plt.grid(True)
    plt.xticks(fontsize=tick_fontsize)
    plt.yticks(fontsize=tick_fontsize)
    plt.scatter(count_per_rsid_gwas_region_df['pos'] / 1000000, count_per_rsid_gwas_region_df[count_col_name], c='blue', s=scatter_dot_size)
    plt.title(title, fontsize=label_fontsize)
    plt.xlabel("Chr{} position [Mbp]".format(chrom), fontsize=label_fontsize)
    plt.ylabel("GWAS category ct.", fontsize=label_fontsize)
    plt.ylim(ylim)
    plt.tight_layout()
    png_path = os.path.join(outdir_path, "count_per_rsid_chr{}_start{}_end{}_categories{}.png".format(chrom, start, end, gwas_category_count))
    plt.savefig(png_path, dpi=dpi)
    plt.clf()
    plt.close()
