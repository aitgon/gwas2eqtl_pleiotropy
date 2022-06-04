from eqtl2gwas_pleiotropy.PathManager import PathManager
from matplotlib import pyplot as plt
from eqtl2gwas_pleiotropy.constants import tick_fontsize, label_fontsize, scatter_dot_size

import os
import pandas
import pathlib

#%% input dir cmpt_count_per_rsid
basename_str = "count_per_rsid_egene.tsv"
indir_path = os.path.join(PathManager.get_project_path(), "out", "cmpt_count_per_rsid.py")
tsv_path = os.path.join(indir_path, basename_str)
count_per_rsid_df = pandas.read_csv(tsv_path, sep="\t")

#%% input regions
region_window_100000_tsv_path = os.path.join(PathManager.get_outdir_path(), "cmpt_pleiotropic_regions.py", "region_window_100000.tsv")
region_window_100000_df = pandas.read_csv(region_window_100000_tsv_path, sep="\t")

#%% Output
if not '__file__' in locals():
    __file__ = "plt_scttr_count_per_rsid_egene.py"
outdir_path = os.path.join(PathManager.get_project_path(), "out", os.path.basename(__file__))
pathlib.Path(outdir_path).mkdir(parents=True, exist_ok=True)

title = "Count Per SNP: eGenes"
count_col_name = "egene_count"
ylim = [0, 35]
c = 'blue'

count_per_rsid_df['pos'] = count_per_rsid_df['pos'].astype('int')

#%% Loop over regions
pleiotropic_regions_df = region_window_100000_df.loc[region_window_100000_df['gwas_category_count'] >= 5, ['chrom', 'start', 'end', 'gwas_category_count']]
for rowi, row in pleiotropic_regions_df.iterrows():
    chrom = row['chrom']
    start = row['start']
    end = row['end']
    gwas_category_count = row['gwas_category_count']
    count_per_rsid_gwas_region_df = count_per_rsid_df.copy()
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
    plt.ylabel("eGene Count", fontsize=label_fontsize)
    plt.ylim(ylim)
    png_path = os.path.join(outdir_path, "count_per_rsid_chr{}_start{}_end{}_categories{}.png".format(chrom, start, end, gwas_category_count))
    plt.savefig(png_path)
    plt.close()
