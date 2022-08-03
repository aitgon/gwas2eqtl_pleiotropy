import sys

import seaborn

from eqtl2gwas_pleiotropy.PathManager import PathManager
from matplotlib import pyplot as plt
from eqtl2gwas_pleiotropy.constants import tick_fontsize, label_fontsize, scatter_dot_size, dpi

import os
import pandas
import pathlib

plt.rcParams["figure.figsize"] = (8, 6)
from eqtl2gwas_pleiotropy.constants import seaborn_theme_dic
seaborn.set_theme(**seaborn_theme_dic)

#%%
help_cmd_str = "todo"
try:
    region_window_100000_tsv_path = sys.argv[1]
    count_per_rsid_egene_tsv_path = sys.argv[2]
    outdir_path = sys.argv[3]
    if len(sys.argv) > 4:
        print("""Two many arguments!
        {}""".format(help_cmd_str))
        sys.exit(1)
except IndexError:
    print("""Argument missing!
    {}""".format(help_cmd_str))
    sys.exit(1)

#%% input dir cmpt_count_per_rsid
basename_str = "count_per_rsid_egene.tsv"
# indir_path = os.path.join(PathManager.get_project_path(), "out", "cmpt_count_per_rsid.py")
# tsv_path = os.path.join(indir_path, basename_str)
count_per_rsid_df = pandas.read_csv(count_per_rsid_egene_tsv_path, sep="\t")

#%% input regions
# region_window_100000_tsv_path = os.path.join(PathManager.get_outdir_path(), "cmpt_pleiotropic_regions.py", "region_window_100000.tsv")
region_window_100000_df = pandas.read_csv(region_window_100000_tsv_path, sep="\t")

# #%% Output
# if not '__file__' in locals():
#     __file__ = "pltsctr_x_per_rsid_y_egene.py"
# outdir_path = os.path.join(PathManager.get_project_path(), "out", os.path.basename(__file__))
pathlib.Path(outdir_path).mkdir(parents=True, exist_ok=True)

title = "Count Per SNP: eGenes"
count_col_name = "egene_count"
ylim = [0, 35]
c = 'blue'

count_per_rsid_df['pos'] = count_per_rsid_df['pos'].astype('int')

#%% Loop over regions
pleiotropic_regions_df = region_window_100000_df.loc[region_window_100000_df['gwas_category_count'] >= 6, ['chrom', 'start', 'end', 'gwas_category_count']]
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
    # plt.scatter(count_per_rsid_gwas_region_df['pos'] / 1000000, count_per_rsid_gwas_region_df[count_col_name], c='blue', s=scatter_dot_size)
    seaborn.scatterplot(x=count_per_rsid_gwas_region_df['pos'] / 1000000, y=count_per_rsid_gwas_region_df[count_col_name], s=scatter_dot_size)

    plt.grid(True)
    plt.title(title, fontsize=label_fontsize)
    plt.xlabel("Chr{} position [Mbp]".format(chrom), fontsize=label_fontsize)
    plt.xticks(fontsize=tick_fontsize, rotation=45)
    plt.ylabel("eGene count", fontsize=label_fontsize)
    plt.ylim(ylim)
    plt.yticks(fontsize=tick_fontsize)

    plt.tight_layout()
    png_path = os.path.join(outdir_path, "count_per_rsid_chr{}_start{}_end{}_categories{}.png".format(chrom, start, end, gwas_category_count))
    plt.savefig(png_path, dpi=dpi)
    plt.clf()
    plt.close()
