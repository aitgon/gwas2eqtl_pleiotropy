"""Scatter of etissue category count vs variant rsid"""

from matplotlib import pyplot as plt
from gwas2eqtl_pleiotropy.constants import tick_fontsize, label_fontsize, scatter_dot_size, dpi

import sys
import seaborn
import os
import pandas
import pathlib

plt.rcParams["figure.figsize"] = (8, 6)
from gwas2eqtl_pleiotropy.constants import seaborn_theme_dic
seaborn.set_theme(**seaborn_theme_dic)

#%%
help_cmd_str = "todo"
try:
    max_gwas_category_count = int(sys.argv[1])
    region_window_100000_ods_path = sys.argv[2]
    count_per_rsid_gwas_ods_path = sys.argv[3]
    outdir_path = sys.argv[4]
    if len(sys.argv) > 5:
        print("""Two many arguments!
        {}""".format(help_cmd_str))
        sys.exit(1)
except IndexError:
    print("""Argument missing!
    {}""".format(help_cmd_str))
    sys.exit(1)

#%% input dir cmpt_count_per_rsid
# basename_str = "count_per_rsid_etissue.tsv"
# indir_path = os.path.join(PathManager.get_project_path(), "out", "cmpt_count_per_rsid.py")
# tsv_path = os.path.join(indir_path, basename_str)
# count_per_rsid_df = pandas.read_csv(count_per_rsid_etissue_tsv_path, sep="\t")
count_per_rsid_df = pandas.read_excel(count_per_rsid_gwas_ods_path, engine='odf')

#%% input regions
# region_window_100000_tsv_path = os.path.join(PathManager.get_outdir_path(), "cmpt_pleiotropic_regions.py", "region_window_100000.tsv")
# region_window_100000_df = pandas.read_csv(region_window_100000_tsv_path, sep="\t")
region_window_100000_df = pandas.read_excel(region_window_100000_ods_path, engine='odf')

# #%% Output
# if not '__file__' in locals():
#     __file__ = "pltsctr_x_per_rsid_y_etissue.py"
# outdir_path = os.path.join(PathManager.get_project_path(), "out", os.path.basename(__file__))
pathlib.Path(outdir_path).mkdir(parents=True, exist_ok=True)

title = "Coloc. eQTL/GWAS variant"
count_col_name = "etissue_category_term_count"
ylim = [0, 40]
c = 'blue'
ylabel = "eTissue count"

count_per_rsid_df['pos38'] = count_per_rsid_df['pos38'].astype('int')

#%% Loop over regions
pleiotropic_regions_df = region_window_100000_df.loc[
    region_window_100000_df['gwas_category_count'] >= max_gwas_category_count, ['chrom', 'start', 'end', 'gwas_category_count']]
for rowi, row in pleiotropic_regions_df.iterrows():
    chrom = row['chrom']
    start = row['start']
    end = row['end']
    gwas_category_count = row['gwas_category_count']
    count_per_rsid_gwas_region_df = count_per_rsid_df.copy()
    count_per_rsid_gwas_region_df = count_per_rsid_gwas_region_df.loc[count_per_rsid_gwas_region_df['chrom'] == chrom, ]
    count_per_rsid_gwas_region_df = count_per_rsid_gwas_region_df.loc[count_per_rsid_gwas_region_df['pos38'] >= start, ]
    count_per_rsid_gwas_region_df = count_per_rsid_gwas_region_df.loc[count_per_rsid_gwas_region_df['pos38'] <= end, ]

    #%%
    # plt.scatter(count_per_rsid_gwas_region_df['pos38'] / 1000000, count_per_rsid_gwas_region_df[count_col_name], c='blue', s=scatter_dot_size)

    seaborn.scatterplot(x=count_per_rsid_gwas_region_df['pos38'] / 1000000,
                        y=count_per_rsid_gwas_region_df[count_col_name], s=scatter_dot_size)

    plt.grid(True)
    plt.xticks(fontsize=tick_fontsize, rotation=45)
    plt.yticks(fontsize=tick_fontsize)
    plt.title(title, fontsize=label_fontsize)
    plt.xlabel("Chr{} position [Mbp]".format(chrom), fontsize=label_fontsize)
    plt.ylabel(ylabel, fontsize=label_fontsize)
    plt.ylim(ylim)

    plt.tight_layout()
    png_path = os.path.join(outdir_path, "count_per_rsid_chr{}_start{}_end{}_categories{}.png".format(chrom, start, end, gwas_category_count))
    plt.savefig(png_path, dpi=dpi)
    plt.clf()
    plt.close()
