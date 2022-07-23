#%%
import sys

from eqtl2gwas_pleiotropy.PathManager import PathManager
from statannot import add_stat_annotation

import matplotlib.pyplot as plt
import os
import pandas
import pathlib
import seaborn as sns

from eqtl2gwas_pleiotropy.constants import label_fontsize, tick_fontsize, dpi

plt.rcParams["figure.figsize"] = (8, 6)

#%%
help_cmd_str = "todo"
try:
    gwas_count_tsv_path = sys.argv[1]
    gwas_count_tsv_path = sys.argv[2]
    egene_count_tsv_path = sys.argv[2]
    etissue_count_tsv_path = sys.argv[3]
    boxplot_gwas_egene_png_path = sys.argv[4]
    boxplot_gwas_etissue_png_path = sys.argv[5]
    if len(sys.argv) > 6:
        print("""Two many arguments!
        {}""".format(help_cmd_str))
        sys.exit(1)
except IndexError:
    print("""Argument missing!
    {}""".format(help_cmd_str))
    sys.exit(1)

# #%% Output
# if not '__file__' in locals():
#     __file__ = "plt_bxplt_gwas_egene_etissue.py"
# outdir_path = os.path.join(PathManager.get_project_path(), "out", os.path.basename(__file__))
pathlib.Path(os.path.dirname(boxplot_gwas_egene_png_path)).mkdir(parents=True, exist_ok=True)
pathlib.Path(os.path.dirname(boxplot_gwas_etissue_png_path)).mkdir(parents=True, exist_ok=True)

#%% input
# indir_path = os.path.join(PathManager.get_project_path(), "out", "cmpt_count_per_rsid.py")
# egene_count_tsv_path = os.path.join(indir_path, "count_per_rsid_egene.tsv")
# etissue_count_tsv_path = os.path.join(indir_path, "count_per_rsid_etissue.tsv")
# gwas_count_tsv_path = os.path.join(indir_path, "count_per_rsid_gwas.tsv")

#%%
gwas_count_df = pandas.read_csv(gwas_count_tsv_path, sep="\t", header=0)
egene_count_df = pandas.read_csv(egene_count_tsv_path, sep="\t", header=0)
etissue_count_df = pandas.read_csv(etissue_count_tsv_path, sep="\t", header=0)

#%%
merged_df = gwas_count_df.merge(egene_count_df, on=['chrom', 'pos', 'rsid'])
merged_df = merged_df.merge(etissue_count_df, on=['chrom', 'pos', 'rsid'])

#%%
sns.set_theme(style="whitegrid")
order = [1, 2, 3, 4, 5]
ax = sns.boxplot(x="gwas_category_count", y="egene_count", data=merged_df, order=order)
test_results = add_stat_annotation(ax, data=merged_df, x="gwas_category_count", y="egene_count", order=order,
                                   box_pairs=[(1, 2), (1, 3), (1, 4), (1, 5)],
                                   test='Mann-Whitney', text_format='star',
                                   loc='inside', verbose=2)

plt.title("GWAS category and eGene count", fontsize=label_fontsize)
plt.xlabel("GWAS category count", fontsize=label_fontsize)
plt.xticks(fontsize=tick_fontsize, rotation=0)
plt.ylabel("eTissue count", fontsize=label_fontsize)
plt.yticks(fontsize=tick_fontsize)

plt.tight_layout()
# png_path = os.path.join(outdir_path, "boxplot_gwas_egene.png")
fig = ax.get_figure()
fig.savefig(boxplot_gwas_egene_png_path)
plt.clf()
plt.close()

#%%
sns.set_theme(style="whitegrid")
ax = sns.boxplot(x="gwas_category_count", y="etissue_subcategory_count", data=merged_df, order=order)
test_results = add_stat_annotation(ax, data=merged_df, x="gwas_category_count", y="egene_count", order=order,
                                   box_pairs=[(1, 2), (1, 3), (1, 4), (1, 5)],
                                   test='Mann-Whitney', text_format='star',
                                   loc='inside', verbose=2)

plt.title("GWAS category and eTissue count", fontsize=label_fontsize)
plt.xlabel("GWAS category count", fontsize=label_fontsize)
plt.xticks(fontsize=tick_fontsize, rotation=0)
plt.ylabel("eGene count", fontsize=label_fontsize)
plt.yticks(fontsize=tick_fontsize)

plt.tight_layout()
fig = ax.get_figure()
# png_path = os.path.join(outdir_path, "boxplot_gwas_etissue.png")
fig.savefig(boxplot_gwas_etissue_png_path, dpi=dpi)
plt.clf()
plt.close()
