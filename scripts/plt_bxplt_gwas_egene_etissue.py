#%%
import pathlib

from statannot import add_stat_annotation

from eqtl2gwas_pleiotropy.PathManager import PathManager

import os
import pandas
import matplotlib.pyplot as plt
import seaborn as sns

#%% Output
if not '__file__' in locals():
    __file__ = "plt_bxplt_gwas_egene_etissue.py"
outdir_path = os.path.join(PathManager.get_project_path(), "out", os.path.basename(__file__))
pathlib.Path(outdir_path).mkdir(parents=True, exist_ok=True)

#%% input
indir_path = os.path.join(PathManager.get_project_path(), "out", "cmpt_count_per_rsid.py")
egene_count_tsv_path = os.path.join(indir_path, "count_per_rsid_egene.tsv")
etissue_count_tsv_path = os.path.join(indir_path, "count_per_rsid_etissue.tsv")
gwas_count_tsv_path = os.path.join(indir_path, "count_per_rsid_gwas.tsv")


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
ax = sns.boxplot(x="gwas_subcategory_count", y="egene_count", data=merged_df, order=order)
test_results = add_stat_annotation(ax, data=merged_df, x="gwas_subcategory_count", y="egene_count", order=order,
                                   box_pairs=[(1, 2), (1, 3), (1, 4), (1, 5)],
                                   test='Mann-Whitney', text_format='star',
                                   loc='inside', verbose=2)
fig = ax.get_figure()
png_path = os.path.join(outdir_path, "boxplot_gwas_egene.png")
fig.savefig(png_path)
plt.clf()
plt.close()

#%%
sns.set_theme(style="whitegrid")
ax = sns.boxplot(x="gwas_subcategory_count", y="etissue_subcategory_count", data=merged_df, order=order)
test_results = add_stat_annotation(ax, data=merged_df, x="gwas_subcategory_count", y="egene_count", order=order,
                                   box_pairs=[(1, 2), (1, 3), (1, 4), (1, 5)],
                                   test='Mann-Whitney', text_format='star',
                                   loc='inside', verbose=2)
fig = ax.get_figure()
png_path = os.path.join(outdir_path, "boxplot_gwas_etissue.png")
fig.savefig(png_path)
plt.clf()
plt.close()
