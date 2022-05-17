#%%
from eqtl2gwas_pleiotropy.PathManager import PathManager

import matplotlib.pyplot as plt
import numpy
import os
import pandas
import pathlib

from eqtl2gwas_pleiotropy.constants import label_fontsize, tick_fontsize

plt.rcParams["figure.figsize"] = (8, 6)

# from matplotlib.pyplot import figure

# figure(figsize=(8, 6), dpi=80)
#%% Output
if not '__file__' in locals():
    __file__ = "plt_hist_gwas_etissue_egene.py"
outdir_path = os.path.join(PathManager.get_project_path(), "out", os.path.basename(__file__))
pathlib.Path(outdir_path).mkdir(parents=True, exist_ok=True)

#%%
indir_path = os.path.join(PathManager.get_project_path(), "out", "cmpt_count_per_rsid.py")
gwas_count_tsv_path = os.path.join(indir_path, "count_per_rsid_gwas.tsv")
etissue_count_tsv_path = os.path.join(indir_path, "count_per_rsid_etissue.tsv")
egene_count_tsv_path = os.path.join(indir_path, "count_per_rsid_egene.tsv")

#%%
ylabel = "Prob. Density"
title = "Colocalized variants{}"
ylim=[1e-10, 1]
edgecolor='k'

linewidth = 2

hist_kwargs = {'density': 1, 'edgecolor': edgecolor, 'linewidth': linewidth}

#%%
egene_count_df = pandas.read_csv(egene_count_tsv_path, sep="\t", header=0)
bins = (numpy.array(range(8)))/2*10
# print(bins)
ax = egene_count_df['egene_count'].hist(**hist_kwargs, bins=bins)
ax.set_xlabel("# eQTL genes", fontsize=label_fontsize)
ax.set_ylabel(ylabel, fontsize=label_fontsize)
ax.set_yscale('log')
ax.set_ylim(ylim)
plt.xticks(fontsize=tick_fontsize)
plt.yticks(fontsize=tick_fontsize)
plt.title(title.format(" and eQTL genes"), fontsize=label_fontsize)
fig = ax.get_figure()
png_path = os.path.join(outdir_path, "hist_egene.png")
fig.savefig(png_path)
plt.clf()
plt.close()

#%%

etissue_count_df = pandas.read_csv(etissue_count_tsv_path, sep="\t", header=0)
bins = (numpy.array(range(11)))*10
# print(bins)
ax = etissue_count_df['etissue_subcategory_count'].hist(**hist_kwargs, bins=bins)
ax.set_xlabel("# eQTL samples", fontsize=label_fontsize)
ax.set_ylabel(ylabel, fontsize=label_fontsize)
ax.set_yscale('log')
ax.set_ylim(ylim)
plt.xticks(fontsize=tick_fontsize)
plt.yticks(fontsize=tick_fontsize)
plt.title(title.format(" and eQTL tissues"), fontsize=label_fontsize)
fig = ax.get_figure()
png_path = os.path.join(outdir_path, "hist_etissue.png")
fig.savefig(png_path)
plt.clf()
plt.close()

#%%
gwas_count_df = pandas.read_csv(gwas_count_tsv_path, sep="\t", header=0)
bins = numpy.array(range(6))
ax = gwas_count_df['gwas_subcategory_count'].hist(**hist_kwargs, bins=bins)

ax.set_yscale('log')
ax.set_xlabel("# GWAS Categories", fontsize=label_fontsize)
ax.set_ylabel(ylabel, fontsize=label_fontsize)
# ax.set_xlim([0, 6])
ax.set_ylim(ylim)
plt.xticks(fontsize=tick_fontsize)
plt.yticks(fontsize=tick_fontsize)
plt.title(title.format(" and GWAS categories"), fontsize=label_fontsize)

fig = ax.get_figure()
png_path = os.path.join(outdir_path, "hist_gwas.png")
fig.savefig(png_path)
plt.clf()
plt.close()
