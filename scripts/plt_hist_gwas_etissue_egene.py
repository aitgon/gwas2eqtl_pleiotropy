#%%
import pdb
import sys

from eqtl2gwas_pleiotropy.PathManager import PathManager

import matplotlib.pyplot as plt
import numpy
import os
import pandas
import pathlib

from eqtl2gwas_pleiotropy.constants import label_fontsize, tick_fontsize, dpi

plt.rcParams["figure.figsize"] = (8, 6)

#%%
help_cmd_str = "todo"
try:
    gwas_count_tsv_path = sys.argv[1]
    egene_count_tsv_path = sys.argv[2]
    etissue_count_tsv_path = sys.argv[3]
    hist_rsid_gwas_path = sys.argv[4]
    hist_rsid_egene_path = sys.argv[5]
    hist_rsid_etissue_path = sys.argv[6]
    if len(sys.argv) > 7:
        print("""Two many arguments!
        {}""".format(help_cmd_str))
        sys.exit(1)
except IndexError:
    print("""Argument missing!
    {}""".format(help_cmd_str))
    sys.exit(1)

pathlib.Path(os.path.dirname(hist_rsid_gwas_path)).mkdir(parents=True, exist_ok=True)
pathlib.Path(os.path.dirname(hist_rsid_egene_path)).mkdir(parents=True, exist_ok=True)
pathlib.Path(os.path.dirname(hist_rsid_etissue_path)).mkdir(parents=True, exist_ok=True)

# from matplotlib.pyplot import figure

# figure(figsize=(8, 6), dpi=80)
# #%% Output
# if not '__file__' in locals():
#     __file__ = "plt_hist_gwas_etissue_egene.py"
# outdir_path = os.path.join(PathManager.get_project_path(), "out", os.path.basename(__file__))
# pathlib.Path(outdir_path).mkdir(parents=True, exist_ok=True)

#%%
# indir_path = os.path.join(PathManager.get_project_path(), "out", "cmpt_count_per_rsid.py")
# gwas_count_tsv_path = os.path.join(indir_path, "count_per_rsid_gwas.tsv")
# egene_count_tsv_path = os.path.join(indir_path, "count_per_rsid_egene.tsv")
# etissue_count_tsv_path = os.path.join(indir_path, "count_per_rsid_etissue.tsv")

#%%
ylabel = "Probability density"
title = "Variants{}"
ylim=[1e-10, 10]
edgecolor='k'
linewidth = 2
hist_kwargs = {'density': 1, 'edgecolor': edgecolor, 'linewidth': linewidth}

#%% gwas
count_df = pandas.read_csv(gwas_count_tsv_path, sep="\t", header=0)
bins = numpy.array(range(6))
data_ser = count_df['gwas_category_count']
plt.hist(data_ser, **hist_kwargs, bins=bins)

plt.grid(True)
plt.title(title.format(" and GWAS categories"), fontsize=label_fontsize)
plt.xlabel("GWAS category count", fontsize=label_fontsize)
plt.xticks(fontsize=tick_fontsize)
plt.xticks(fontsize=tick_fontsize, rotation=0)
plt.ylabel(ylabel, fontsize=label_fontsize)
plt.ylim(ylim)
plt.yscale('log')
plt.yticks(fontsize=tick_fontsize)

plt.tight_layout()
# png_path = os.path.join(outdir_path, "hist_gwas.png")
plt.savefig(hist_rsid_gwas_path, dpi=dpi)
plt.clf()
plt.close()

#%% egene
count_df = pandas.read_csv(egene_count_tsv_path, sep="\t", header=0)
data_ser = count_df['egene_count']
plt.hist(data_ser, **hist_kwargs)

plt.grid(True)
plt.title(title.format(" and eGenes"), fontsize=label_fontsize)
plt.xlabel("eGene count", fontsize=label_fontsize)
plt.xticks(fontsize=tick_fontsize)
plt.xticks(fontsize=tick_fontsize, rotation=0)
plt.ylabel(ylabel, fontsize=label_fontsize)
plt.ylim(ylim)
plt.yscale('log')
plt.yticks(fontsize=tick_fontsize)

plt.tight_layout()
# png_path = os.path.join(outdir_path, "hist_egene.png")
plt.savefig(hist_rsid_egene_path, dpi=dpi)
plt.clf()
plt.close()

#%% etissue
count_df = pandas.read_csv(etissue_count_tsv_path, sep="\t", header=0)
data_ser = count_df['etissue_subcategory_count']
plt.hist(data_ser, **hist_kwargs)

plt.grid(True)
plt.title(title.format(" and eTissues"), fontsize=label_fontsize)
plt.xlabel("eTissue subcategory count", fontsize=label_fontsize)
plt.xticks(fontsize=tick_fontsize)
plt.xticks(fontsize=tick_fontsize, rotation=0)
plt.ylabel(ylabel, fontsize=label_fontsize)
plt.ylim(ylim)
plt.yscale('log')
plt.yticks(fontsize=tick_fontsize)

plt.tight_layout()
# png_path = os.path.join(outdir_path, "hist_etissue.png")
plt.savefig(hist_rsid_etissue_path, dpi=dpi)
plt.clf()
plt.close()
