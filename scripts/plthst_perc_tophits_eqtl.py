#%%
import sys
import matplotlib.pyplot as plt
import os
import pandas
import pathlib
import seaborn

from gwas2eqtl_pleiotropy.constants import label_fontsize, tick_fontsize, dpi

plt.rcParams["figure.figsize"] = (8, 6)
from gwas2eqtl_pleiotropy.constants import seaborn_theme_dic
seaborn.set_theme(**seaborn_theme_dic)

#%%
help_cmd_str = "todo"
try:
    perc_tophits_eqtl_tsv = sys.argv[1]
    hist_perc_tophits_eqtl_png = sys.argv[2]
    if len(sys.argv) > 4:
        print("""Two many arguments!
        {}""".format(help_cmd_str))
        sys.exit(1)
except IndexError:
    print("""Argument missing!
    {}""".format(help_cmd_str))
    sys.exit(1)

pathlib.Path(os.path.dirname(hist_perc_tophits_eqtl_png)).mkdir(parents=True, exist_ok=True)
# import pdb; pdb.set_trace()
# count_df = pandas.read_csv(perc_tophits_eqtl_tsv, sep="\t")
perc_tophits_eqtl_df = pandas.read_csv(perc_tophits_eqtl_tsv, sep="\t", usecols=['gwas_id', 'loci_explained_perc'], index_col='gwas_id')
# import pdb; pdb.set_trace()

#%%
ylabel = "Percentage of GWAS"
# ylim=[0, 100]
edgecolor='k'
linewidth = 2
grid_axis = 'y'
hist_kwargs = {'density': 1, 'edgecolor': edgecolor, 'linewidth': linewidth}
stat = 'percent'
label_fontsize = 28

#%% loci_coloc_perc
title = "eQTL colocalized tophits (incl. MHC)"
data_ser = perc_tophits_eqtl_df['loci_explained_perc']
shplt = seaborn.histplot(data_ser, stat=stat, discrete=True)

plt.grid(visible=True, axis='y')
plt.title(title, fontsize=label_fontsize)
plt.xlabel("Percentage of tophits", fontsize=label_fontsize)
plt.xticks(fontsize=tick_fontsize)
plt.xticks(fontsize=tick_fontsize, rotation=0)
plt.ylabel(ylabel, fontsize=label_fontsize)
plt.yticks(fontsize=tick_fontsize)

plt.tight_layout()
plt.savefig(hist_perc_tophits_eqtl_png, dpi=dpi)
plt.close()

# #%% loci_nomhc_coloc_perc
# title = "eQTL colocalized tophits (excl. MHC)"
# data_ser = perc_tophits_eqtl_df['loci_nomhc_coloc_perc']
# shplt = seaborn.histplot(data_ser, stat=stat, discrete=True)
#
# plt.grid(visible=True, axis='y')
# plt.title(title, fontsize=label_fontsize)
# plt.xlabel("Percentage of tophits", fontsize=label_fontsize)
# plt.xticks(fontsize=tick_fontsize)
# plt.xticks(fontsize=tick_fontsize, rotation=0)
# plt.ylabel(ylabel, fontsize=label_fontsize)
# plt.yticks(fontsize=tick_fontsize)
#
# plt.tight_layout()
# plt.savefig(hist_perc_tophits_nomhc_eqtl_png, dpi=dpi)
# plt.close()


