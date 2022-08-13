#%%
import sys
import matplotlib.pyplot as plt
import numpy
import os
import pandas
import pathlib
import seaborn

from eqtl2gwas_pleiotropy.constants import label_fontsize, tick_fontsize, dpi

plt.rcParams["figure.figsize"] = (8, 6)
from eqtl2gwas_pleiotropy.constants import seaborn_theme_dic
seaborn.set_theme(**seaborn_theme_dic)

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

#%%
ylabel = "Percent"
title = "Colocalized eQTL/GWAS variants"
ylim=[0, 100]
edgecolor='k'
linewidth = 2
grid_axis = 'y'
hist_kwargs = {'density': 1, 'edgecolor': edgecolor, 'linewidth': linewidth}
stat = 'percent'
label_fontsize = 28

#%% gwas
count_df = pandas.read_csv(gwas_count_tsv_path, sep="\t", header=0)
# bins = numpy.array(range(6))
data_ser = count_df['gwas_category_count']
shplt = seaborn.histplot(data_ser, stat=stat, discrete=True)

plt.grid(visible=True, axis='y')
plt.title(title.format(" and GWAS categories"), fontsize=label_fontsize)
plt.xlabel("GWAS category count", fontsize=label_fontsize)
plt.xticks(fontsize=tick_fontsize)
plt.xticks(fontsize=tick_fontsize, rotation=0)
plt.ylabel(ylabel, fontsize=label_fontsize)
plt.ylim(ylim)
# plt.yscale('log')
plt.yticks(fontsize=tick_fontsize)

plt.tight_layout()
plt.savefig(hist_rsid_gwas_path, dpi=dpi)
plt.close()

with open(hist_rsid_gwas_path + ".txt", 'w') as fout:
    fout.write('\n'.join([str(h.get_height()) for h in shplt.patches]))

#%% egene
count_df = pandas.read_csv(egene_count_tsv_path, sep="\t", header=0)
data_ser = count_df['egene_count']
# plt.hist(data_ser, **hist_kwargs)
shplt = seaborn.histplot(data_ser, stat=stat, discrete=True)

plt.grid(visible=True, axis='y')
plt.title(title.format(" and egenes"), fontsize=label_fontsize)
plt.xlabel("eGene count", fontsize=label_fontsize)
plt.xticks(fontsize=tick_fontsize)
plt.xticks(fontsize=tick_fontsize, rotation=0)
plt.ylabel(ylabel, fontsize=label_fontsize)
plt.ylim(ylim)
# plt.yscale('log')
plt.yticks(fontsize=tick_fontsize)

plt.tight_layout()
# png_path = os.path.join(outdir_path, "hist_egene.png")
plt.savefig(hist_rsid_egene_path, dpi=dpi)
plt.close()

with open(hist_rsid_egene_path + ".txt", 'w') as fout:
    fout.write('\n'.join([str(h.get_height()) for h in shplt.patches]))

#%% etissue
count_df = pandas.read_csv(etissue_count_tsv_path, sep="\t", header=0)
data_ser = count_df['etissue_label_count']
# plt.hist(data_ser, **hist_kwargs)
shplt = seaborn.histplot(data_ser, stat=stat, discrete=True)

plt.grid(visible=True, axis='y')
plt.title(title.format(" and etissues"), fontsize=label_fontsize)
plt.xlabel("eTissue category count", fontsize=label_fontsize)
plt.xticks(fontsize=tick_fontsize)
plt.xticks(fontsize=tick_fontsize, rotation=0)
plt.ylabel(ylabel, fontsize=label_fontsize)
plt.ylim(ylim)
# plt.yscale('log')
plt.yticks(fontsize=tick_fontsize)

plt.tight_layout()
# png_path = os.path.join(outdir_path, "hist_etissue.png")
plt.savefig(hist_rsid_etissue_path, dpi=dpi)
plt.close()

with open(hist_rsid_etissue_path + ".txt", 'w') as fout:
    fout.write('\n'.join([str(h.get_height()) for h in shplt.patches]))
