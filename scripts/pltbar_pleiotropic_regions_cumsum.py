import sys

import seaborn

from eqtl2gwas_pleiotropy.PathManager import PathManager
from eqtl2gwas_pleiotropy.constants import region_bin, label_fontsize, \
    tick_fontsize, dpi
from matplotlib import pyplot as plt

import math
import numpy
import os
import pandas
import pathlib


# #%% Outdir
# if not '__file__' in locals():
#     __file__ = "cmpt_pleiotropic_regions.py"
# outdir_path = os.path.join(PathManager.get_outdir_path(), os.path.basename(__file__))
# pathlib.Path(outdir_path).mkdir(parents=True, exist_ok=True)

from eqtl2gwas_pleiotropy.constants import seaborn_theme_dic
seaborn.set_theme(**seaborn_theme_dic)

#%%
help_cmd_str = "todo"
try:
    pleio_regions_tsv_path = sys.argv[1]
    png_path = sys.argv[2]
    if len(sys.argv) > 3:
        print("""Two many arguments!
        {}""".format(help_cmd_str))
        sys.exit(1)
except IndexError:
    print("""Argument missing!
    {}""".format(help_cmd_str))
    sys.exit(1)
outdir_path = os.path.dirname(png_path)
pathlib.Path(outdir_path).mkdir(parents=True, exist_ok=True)

#%%
# count_per_rsid_gwas_tsv_path = os.path.join(PathManager.get_outdir_path(), "cmpt_count_per_rsid.py", "count_per_rsid_gwas.tsv")
regions_df = pandas.read_csv(pleio_regions_tsv_path, sep="\t")

#%% barplot cumulated covering region
cumsum_df = regions_df.copy()
cumsum_df['cumsum'] = cumsum_df['end'] - cumsum_df['start']
cumsum_df = cumsum_df[['gwas_category_count', 'cumsum']].groupby('gwas_category_count').sum().reset_index()
cumsum_df.sort_values('gwas_category_count', ascending=False, inplace=True)
cumsum_df['cumsum'] = cumsum_df['cumsum'].cumsum()
cumsum_df['cumsum'] = cumsum_df['cumsum']/10e6

# import pdb; pdb.set_trace()
#%%
order = cumsum_df['gwas_category_count'].tolist()
# xticklabels = order.copy()
# xticklabels[-1] = '≥{}'.format(order[-1])
title = "Regions"
xlabel = "GWAS category count"
ylabel = "Cumulative sum [Mb]"
y = "cumsum"
x = "gwas_category_count"

#%%
ax = seaborn.barplot(x=x, y=y, data=cumsum_df, order=order, palette="rocket")

#%%
# pairs = [('1', x) for x in out_df['gwas_category_count'] if x != "1"]
# formatted_pvalues = out_df['signif'].tolist()[1:]
#
# annotator = Annotator(ax, pairs, data=out_df, x=x, y=y, order=order, size=label_fontsize)
# annotator.set_custom_annotations(formatted_pvalues)
# annotator.configure(**annotator_config_dic)
# annotator.annotate()

# import pdb; pdb.set_trace()
# ax.set_xticklabels(xticklabels)
plt.title(title, fontsize=label_fontsize)
plt.xlabel(xlabel, fontsize=label_fontsize)
plt.ylabel(ylabel, fontsize=label_fontsize)
plt.yticks(ax.get_yticks().tolist(), fontsize=tick_fontsize)
# import pdb; pdb.set_trace()
# ax.set_yticks(ax.get_yticks().tolist())

plt.tight_layout()
plt.savefig(png_path, dpi=dpi)
plt.close()

# import pdb; pdb.set_trace()