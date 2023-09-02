from gwas2eqtl_pleiotropy.constants import label_fontsize, tick_fontsize, dpi
from matplotlib import pyplot as plt
from gwas2eqtl_pleiotropy.constants import seaborn_theme_dic

import os
import pandas
import pathlib
import seaborn
import sys

seaborn.set_theme(**seaborn_theme_dic)

#%%
help_cmd_str = "todo"
try:
    count_per_rsid_gwas_ods_path = sys.argv[1]
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
# regions_df = pandas.read_csv(pleio_regions_tsv_path, sep="\t")
regions_df = pandas.read_excel(count_per_rsid_gwas_ods_path, engine='odf')

#%% barplot cumulated covering region
cumsum_df = regions_df.copy()
cumsum_df['cumsum'] = cumsum_df['end'] - cumsum_df['start']
cumsum_df = cumsum_df[['gwas_category_count', 'cumsum']].groupby('gwas_category_count').sum().reset_index()
cumsum_df.sort_values('gwas_category_count', ascending=False, inplace=True)
cumsum_df['cumsum'] = cumsum_df['cumsum'].cumsum()
cumsum_df['cumsum'] = cumsum_df['cumsum']/10e6

#%%
order = cumsum_df['gwas_category_count'].tolist()
title = "Coverage of pleiotropic regions"
xlabel = "Trait category count"
ylabel = "Cumulative sum [Mbp]"
y = "cumsum"
x = "gwas_category_count"

#%%
ax = seaborn.barplot(x=x, y=y, data=cumsum_df, order=order, palette="rocket")

#%%
label_fontsize = 26
plt.title(title, fontsize=label_fontsize)
plt.xlabel(xlabel, fontsize=label_fontsize)
plt.ylabel(ylabel, fontsize=label_fontsize)
plt.yticks(fontsize=tick_fontsize)
plt.xticks(fontsize=tick_fontsize)
plt.yscale('log')

plt.tight_layout()
plt.savefig(png_path, dpi=dpi)
plt.close()