import seaborn
from statannot import add_stat_annotation

from eqtl2gwas_pleiotropy.constants import label_fontsize, dpi, tick_fontsize

import os
import pandas
import pathlib
import sys
import matplotlib.pyplot as plt


# Plot parameters
plt.rcParams["figure.figsize"] = (8, 6)

#%%
help_cmd_str = "todo"
try:
    h4_annot_tsv_path = sys.argv[1]
    count_per_rsid_gwas_tsv_path = sys.argv[2]
    upper_var_gwas_cat_count = int(sys.argv[3])
    vlnplt_png_path = sys.argv[4]
    if len(sys.argv) > 5:
        print("""Two many arguments!
        {}""".format(help_cmd_str))
        sys.exit(1)
except IndexError:
    print("""Argument missing!
    {}""".format(help_cmd_str))
    sys.exit(1)

#%% Input1
if not os.path.isfile(h4_annot_tsv_path):
    print("input file does not exit")
    sys.exit(1)

#%% Input2
if not os.path.isfile(count_per_rsid_gwas_tsv_path):
    print("input file does not exit")
    sys.exit(1)

outdir_path = os.path.dirname(vlnplt_png_path)
pathlib.Path(outdir_path).mkdir(parents=True, exist_ok=True)

#%%
h4_df = pandas.read_csv(h4_annot_tsv_path, sep="\t")

#%%
count_per_rsid_gwas_df = pandas.read_csv(count_per_rsid_gwas_tsv_path, sep="\t")
gwas_category_count_max_int = count_per_rsid_gwas_df['gwas_category_count'].max()

#%%
m_df = h4_df.merge(count_per_rsid_gwas_df, on=['chrom', 'pos', 'rsid'])

# %%
sel_cols = ['rsid', 'egene', 'etissue_category', 'gwas_trait']  # gwas per egene

#%% set upper_var_gwas_cat_count
m_df = m_df[sel_cols + ['gwas_category_count']]
m_df.loc[m_df['gwas_category_count'] >= upper_var_gwas_cat_count, "gwas_category_count"] = upper_var_gwas_cat_count

#%%
m_df = m_df.drop_duplicates()
m_df = m_df.groupby(['rsid', 'egene', 'etissue_category', 'gwas_category_count']).count()  # count gwas traits
m_df = m_df.reset_index()
m_df.columns = ['rsid', 'egene', 'etissue_category', 'gwas_category_count', 'gwas_trait_count']

#%%
describe_tsv_path = os.path.join(outdir_path, "describe.tsv")
m_df.groupby('gwas_category_count')['gwas_trait_count'].apply(lambda x: x.describe()).to_csv(describe_tsv_path, sep="\t")

#%%
order = [*range(1, upper_var_gwas_cat_count+1)]
seaborn.set_theme(style="whitegrid")
xticklabels = order.copy()
xticklabels[-1] = '≥{}'.format(order[-1])
title = "GWAS per var.-gene-tissue"
xlabel = "GWAS category count"
ylabel = "GWAS trait count"
y = "gwas_trait_count"

#%%
box_pairs = [(1, i) for i in range(2, upper_var_gwas_cat_count+1)]
ax = seaborn.violinplot(x="gwas_category_count", y=y, data=m_df, order=order, palette="rocket_r", alpha=0.5)
test_results = add_stat_annotation(ax, data=m_df, x="gwas_category_count", y=y, order=order,
                                   box_pairs=box_pairs,
                                   test='Mann-Whitney', text_format='star',
                                   loc='inside', verbose=2)

plt.title(title, fontsize=label_fontsize)
plt.xlabel(xlabel, fontsize=label_fontsize)
plt.xticks(fontsize=tick_fontsize, rotation=0)
plt.ylabel(ylabel, fontsize=label_fontsize)
plt.yticks(fontsize=tick_fontsize)
ax.set_xticklabels(xticklabels)

plt.tight_layout()
plt.savefig(vlnplt_png_path)
plt.close()
