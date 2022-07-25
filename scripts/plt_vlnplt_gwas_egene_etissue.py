from statannot import add_stat_annotation
from eqtl2gwas_pleiotropy.constants import label_fontsize, tick_fontsize, dpi

import sys
import matplotlib.pyplot as plt
import os
import pandas
import pathlib
import seaborn as sns


plt.rcParams["figure.figsize"] = (8, 6)

#%%
help_cmd_str = "todo"
try:
    gwas_count_tsv_path = sys.argv[1]
    egene_count_tsv_path = sys.argv[2]
    etissue_count_tsv_path = sys.argv[3]
    upper_var_gwas_cat_count = int(sys.argv[4])
    boxplot_gwas_egene_png_path = sys.argv[5]
    boxplot_gwas_etissue_png_path = sys.argv[6]
    if len(sys.argv) > 7:
        print("""Two many arguments!
        {}""".format(help_cmd_str))
        sys.exit(1)
except IndexError:
    print("""Argument missing!
    {}""".format(help_cmd_str))
    sys.exit(1)

pathlib.Path(os.path.dirname(boxplot_gwas_egene_png_path)).mkdir(parents=True, exist_ok=True)
pathlib.Path(os.path.dirname(boxplot_gwas_etissue_png_path)).mkdir(parents=True, exist_ok=True)

#%%
gwas_count_df = pandas.read_csv(gwas_count_tsv_path, sep="\t", header=0)
egene_count_df = pandas.read_csv(egene_count_tsv_path, sep="\t", header=0)
etissue_count_df = pandas.read_csv(etissue_count_tsv_path, sep="\t", header=0)

#%%
merged_df = gwas_count_df.merge(egene_count_df, on=['chrom', 'pos', 'rsid'])
merged_df = merged_df.merge(etissue_count_df, on=['chrom', 'pos', 'rsid'])

#%% prepare comparisons
merged_df.loc[merged_df['gwas_category_count'] >= upper_var_gwas_cat_count, 'gwas_category_count'] = upper_var_gwas_cat_count
order = [*range(1, upper_var_gwas_cat_count+1)]
box_pairs = [(1, i) for i in range(2, upper_var_gwas_cat_count+1)]

#%%
sns.set_theme(style="whitegrid")
xticklabels = order.copy()
xticklabels[-1] = '≥{}'.format(order[-1])

#%%
ax = sns.violinplot(x="gwas_category_count", y="egene_count", data=merged_df, order=order, palette="rocket_r")
test_results = add_stat_annotation(ax, data=merged_df, x="gwas_category_count", y="egene_count", order=order,
                                   box_pairs=box_pairs,
                                   test='Mann-Whitney', text_format='star',
                                   loc='inside', verbose=2)

plt.title("Colocalized variants", fontsize=label_fontsize)
plt.xlabel("GWAS category count", fontsize=label_fontsize)
plt.xticks(fontsize=tick_fontsize, rotation=0)
plt.ylabel("eGene count", fontsize=label_fontsize)
plt.yticks(fontsize=tick_fontsize)
ax.set_xticklabels(xticklabels)

plt.tight_layout()
fig = ax.get_figure()
fig.savefig(boxplot_gwas_egene_png_path)
plt.clf()
plt.close()

#%%
ax = sns.violinplot(x="gwas_category_count", y="etissue_label_count", data=merged_df, order=order, palette="rocket_r")
test_results = add_stat_annotation(ax, data=merged_df, x="gwas_category_count", y="egene_count", order=order,
                                   box_pairs=box_pairs,
                                   test='Mann-Whitney', text_format='star',
                                   loc='inside', verbose=2)

plt.title("Colocalized variants", fontsize=label_fontsize)
plt.xlabel("GWAS category count", fontsize=label_fontsize)
plt.xticks(fontsize=tick_fontsize, rotation=0)
plt.ylabel("eTissue count", fontsize=label_fontsize)
plt.yticks(fontsize=tick_fontsize)
ax.set_xticklabels(xticklabels)

plt.tight_layout()
fig = ax.get_figure()
plt.savefig(boxplot_gwas_etissue_png_path, dpi=dpi)
plt.clf()
plt.close()
