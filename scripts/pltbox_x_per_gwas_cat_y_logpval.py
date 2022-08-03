import numpy

from eqtl2gwas_pleiotropy.constants import label_fontsize, tick_fontsize, alpha, boxplot_kwargs, dpi, \
    annotator_config_dic

import matplotlib.pyplot as plt
import os
import pandas
import pathlib
import seaborn
import sys

from statannotations.Annotator import Annotator
# Plot parameters
plt.rcParams["figure.figsize"] = (8, 6)
from eqtl2gwas_pleiotropy.constants import seaborn_theme_dic
seaborn.set_theme(**seaborn_theme_dic)

#%%
help_cmd_str = "todo"
try:
    h4_annot_tsv_path = sys.argv[1]
    count_per_rsid_gwas_tsv_path = sys.argv[2]
    upper_var_gwas_cat_count = int(sys.argv[3])
    eqtl_pval_png_path = sys.argv[4]
    gwas_pval_png_path = sys.argv[5]
    if len(sys.argv) > 6:
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

outdir_path = os.path.dirname(eqtl_pval_png_path)
pathlib.Path(outdir_path).mkdir(parents=True, exist_ok=True)

#%%
h4_df = pandas.read_csv(h4_annot_tsv_path, sep="\t")

#%%
count_per_rsid_gwas_df = pandas.read_csv(count_per_rsid_gwas_tsv_path, sep="\t")
gwas_category_count_max_int = count_per_rsid_gwas_df['gwas_category_count'].max()

#%%
m_df = h4_df.merge(count_per_rsid_gwas_df, on=['chrom', 'pos', 'rsid'])

#%%
m_df = m_df[['rsid', 'eqtl_beta', 'eqtl_pvalue', 'egene', 'etissue_category', 'gwas_beta', 'gwas_pvalue', 'gwas_identifier', 'gwas_category_count']].drop_duplicates()
m_df['eqtl_logpval'] = m_df['eqtl_pvalue'].apply(lambda x: -numpy.log(x))
m_df['gwas_logpval'] = m_df['gwas_pvalue'].apply(lambda x: -numpy.log(x))

#%%
m_df.loc[m_df['gwas_category_count'] >= upper_var_gwas_cat_count, "gwas_category_count"] = upper_var_gwas_cat_count
order = [str(x) for x in range(1, upper_var_gwas_cat_count+1)]
xticklabels = order.copy()
xticklabels[-1] = '≥{}'.format(order[-1])
pairs = [(str(1), str(i)) for i in range(2, upper_var_gwas_cat_count + 1)]
x = 'gwas_category_count'
xlabel = "GWAS category count"

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
y = "eqtl_logpval"
title = y
ylabel = y

#%%
plt_df = m_df[['gwas_category_count', 'rsid', 'egene', 'etissue_category', y]].drop_duplicates()
plt_df[y] = plt_df[y].abs()

#%%
describe_tsv_path = os.path.join(outdir_path, y + "_describe.tsv")
describe_df = plt_df.groupby(['gwas_category_count'])[y].apply(lambda x: x.describe()).to_csv(describe_tsv_path, sep="\t")

#%%
plt_df[x] = plt_df[x].astype(str)
ax = seaborn.boxplot(x=x, y=y, data=plt_df, order=order, **boxplot_kwargs)
annotator = Annotator(ax, pairs, data=plt_df, x=x, y=y, order=order)
annotator.configure(test='Mann-Whitney', text_format='star', **annotator_config_dic)
annotator.apply_and_annotate()

plt.title(title, fontsize=label_fontsize)
plt.xlabel(xlabel, fontsize=label_fontsize)
plt.xticks(fontsize=tick_fontsize, rotation=0)
plt.ylabel(ylabel, fontsize=label_fontsize)
plt.yticks(fontsize=tick_fontsize)
ax.set_xticklabels(xticklabels)

plt.tight_layout()
plt.savefig(eqtl_pval_png_path)
plt.close()

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
y = "gwas_logpval"
title = y
ylabel = y

#%%
plt_df = m_df[['gwas_category_count', 'rsid', 'gwas_identifier', y]].drop_duplicates()
plt_df[y] = plt_df[y].abs()

#%%
describe_tsv_path = os.path.join(outdir_path, y + "_describe.tsv")
describe_df = plt_df.groupby(['gwas_category_count'])[y].apply(lambda x: x.describe()).to_csv(describe_tsv_path, sep="\t")

#%%
plt_df[x] = plt_df[x].astype(str)
ax = seaborn.boxplot(x=x, y=y, data=plt_df, order=order, **boxplot_kwargs)
annotator = Annotator(ax, pairs, data=plt_df, x=x, y=y, order=order)
annotator.configure(test='Mann-Whitney', text_format='star', loc='inside')
annotator.apply_and_annotate()

plt.title(title, fontsize=label_fontsize)
plt.xlabel(xlabel, fontsize=label_fontsize)
plt.xticks(fontsize=tick_fontsize, rotation=0)
plt.ylabel(ylabel, fontsize=label_fontsize)
plt.yticks(fontsize=tick_fontsize)
ax.set_xticklabels(xticklabels)

plt.tight_layout()
plt.savefig(gwas_pval_png_path, dpi=dpi)
plt.close()
