from statannotations.Annotator import Annotator
from gwas2eqtl_pleiotropy.constants import seaborn_theme_dic
from gwas2eqtl_pleiotropy.constants import label_fontsize, tick_fontsize, dpi, boxplot_kwargs, annotator_config_dic

import matplotlib.pyplot as plt
import numpy
import os
import pandas
import pathlib
import seaborn
import sys


#%%
# Plot parameters
plt.rcParams["figure.figsize"] = (8, 6)
seaborn.set_theme(**seaborn_theme_dic)

#%%
help_cmd_str = "todo"
try:
    snp_pp_h4 = float(sys.argv[1])
    max_gwas_class_count = int(sys.argv[2])
    url = sys.argv[3]
    count_per_rsid_gwas_tsv_path = sys.argv[4]
    vlnplt_png_path = sys.argv[5]
    if len(sys.argv) > 6:
        print("""Two many arguments!
        {}""".format(help_cmd_str))
        sys.exit(1)
except IndexError:
    print("""Argument missing!
    {}""".format(help_cmd_str))
    sys.exit(1)

#%%
if not os.path.isfile(count_per_rsid_gwas_tsv_path):
    print("input file does not exit")
    sys.exit(1)

outdir_path = os.path.dirname(vlnplt_png_path)
pathlib.Path(outdir_path).mkdir(parents=True, exist_ok=True)

#%%
sql = 'select * from colocpleio where snp_pp_h4>={}'.format(snp_pp_h4)
h4_df = pandas.read_sql(sql, con=url).drop_duplicates()

#%%
count_per_rsid_gwas_df = pandas.read_csv(count_per_rsid_gwas_tsv_path, sep="\t")
gwas_class_count_max_int = count_per_rsid_gwas_df['gwas_class_count'].max()

#%%
m_df = h4_df.merge(count_per_rsid_gwas_df, on=['chrom', 'pos38', 'rsid'])

# %%
sel_cols = ['rsid', 'eqtl_gene_id', 'etissue_class']  # eqtl_gene_id per variant-etissuecategory

#%%
m2_df = m_df[['chrom', 'pos38'] + sel_cols + ['eqtl_gene_symbol', 'gwas_class_count']].drop_duplicates()
m2_df.sort_values(['gwas_class_count', 'chrom', 'pos38'], inplace=True, ascending=[False, True, True])
tsv_path = os.path.join(outdir_path, 'variants2egenes.tsv')
m2_df.to_csv(tsv_path, header=True, index=False, sep='\t')

#%% set max_gwas_class_count
m_df = m_df[sel_cols + ['gwas_class_count']]
m_df.loc[m_df['gwas_class_count'] >= max_gwas_class_count, "gwas_class_count"] = max_gwas_class_count

#%% keep unique rsid-etissue_class pairs with max. gwas category
m_df.sort_values('gwas_class_count', ascending=False, inplace=True)
m_df = m_df.drop_duplicates(subset=['rsid', 'etissue_class', 'eqtl_gene_id'], keep='first')

#%%
m_df = m_df.groupby(['rsid', 'etissue_class', 'gwas_class_count']).count()
m_df = m_df.reset_index()
m_df.columns = ['rsid', 'etissue_class', 'gwas_class_count', 'egene_count']

#%%
describe_tsv_path = os.path.join(outdir_path, "describe.tsv")
m_df.groupby('gwas_class_count')['egene_count'].apply(lambda x: x.describe()).to_csv(describe_tsv_path, sep="\t")

#%%
order = [str(x) for x in range(1, max(m_df['gwas_class_count'].unique())+1)]
xticklabels = order.copy()
xticklabels[-1] = '≥{}'.format(order[-1])
title = "Genes per eQTL-tissue"
xlabel = "GWAS class count"
ylabel = "Gene count mean"
y = "egene_count"
x = "gwas_class_count"

#%%
pairs = [(str(1), str(i)) for i in range(2, max(m_df['gwas_class_count'].unique()) + 1)]
m_df[x] = m_df[x].astype(str)
# ax = seaborn.boxplot(x=x, y=y, data=m_df, order=order, **boxplot_kwargs)
ax = seaborn.barplot(x=x, y=y, data=m_df, order=order, estimator=numpy.mean, palette="rocket_r")
annotator = Annotator(ax, pairs, data=m_df, x=x, y=y, order=order)
annotator.configure(test='Mann-Whitney', text_format='star', **annotator_config_dic)
annotator.apply_and_annotate()

plt.title(title, fontsize=label_fontsize)
plt.xlabel(xlabel, fontsize=label_fontsize)
plt.xticks(fontsize=tick_fontsize, rotation=0)
plt.ylabel(ylabel, fontsize=label_fontsize)
plt.yticks(fontsize=tick_fontsize)
ax.set_xticklabels(xticklabels)

plt.tight_layout()
plt.savefig(vlnplt_png_path, dpi=dpi)
plt.close()
