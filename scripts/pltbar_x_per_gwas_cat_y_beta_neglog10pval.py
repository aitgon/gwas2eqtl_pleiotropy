import pdb
import sys

from gwas2eqtl_pleiotropy.constants import label_fontsize, tick_fontsize, annotator_config_dic
from statannotations.Annotator import Annotator

import matplotlib.pyplot as plt
import numpy
import os
import pandas
import pathlib
import seaborn
import sqlalchemy


# Plot parameters
plt.rcParams["figure.figsize"] = (8, 6)
from gwas2eqtl_pleiotropy.constants import seaborn_theme_dic
seaborn.set_theme(**seaborn_theme_dic)

#%%
help_cmd_str = "todo"
try:
    snp_pp_h4 = float(sys.argv[1])
    url = sys.argv[2]
    count_per_rsid_gwas_ods_path = sys.argv[3]
    eqtl_beta_png_path = sys.argv[4]
    gwas_beta_png_path = sys.argv[5]
    eqtl_neglogpval_png_path = sys.argv[6]
    gwas_neglogpval_png_path = sys.argv[7]
    if len(sys.argv) > 8:
        print("""Two many arguments!
        {}""".format(help_cmd_str))
        sys.exit(1)
except IndexError:
    print("""Argument missing!
    {}""".format(help_cmd_str))
    sys.exit(1)

#%% Input
if not os.path.isfile(count_per_rsid_gwas_ods_path):
    print("input file does not exit")
    sys.exit(1)

outdir_path = os.path.dirname(eqtl_beta_png_path)
pathlib.Path(outdir_path).mkdir(parents=True, exist_ok=True)

#%%
sql = 'select * from colocpleio where snp_pp_h4>={}'.format(snp_pp_h4)
engine = sqlalchemy.create_engine(url)
with engine.begin() as conn:
    h4_df = pandas.read_sql(sqlalchemy.text(sql), con=conn).drop_duplicates()

#%%
count_per_rsid_gwas_df = pandas.read_excel(count_per_rsid_gwas_ods_path, engine='odf')
max_gwas_category_count = count_per_rsid_gwas_df['gwas_category_count'].max()

#%%
m_df = h4_df.merge(count_per_rsid_gwas_df, on=['chrom', 'pos38', 'rsid','alt'])

#%%
m_gwas_df = m_df[['chrom', 'pos38', 'rsid', 'ref', 'alt', 'gwas_beta', 'gwas_pval', 'gwas_id', 'gwas_category_count']].drop_duplicates()
m_eqtl_df = m_df[['chrom', 'pos38', 'rsid', 'ref', 'alt', 'eqtl_beta', 'eqtl_pval', 'eqtl_gene_id', 'eqtl_id', 'gwas_category_count']].drop_duplicates()

#%%
m_gwas_df = m_gwas_df.loc[m_gwas_df['gwas_pval']!=0]  # remove pval with zeros
m_gwas_df[['gwas_beta_abs']] = m_gwas_df[['gwas_beta']].abs()
m_gwas_df[['gwas_neglog10pval']] = -numpy.log10(m_gwas_df[['gwas_pval']])
m_gwas_df = m_gwas_df[['gwas_category_count', 'gwas_beta_abs', 'gwas_neglog10pval']]

m_eqtl_df[['eqtl_beta_abs']] = m_eqtl_df[['eqtl_beta']].abs()
m_eqtl_df[['eqtl_neglog10pval']] = -numpy.log10(m_eqtl_df[['eqtl_pval']])
m_eqtl_df = m_eqtl_df[['gwas_category_count', 'eqtl_beta_abs', 'eqtl_neglog10pval']]

#%%
order = [str(x) for x in range(1, max(m_gwas_df['gwas_category_count'].unique())+1)]
xticklabels = order.copy()
# xticklabels[-1] = '≥{}'.format(order[-1])
box_pairs = [(1, i) for i in range(1, max_gwas_category_count+1) ]
x = 'gwas_category_count'
xlabel = "GWAS category count"

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# boxplot and mann-whitney
y = "gwas_neglog10pval"
title = "GWAS neg. log10 p-val"
ylabel = "Neg. log10 p-val mean"

#%%
pairs = [(str(1), str(i)) for i in range(2, max_gwas_category_count+1)]
m_gwas_df[x] = m_gwas_df[x].astype(str)
ax = seaborn.barplot(x=x, y=y, data=m_gwas_df, order=order, estimator=numpy.mean, palette="rocket_r")

annotator = Annotator(ax, pairs, data=m_gwas_df, x=x, y=y, order=order)
annotator.configure(test='Mann-Whitney', text_format='star', **annotator_config_dic)
annotator.apply_and_annotate()

plt.title(title, fontsize=label_fontsize)
plt.xlabel(xlabel, fontsize=label_fontsize)
plt.xticks(fontsize=tick_fontsize, rotation=0)
plt.ylabel(ylabel, fontsize=label_fontsize)
plt.yticks(fontsize=tick_fontsize)
ax.set_xticklabels(xticklabels)

plt.tight_layout()
plt.savefig(gwas_neglogpval_png_path)
plt.close()

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# boxplot and mann-whitney
y = "eqtl_neglog10pval"
title = "eQTL neg. log10 p-val"
ylabel = "Neg. log10 p-val mean"

#%%
pairs = [(str(1), str(i)) for i in range(2, max_gwas_category_count+1)]
m_eqtl_df[x] = m_eqtl_df[x].astype(str)

ax = seaborn.barplot(x=x, y=y, data=m_eqtl_df, order=order, estimator=numpy.mean, palette="rocket_r")

annotator = Annotator(ax, pairs, data=m_eqtl_df, x=x, y=y, order=order)
annotator.configure(test='Mann-Whitney', text_format='star', **annotator_config_dic)
annotator.apply_and_annotate()

plt.title(title, fontsize=label_fontsize)
plt.xlabel(xlabel, fontsize=label_fontsize)
plt.xticks(fontsize=tick_fontsize, rotation=0)
plt.ylabel(ylabel, fontsize=label_fontsize)
plt.yticks(fontsize=tick_fontsize)
ax.set_xticklabels(xticklabels)

plt.tight_layout()
plt.savefig(eqtl_neglogpval_png_path)
plt.close()

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# boxplot and mann-whitney
y = "eqtl_beta_abs"
title = "eQTL beta"
ylabel = "Absolute beta mean"

#%%
pairs = [(str(1), str(i)) for i in range(2, max_gwas_category_count+1)]
m_eqtl_df[x] = m_eqtl_df[x].astype(str)

ax = seaborn.barplot(x=x, y=y, data=m_eqtl_df, order=order, estimator=numpy.mean, palette="rocket_r")

annotator = Annotator(ax, pairs, data=m_eqtl_df, x=x, y=y, order=order)
annotator.configure(test='Mann-Whitney', text_format='star', **annotator_config_dic)
annotator.apply_and_annotate()

plt.title(title, fontsize=label_fontsize)
plt.xlabel(xlabel, fontsize=label_fontsize)
plt.xticks(fontsize=tick_fontsize, rotation=0)
plt.ylabel(ylabel, fontsize=label_fontsize)
plt.yticks(fontsize=tick_fontsize)
ax.set_xticklabels(xticklabels)

plt.tight_layout()
plt.savefig(eqtl_beta_png_path)
plt.close()

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# boxplot and mann-whitney
y = "gwas_beta_abs"
title = "GWAS beta"
ylabel = "Absolute beta mean"

#%%
pairs = [(str(1), str(i)) for i in range(2, max_gwas_category_count+1)]
m_gwas_df[x] = m_gwas_df[x].astype(str)

ax = seaborn.barplot(x=x, y=y, data=m_gwas_df, order=order, estimator=numpy.mean, palette="rocket_r")

annotator = Annotator(ax, pairs, data=m_gwas_df, x=x, y=y, order=order)
annotator.configure(test='Mann-Whitney', text_format='star', **annotator_config_dic)
annotator.apply_and_annotate()

plt.title(title, fontsize=label_fontsize)
plt.xlabel(xlabel, fontsize=label_fontsize)
plt.xticks(fontsize=tick_fontsize, rotation=0)
plt.ylabel(ylabel, fontsize=label_fontsize)
plt.yticks(fontsize=tick_fontsize)
ax.set_xticklabels(xticklabels)

plt.tight_layout()
plt.savefig(gwas_beta_png_path)
plt.close()
