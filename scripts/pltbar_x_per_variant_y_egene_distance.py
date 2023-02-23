import math

import sqlalchemy
from statannotations.Annotator import Annotator
from gwas2eqtl_pleiotropy.constants import label_fontsize, dpi, tick_fontsize, boxplot_kwargs, annotator_config_dic

import numpy
import os
import pandas
import pathlib
import seaborn
import sys
import matplotlib.pyplot as plt


#%%
plt.rcParams["figure.figsize"] = (8, 6)
from gwas2eqtl_pleiotropy.constants import seaborn_theme_dic
seaborn.set_theme(**seaborn_theme_dic)

#%%
help_cmd_str = "todo"
try:
    snp_pp_h4 = float(sys.argv[1])
    sa_url = sys.argv[2]
    count_per_rsid_gwas_ods_path = sys.argv[3]
    vlnplt_png_path = sys.argv[4]
    if len(sys.argv) > 6:
        print("""Two many arguments!
        {}""".format(help_cmd_str))
        sys.exit(1)
except IndexError:
    print("""Argument missing!
    {}""".format(help_cmd_str))
    sys.exit(1)


#%%
if not os.path.isfile(count_per_rsid_gwas_ods_path):
    print("input file does not exit")
    sys.exit(1)

outdir_path = os.path.dirname(vlnplt_png_path)
pathlib.Path(outdir_path).mkdir(parents=True, exist_ok=True)

#%%
sql = 'select distinct * from colocpleio where snp_pp_h4>={}'.format(snp_pp_h4)
# columns = ['rsid', 'eqtl_beta', 'eqtl_gene_id', 'gwas_id', 'eqtl_id']
# h4_df = pandas.read_sql(sql, con=url).drop_duplicates()
engine = sqlalchemy.create_engine(sa_url)
with engine.begin() as conn:
    h4_df = pandas.read_sql(sqlalchemy.text(sql), con=conn).drop_duplicates()

#%%
# count_per_rsid_gwas_df = pandas.read_csv(count_per_rsid_gwas_ods_path, sep="\t")
count_per_rsid_gwas_df = pandas.read_excel(count_per_rsid_gwas_ods_path, engine='odf')
gwas_category_count_max_int = count_per_rsid_gwas_df['gwas_category_count'].max()

#%%
m_df = h4_df.merge(count_per_rsid_gwas_df, on=['chrom', 'pos38', 'rsid'])

#%%
# distance to transcript middle
# m_df['eqtl_gene_position'] = m_df['eqtl_refseq_transcript_start38'] + ((m_df['eqtl_refseq_transcript_end38'] - m_df['eqtl_refseq_transcript_start38'])/2)
# m_df['eqtl_gene_distance'] = (m_df['pos38'] - m_df['eqtl_gene_position']).abs()
# # variant within transcript
# m_df.loc[(m_df['pos38'] >= m_df['eqtl_refseq_transcript_start38']) & (m_df['pos38'] <= m_df['eqtl_refseq_transcript_end38']), 'eqtl_gene_distance'] = 0

# distance to transcript middle
m_df['eqtl_gene_distance'] = math.nan
# import pdb; pdb.set_trace()
m_df.loc[m_df['eqtl_refseq_transcript_strand'] == '+', 'eqtl_gene_distance'] = (m_df['pos38'] - m_df['eqtl_refseq_transcript_start38']).abs()
m_df.loc[m_df['eqtl_refseq_transcript_strand'] == '-', 'eqtl_gene_distance'] = (m_df['pos38'] - m_df['eqtl_refseq_transcript_end38']).abs()
m_df.loc[(m_df['pos38'] >= m_df['eqtl_refseq_transcript_start38']) & (m_df['pos38'] <= m_df['eqtl_refseq_transcript_end38']), 'eqtl_gene_distance'] = 0

sel_cols = ['rsid', 'eqtl_gene_id', 'gwas_category_count', 'eqtl_gene_distance']
m_df = m_df[sel_cols].drop_duplicates()

#%%
describe_tsv_path = os.path.join(outdir_path, "describe.tsv")
m_df.groupby('gwas_category_count')['eqtl_gene_distance'].apply(lambda x: x.describe()).to_csv(describe_tsv_path, sep="\t")

#%%
order = [str(x) for x in range(1, max(m_df['gwas_category_count'].unique())+1)]
xticklabels = order.copy()
title = "Tissues per eQTL-gene "
xlabel = "GWAS category count"
ylabel = "Tissue count mean"
y = "eqtl_gene_distance"
x = "gwas_category_count"

#%%
pairs = [(str(1), str(i)) for i in range(2, max(m_df['gwas_category_count'].unique()) + 1)]

#%% Histogram
m_df[x] = m_df[x].astype(str)
# ax = seaborn.boxplot(x=x, y=y, data=m_df, order=order, **boxplot_kwargs)
# ax = seaborn.barplot(x=x, y=y, data=m_df, order=order, estimator=numpy.mean, palette="rocket_r")
# ax = seaborn.histplot(x=x, y=y, data=m_df, palette="rocket_r")
import pdb; pdb.set_trace()
# ax = (m_df.loc[m_df['gwas_category_count'] == '1', 'eqtl_gene_distance']).hist(density=True)
# (m_df.loc[m_df['gwas_category_count'] == '2', 'eqtl_gene_distance']).hist(density=True)
# (m_df.loc[m_df['gwas_category_count'] == '3', 'eqtl_gene_distance']).hist(density=True)
# (m_df.loc[m_df['gwas_category_count'] == '4', 'eqtl_gene_distance']).hist(density=True)
# ax = seaborn.histplot(data=(m_df.loc[m_df['gwas_category_count'] == '1', 'eqtl_gene_distance']), stat="proportion", color='red')
# seaborn.histplot(data=(m_df.loc[m_df['gwas_category_count'] == '2', 'eqtl_gene_distance']), stat="proportion", color='green')
# seaborn.histplot(data=(m_df.loc[m_df['gwas_category_count'] == '3', 'eqtl_gene_distance']), stat="proportion", color='blue')
# seaborn.histplot(data=(m_df.loc[m_df['gwas_category_count'] == '4', 'eqtl_gene_distance']), stat="proportion", color='gray')
# seaborn.boxplot(x=x, y=y, data=m_df, order=order, palette="rocket_r")
# annotator = Annotator(ax, pairs, data=m_df, x=x, y=y, order=order)
# annotator.configure(test='Mann-Whitney', text_format='star', **annotator_config_dic)
# annotator.apply_and_annotate()
# ax = seaborn.histplot(data=m_df, x=y, hue=x, stat='percent', common_norm=False, element="step")
ax = seaborn.histplot(data=m_df, x=y, hue=x, stat='percent', multiple="dodge", common_norm=False, common_bins=True, shrink=.8, bins=10)

# plt.title(title, fontsize=label_fontsize)
# plt.xlabel(xlabel, fontsize=label_fontsize)
# plt.xticks(fontsize=tick_fontsize, rotation=0)
# plt.ylabel(ylabel, fontsize=label_fontsize)
# plt.yticks(fontsize=tick_fontsize)
# ax.set_xticklabels(xticklabels)

plt.tight_layout()
plt.savefig(vlnplt_png_path)
plt.close()

