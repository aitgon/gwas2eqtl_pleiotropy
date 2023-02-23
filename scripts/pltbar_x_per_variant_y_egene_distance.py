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
sql = 'select distinct * from colocpleio2 where snp_pp_h4>={}'.format(snp_pp_h4)
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
# m_df['eqtl_gene_position'] = m_df['refseq_transcript_start38'] + ((m_df['refseq_transcript_end38'] - m_df['refseq_transcript_start38'])/2)
# m_df['eqtl_gene_distance'] = (m_df['pos38'] - m_df['eqtl_gene_position']).abs()
# # variant within transcript
# m_df.loc[(m_df['pos38'] >= m_df['refseq_transcript_start38']) & (m_df['pos38'] <= m_df['refseq_transcript_end38']), 'eqtl_gene_distance'] = 0

# distance to transcript middle
m_df['eqtl_gene_distance'] = math.nan
# import pdb; pdb.set_trace()
m_df.loc[m_df['refseq_transcript_strand'] == '+', 'eqtl_gene_distance'] = (m_df['pos38'] - m_df['refseq_transcript_start38']).abs()
m_df.loc[m_df['refseq_transcript_strand'] == '-', 'eqtl_gene_distance'] = (m_df['pos38'] - m_df['refseq_transcript_end38']).abs()
m_df.loc[(m_df['pos38'] >= m_df['refseq_transcript_start38']) & (m_df['pos38'] <= m_df['refseq_transcript_end38']), 'eqtl_gene_distance'] = 0

# import pdb; pdb.set_trace()
# import pdb; pdb.set_trace()
sel_cols = ['rsid', 'eqtl_gene_id', 'gwas_category_count', 'eqtl_gene_distance']
m_df = m_df[sel_cols].drop_duplicates()

# import pdb; pdb.set_trace()

# # %%
# sel_cols = ['rsid', 'eqtl_gene_id', 'etissue_category_term']  # tissue per variant-eqtl_gene_id
#
# #%% set max_gwas_category_count
# m_df = m_df[sel_cols + ['gwas_category_count']]
# # m_df.loc[m_df['gwas_category_count'] >= max_gwas_category_count, "gwas_category_count"] = max_gwas_category_count
#
# #%% keep unique rsid-etissue_category_term pairs with max. gwas category
# m_df.sort_values('gwas_category_count', ascending=False, inplace=True)
# m_df = m_df.drop_duplicates(subset=['rsid', 'etissue_category_term', 'eqtl_gene_id'], keep='first')

#%%
# m_df = m_df.groupby(['rsid', 'eqtl_gene_id', 'gwas_category_count']).count()
# m_df = m_df.reset_index()
# m_df.columns = ['rsid', 'eqtl_gene_id', 'gwas_category_count', 'etissue_category_term_count']

#%%
describe_tsv_path = os.path.join(outdir_path, "describe.tsv")
m_df.groupby('gwas_category_count')['eqtl_gene_distance'].apply(lambda x: x.describe()).to_csv(describe_tsv_path, sep="\t")

#%%
order = [str(x) for x in range(1, max(m_df['gwas_category_count'].unique())+1)]
xticklabels = order.copy()
# xticklabels[-1] = '≥{}'.format(order[-1])
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
ax = seaborn.histplot(x=x, y=y, data=m_df, order=order, palette="rocket_r")
# seaborn.boxplot(x=x, y=y, data=m_df, order=order, palette="rocket_r")
# annotator = Annotator(ax, pairs, data=m_df, x=x, y=y, order=order)
# annotator.configure(test='Mann-Whitney', text_format='star', **annotator_config_dic)
# annotator.apply_and_annotate()

plt.title(title, fontsize=label_fontsize)
plt.xlabel(xlabel, fontsize=label_fontsize)
plt.xticks(fontsize=tick_fontsize, rotation=0)
plt.ylabel(ylabel, fontsize=label_fontsize)
plt.yticks(fontsize=tick_fontsize)
ax.set_xticklabels(xticklabels)

plt.tight_layout()
plt.savefig(vlnplt_png_path)
plt.close()

