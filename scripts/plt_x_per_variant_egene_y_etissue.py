import sqlalchemy
from statannotations.Annotator import Annotator
from gwas2eqtl_pleiotropy.constants import label_fontsize, dpi, tick_fontsize, boxplot_kwargs, annotator_config_dic
from matplotlib.ticker import MaxNLocator

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
    pleio_high_cutoff = int(sys.argv[2])
    sa_url = sys.argv[3]
    count_per_rsid_gwas_ods_path = sys.argv[4]
    png_path = sys.argv[5]
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

outdir_path = os.path.dirname(png_path)
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

# %%
sel_cols = ['rsid', 'eqtl_gene_id', 'etissue_category_term']  # tissue per variant-eqtl_gene_id

#%% set max_gwas_category_count
m_df = m_df[sel_cols + ['gwas_category_count']]
m_df.loc[m_df['gwas_category_count'] >= pleio_high_cutoff, "gwas_category_count"] = pleio_high_cutoff

#%% keep unique rsid-etissue_category_term pairs with max. gwas category
m_df.sort_values('gwas_category_count', ascending=False, inplace=True)
m_df = m_df.drop_duplicates(subset=['rsid', 'etissue_category_term', 'eqtl_gene_id'], keep='first')

#%%
m_df = m_df.groupby(['rsid', 'eqtl_gene_id', 'gwas_category_count']).count()
m_df = m_df.reset_index()
m_df.columns = ['rsid', 'eqtl_gene_id', 'gwas_category_count', 'etissue_category_term_count']

#%%
describe_tsv_path = os.path.join(outdir_path, "describe.tsv")
m_df.groupby('gwas_category_count')['etissue_category_term_count'].apply(lambda x: x.describe()).to_csv(describe_tsv_path, sep="\t")

#%%
m_df.loc[m_df['gwas_category_count'] >= pleio_high_cutoff, 'gwas_category_count'] = '≥' + str(pleio_high_cutoff)

# order = [str(x) for x in range(1, max(m_df['gwas_category_count'].unique())+1)]
order = [str(x) for x in [*range(1, pleio_high_cutoff)] + ['≥' + str(pleio_high_cutoff)]]

# xticklabels = order.copy()
# xticklabels[-1] = '≥{}'.format(order[-1])
title = "Tissues per eQTL-gene "
xlabel = "Trait category count"
ylabel = "Tissue count mean"
y = "etissue_category_term_count"
x = "gwas_category_count"

#%%
pairs = [(str(1), str(i)) for i in order[1:]]
m_df[x] = m_df[x].astype(str)

#%% histplot
ax = seaborn.histplot(data=m_df, hue=x, x=y, hue_order=order, stat="density", cumulative=True, common_norm=False, fill=False, element="step", palette="icefire", lw=3)

ax.xaxis.set_major_locator(MaxNLocator(integer=True))
plt.title(title, fontsize=label_fontsize)
plt.xlabel("Tissue count", fontsize=label_fontsize)
plt.xticks(fontsize=tick_fontsize, rotation=0)
plt.ylabel('Cumulative density', fontsize=label_fontsize)
plt.yticks(fontsize=tick_fontsize)
plt.xlim([1, m_df['etissue_category_term_count'].max()])
plt.xlim([1, 20])

# for legend text
plt.setp(ax.get_legend().get_texts(), fontsize=tick_fontsize)
# for legend title
plt.setp(ax.get_legend().get_title(), fontsize=tick_fontsize)

frame = ax.get_legend().get_frame()
frame.set_facecolor('white')

plt.tight_layout()
plt.savefig(png_path)
plt.close()
