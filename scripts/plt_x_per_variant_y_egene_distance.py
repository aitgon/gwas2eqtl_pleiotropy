import math
import sqlalchemy
import os
import pandas
import pathlib
import seaborn
import sys
import matplotlib.pyplot as plt
from scipy.stats import kstest
from statannotations.Annotator import Annotator
from statannotations.stats.StatTest import StatTest

#%%
plt.rcParams["figure.figsize"] = (8, 6)
from gwas2eqtl_pleiotropy.constants import seaborn_theme_dic, annotator_config_dic

seaborn.set_theme(**seaborn_theme_dic)

#%%
help_cmd_str = "todo"
try:
    snp_pp_h4 = float(sys.argv[1])
    sa_url = sys.argv[2]
    count_per_rsid_gwas_ods_path = sys.argv[3]
    hist_png_path = sys.argv[4]
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

outdir_path = os.path.dirname(hist_png_path)
pathlib.Path(outdir_path).mkdir(parents=True, exist_ok=True)

#%%
sql = 'select distinct * from colocpleio where snp_pp_h4>={}'.format(snp_pp_h4)
columns = ['chrom', 'pos38', 'rsid', 'ref', 'alt', 'eqtl_refseq_transcript_id', 'eqtl_refseq_transcript_start38', 'eqtl_refseq_transcript_end38', 'eqtl_refseq_transcript_strand']
# columns = ['rsid', 'eqtl_beta', 'eqtl_gene_id', 'gwas_id', 'eqtl_id']
# h4_df = pandas.read_sql(sql, con=url).drop_duplicates()
engine = sqlalchemy.create_engine(sa_url)
with engine.begin() as conn:
    h4_df = pandas.read_sql(sqlalchemy.text(sql), con=conn)
h4_df = h4_df[columns].drop_duplicates()

#%%
# count_per_rsid_gwas_df = pandas.read_csv(count_per_rsid_gwas_ods_path, sep="\t")
count_per_rsid_gwas_df = pandas.read_excel(count_per_rsid_gwas_ods_path, engine='odf', usecols=['chrom', 'pos38', 'rsid', 'ref', 'alt', 'gwas_category_count'])
count_per_rsid_gwas_df = count_per_rsid_gwas_df.drop_duplicates()
gwas_category_count_max_int = count_per_rsid_gwas_df['gwas_category_count'].max()

#%%
m_df = h4_df.merge(count_per_rsid_gwas_df, on=['chrom', 'pos38', 'rsid', 'ref', 'alt'])

# distance to tss
m_df['eqtl_gene_distance'] = math.nan

m_df.loc[m_df['eqtl_refseq_transcript_strand'] == '+', 'eqtl_gene_distance'] = (m_df['pos38'] - m_df['eqtl_refseq_transcript_start38']).abs()
m_df.loc[m_df['eqtl_refseq_transcript_strand'] == '-', 'eqtl_gene_distance'] = (m_df['pos38'] - m_df['eqtl_refseq_transcript_end38']).abs()

########################################################################################################################
#
# Will select the MINIMAL distance per variant
#
########################################################################################################################

m_df.sort_values(['chrom', 'pos38', 'rsid', 'eqtl_gene_distance'], ascending=True, inplace=True)
sel_cols = ['rsid', 'gwas_category_count', 'eqtl_gene_distance']
m2_df = m_df[sel_cols].drop_duplicates(subset=['rsid', 'gwas_category_count'])

#%%
describe_tsv_path = os.path.join(outdir_path, "min_describe.tsv")
m2_df.groupby('gwas_category_count')['eqtl_gene_distance'].apply(lambda x: x.describe()).to_csv(describe_tsv_path, sep="\t")

#%%
order = [str(x) for x in range(1, max(m2_df['gwas_category_count'].unique())+1)]
xticklabels = order.copy()
y = "eqtl_gene_distance"
x = "gwas_category_count"

#%%
pairs = [(str(1), str(i)) for i in range(2, max(m2_df['gwas_category_count'].unique()) + 1)]
m2_df[x] = m2_df[x].astype(str)

#%% Histogram
seaborn.histplot(data=m2_df, x=y, hue=x, stat='percent', multiple="dodge", common_norm=False, common_bins=True, shrink=.8, bins=10, palette="rocket_r")
plt.tight_layout()
plt.savefig(hist_png_path)
plt.close()

#%% Cumulative distribution
seaborn.histplot(data=m2_df, x=y, hue=x, stat='percent', common_norm=False, element="step", fill=False, cumulative=True, palette="rocket_r")
plt.tight_layout()
hist_png_path = os.path.join(outdir_path, "cumulative.png")
plt.savefig(hist_png_path)
plt.close()

#%% Violin plot
ax = seaborn.violinplot(data=m2_df, x=x, y=y, palette="rocket_r")

annotator = Annotator(ax, pairs, data=m2_df, x=x, y=y, order=order)
annotator.configure(test='Mann-Whitney', text_format='star', **annotator_config_dic)
annotator.apply_and_annotate()

plt.tight_layout()
plt.xlabel("GWAS category count")
plt.ylabel("Minimal eQTL gene distance")
hist_png_path = os.path.join(outdir_path, "violin.png")
plt.savefig(hist_png_path)
plt.close()

########################################################################################################################
#
# Will select the MAXIMAL distance per variant
#
########################################################################################################################

m_df.sort_values(['chrom', 'pos38', 'rsid', 'eqtl_gene_distance'], ascending=False, inplace=True)
sel_cols = ['rsid', 'gwas_category_count', 'eqtl_gene_distance']
m2_df = m_df[sel_cols].drop_duplicates(subset=['rsid', 'gwas_category_count'])

#%%
describe_tsv_path = os.path.join(outdir_path, "max_describe.tsv")
m2_df.groupby('gwas_category_count')['eqtl_gene_distance'].apply(lambda x: x.describe()).to_csv(describe_tsv_path, sep="\t")

#%%
order = [str(x) for x in range(1, max(m2_df['gwas_category_count'].unique())+1)]
xticklabels = order.copy()
y = "eqtl_gene_distance"
x = "gwas_category_count"

#%%
pairs = [(str(1), str(i)) for i in range(2, max(m2_df['gwas_category_count'].unique()) + 1)]
m2_df[x] = m2_df[x].astype(str)

#%% Histogram
seaborn.histplot(data=m2_df, x=y, hue=x, stat='percent', multiple="dodge", common_norm=False, common_bins=True, shrink=.8, bins=10, palette="rocket_r")
plt.tight_layout()
hist_png_path = os.path.join(outdir_path, "max_histogram.png")
plt.savefig(hist_png_path)
plt.close()

#%% Cumulative distribution
seaborn.histplot(data=m2_df, x=y, hue=x, stat='percent', common_norm=False, element="step", fill=False, cumulative=True, palette="rocket_r")
plt.tight_layout()
hist_png_path = os.path.join(outdir_path, "max_cumulative.png")
plt.savefig(hist_png_path)
plt.close()

#%% Violin plot
ax = seaborn.violinplot(data=m2_df, x=x, y=y, palette="rocket_r")

annotator = Annotator(ax, pairs, data=m2_df, x=x, y=y, order=order)
annotator.configure(test='Mann-Whitney', text_format='star', **annotator_config_dic)
annotator.apply_and_annotate()

plt.tight_layout()
plt.xlabel("GWAS category count")
plt.ylabel("Maximal eQTL gene distance")
hist_png_path = os.path.join(outdir_path, "max_violin.png")
plt.savefig(hist_png_path)
plt.close()

