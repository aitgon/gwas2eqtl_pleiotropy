import math
import numpy
import sqlalchemy
import os
import pandas
import pathlib
import seaborn
import sys
import matplotlib.pyplot as plt
from gwas2eqtl_pleiotropy import boxenplot_with_mannwhitneyu

from statannotations.Annotator import Annotator
from gwas2eqtl_pleiotropy.constants import seaborn_theme_dic, annotator_config_dic, label_fontsize, tick_fontsize, \
    palette, boxenplot_line_kws, boxenplot_kws

#%%
plt.rcParams["figure.figsize"] = (8, 6)
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
columns = ['chrom', 'pos38', 'rsid', 'ref', 'alt', 'eqtl_refseq_transcript_id', 'eqtl_refseq_transcript_start38', 'eqtl_refseq_transcript_end38', 'eqtl_refseq_transcript_strand']
engine = sqlalchemy.create_engine(sa_url)
with engine.begin() as conn:
    h4_df = pandas.read_sql(sqlalchemy.text(sql), con=conn)
h4_df = h4_df[columns].drop_duplicates()

#%%
count_per_rsid_gwas_df = pandas.read_excel(count_per_rsid_gwas_ods_path, engine='odf', usecols=['chrom', 'pos38', 'rsid', 'ref', 'alt', 'gwas_category_count'])
count_per_rsid_gwas_df = count_per_rsid_gwas_df.drop_duplicates()
gwas_category_count_max_int = count_per_rsid_gwas_df['gwas_category_count'].max()

#%%
m_df = h4_df.merge(count_per_rsid_gwas_df, on=['chrom', 'pos38', 'rsid', 'ref', 'alt'])
m_df.loc[m_df['gwas_category_count'] >= pleio_high_cutoff, "gwas_category_count"] = pleio_high_cutoff

# distance to tss
m_df['eqtl_gene_distance'] = math.nan

m_df.loc[m_df['eqtl_refseq_transcript_strand'] == '+', 'eqtl_gene_distance'] = (m_df['pos38'] - m_df['eqtl_refseq_transcript_start38']).abs()
m_df.loc[m_df['eqtl_refseq_transcript_strand'] == '-', 'eqtl_gene_distance'] = (m_df['pos38'] - m_df['eqtl_refseq_transcript_end38']).abs()
m_df['eqtl_gene_distance_kb'] = m_df['eqtl_gene_distance'] / 1000

########################################################################################################################
#
# Will select the CLOSEST distance per variant
#
########################################################################################################################

m2_df = m_df.sort_values(['eqtl_gene_distance', 'chrom', 'pos38', 'rsid'], ascending=True)
sel_cols = ['rsid', 'eqtl_gene_distance_kb', 'eqtl_gene_distance', 'gwas_category_count']
m3_df = m2_df[sel_cols].drop_duplicates(subset=['rsid', 'gwas_category_count'])

#%%
describe_tsv_path = os.path.join(outdir_path, "describe_closest_distance.tsv")
m3_df.groupby('gwas_category_count')['eqtl_gene_distance_kb'].apply(lambda x: x.describe()).to_csv(describe_tsv_path, sep="\t")

#%%
order = [str(x) for x in range(1, max(m3_df['gwas_category_count'].unique()) + 1)]
xticklabels = order.copy()
y = "eqtl_gene_distance_kb"
x = "gwas_category_count"

#%%
pairs = [(str(1), str(i)) for i in range(2, max(m3_df['gwas_category_count'].unique()) + 1)]
m3_df[x] = m3_df[x].astype(str)

#%% boxenplot custom stats
ax = seaborn.boxenplot(data=m3_df, x=x, y=y, palette=palette, showfliers=False, line_kws=boxenplot_line_kws, saturation=1)

ylim = [0, 300]
x1_annot1 = 0.
delta_h = 0.03

y_annot1 = ylim[1] - 0.15 * ylim[1]
y_annot2 = y_annot1 - 0.15 * ylim[1]
h_annot = ylim[1] * delta_h
x2_annot1 = x1_annot1 + 1
x2_annot2 = x1_annot1 + 2
#
group1 = m3_df.where(m3_df.gwas_category_count == '1').dropna()[y]
group2 = m3_df.where(m3_df.gwas_category_count == '3').dropna()[y]
boxenplot_with_mannwhitneyu(group1, group2, x1_annot1, x2_annot2, y_annot1, h_annot)

group1 = m3_df.where(m3_df.gwas_category_count == '1').dropna()[y]
group2 = m3_df.where(m3_df.gwas_category_count == '2').dropna()[y]
boxenplot_with_mannwhitneyu(group1, group2, x1_annot1, x2_annot1, y_annot2, h_annot)

plt.title("Distance to closest gene", fontsize=label_fontsize)
plt.xlabel("Trait category count", fontsize=label_fontsize)
plt.ylabel("Distance [kbp]", fontsize=label_fontsize)
xticks_labels = [str(x) for x in (plt.xticks()[0] + 1)]
xticks_labels[-1] = '≥' + str(xticks_labels[-1])
plt.xticks(ticks=(plt.xticks()[0]), labels=xticks_labels, fontsize=tick_fontsize, rotation=0)
plt.yticks(fontsize=tick_fontsize)
plt.ylim(ylim)

plt.tight_layout()
png_path = os.path.join(outdir_path, "boxenplot_custom_closest_distance.png")
plt.savefig(png_path)
plt.close()

#%% boxenplot stats
ax = seaborn.boxenplot(data=m3_df, x=x, y=y, palette=palette, showfliers=False, line_kws=boxenplot_line_kws, saturation=1)

annotator = Annotator(ax, pairs, data=m3_df, x=x, y=y, order=order)
annotator.configure(test='Mann-Whitney', text_format='star', **annotator_config_dic)
annotator.apply_and_annotate()

plt.title("Distance to closest gene", fontsize=label_fontsize)
plt.xlabel("Trait category count", fontsize=label_fontsize)
plt.ylabel("Distance [kbp]", fontsize=label_fontsize)
xticks_labels = [str(x) for x in (plt.xticks()[0] + 1)]
xticks_labels[-1] = '≥' + str(xticks_labels[-1])
plt.xticks(ticks=(plt.xticks()[0]), labels=xticks_labels, fontsize=tick_fontsize, rotation=0)
plt.yticks(fontsize=tick_fontsize)

plt.tight_layout()
png_path = os.path.join(outdir_path, "boxenplot_auto_closest_distance.png")
plt.savefig(png_path)
plt.close()

########################################################################################################################
#
# Will select the LARGEST distance per variant
#
########################################################################################################################

m2_df = m_df.sort_values(['eqtl_gene_distance', 'chrom', 'pos38', 'rsid'], ascending=[False, True, True, True])
sel_cols = ['rsid', 'eqtl_gene_distance_kb', 'eqtl_gene_distance', 'gwas_category_count']
m3_df = m2_df[sel_cols].drop_duplicates(subset=['rsid', 'gwas_category_count'])

#%%
describe_tsv_path = os.path.join(outdir_path, "describe_furthest_distance.tsv")
m3_df.groupby('gwas_category_count')['eqtl_gene_distance_kb'].apply(lambda x: x.describe()).to_csv(describe_tsv_path, sep="\t")

#%%
order = [str(x) for x in range(1, max(m3_df['gwas_category_count'].unique()) + 1)]
xticklabels = order.copy()
y = "eqtl_gene_distance_kb"
x = "gwas_category_count"

#%%
pairs = [(str(1), str(i)) for i in range(2, max(m3_df['gwas_category_count'].unique()) + 1)]
m3_df[x] = m3_df[x].astype(str)

#%% boxenplot custom stats
ax = seaborn.boxenplot(data=m3_df, x=x, y=y, **boxenplot_kws, line_kws=boxenplot_line_kws)

ylim = [0, 600]
x1_annot1 = 0.
delta_h = 0.03

y_annot1 = ylim[1] - 0.15 * ylim[1]
y_annot2 = y_annot1 - 0.15 * ylim[1]
h_annot = ylim[1] * delta_h
x2_annot1 = x1_annot1 + 1
x2_annot2 = x1_annot1 + 2

group1 = m3_df.where(m3_df.gwas_category_count == '1').dropna()[y]
group2 = m3_df.where(m3_df.gwas_category_count == '3').dropna()[y]
boxenplot_with_mannwhitneyu(group1, group2, x1_annot1, x2_annot2, y_annot1, h_annot)

group1 = m3_df.where(m3_df.gwas_category_count == '1').dropna()[y]
group2 = m3_df.where(m3_df.gwas_category_count == '2').dropna()[y]
boxenplot_with_mannwhitneyu(group1, group2, x1_annot1, x2_annot1, y_annot2, h_annot)

plt.title("Distance to furthest gene", fontsize=label_fontsize)
plt.xlabel("Trait category count", fontsize=label_fontsize)
plt.ylabel("Distance [kbp]", fontsize=label_fontsize)
xticks_labels = [str(x) for x in (plt.xticks()[0] + 1)]
xticks_labels[-1] = '≥' + str(xticks_labels[-1])
plt.xticks(ticks=(plt.xticks()[0]), labels=xticks_labels, fontsize=tick_fontsize, rotation=0)
plt.yticks(fontsize=tick_fontsize)
plt.ylim(ylim)

plt.tight_layout()
png_path = os.path.join(outdir_path, "boxenplot_custom_furthest_distance.png")
plt.savefig(png_path)
plt.close()

#%% boxenplot stats
ax = seaborn.boxenplot(data=m3_df, x=x, y=y, **boxenplot_kws, line_kws=boxenplot_line_kws)

annotator = Annotator(ax, pairs, data=m3_df, x=x, y=y, order=order)
annotator.configure(test='Mann-Whitney', text_format='star', **annotator_config_dic)
annotator.apply_and_annotate()

plt.title("Distance to furthest gene", fontsize=label_fontsize)
plt.xlabel("Trait category count", fontsize=label_fontsize)
plt.ylabel("Distance [kbp]", fontsize=label_fontsize)
xticks_labels = [str(x) for x in (plt.xticks()[0] + 1)]
xticks_labels[-1] = '≥' + str(xticks_labels[-1])
plt.xticks(ticks=(plt.xticks()[0]), labels=xticks_labels, fontsize=tick_fontsize, rotation=0)
plt.yticks(fontsize=tick_fontsize)

plt.tight_layout()
png_path = os.path.join(outdir_path, "boxenplot_auto_furthest_distance.png")
plt.savefig(png_path)
plt.close()
