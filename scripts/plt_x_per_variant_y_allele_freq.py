import scipy
import sqlalchemy

from gwas2eqtl_pleiotropy import boxenplot_with_mannwhitneyu
from gwas2eqtl_pleiotropy.constants import label_fontsize, tick_fontsize, boxplot_kwargs, annotator_config_dic, \
    boxenplot_kws
from statannotations.Annotator import Annotator

import matplotlib.pyplot as plt
import numpy
import os
import pandas
import pathlib
import seaborn
import sys


# Plot parameters
plt.rcParams["figure.figsize"] = (8, 6)
from gwas2eqtl_pleiotropy.constants import seaborn_theme_dic
seaborn.set_theme(**seaborn_theme_dic)

#%%
help_cmd_str = "todo"
try:
    snp_pp_h4 = float(sys.argv[1])
    pleio_high_cutoff = int(sys.argv[2])
    url = sys.argv[3]
    count_per_rsid_gwas_ods_path = sys.argv[4]
    eur_af_png_path = sys.argv[5]
    if len(sys.argv) > 6:
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

outdir_path = os.path.dirname(eur_af_png_path)
pathlib.Path(outdir_path).mkdir(parents=True, exist_ok=True)

#%%
sql = 'select * from colocpleio where snp_pp_h4>={}'.format(snp_pp_h4)
engine = sqlalchemy.create_engine(url)
with engine.begin() as conn:
    h4_df = pandas.read_sql(sqlalchemy.text(sql), con=conn).drop_duplicates()

#%%
count_per_rsid_gwas_df = pandas.read_excel(count_per_rsid_gwas_ods_path, engine='odf')
# max_gwas_category_count = count_per_rsid_gwas_df['gwas_category_count'].max()

#%%
m_df = h4_df.merge(count_per_rsid_gwas_df[['chrom', 'pos38', 'rsid', 'alt', 'gwas_category_count']], on=['chrom', 'pos38', 'rsid', 'alt'])
m_df = m_df[['chrom', 'pos38', 'rsid', 'ref', 'alt', 'afr_af', 'amr_af', 'eas_af', 'eur_af', 'sas_af', 'gwas_category_count']].drop_duplicates()
m_df.loc[m_df['gwas_category_count'] >= pleio_high_cutoff, "gwas_category_count"] = pleio_high_cutoff

#%%
order = [str(x) for x in range(1, max(m_df['gwas_category_count'].unique()) + 1)]
xticklabels = order.copy()
box_pairs = [(1, i) for i in range(1, pleio_high_cutoff + 1)]
x = 'gwas_category_count'
xlabel = "Trait category count"

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# boxplot and mann-whitney
ylabel = "Alternative allele freq."

#%%
pairs = [(str(1), str(i)) for i in range(2, pleio_high_cutoff + 1)]
m_df[x] = m_df[x].astype(str)

y_labels = ['afr_af', 'amr_af', 'eas_af', 'eur_af', 'sas_af']
y_titles = ['African population', 'American population', 'East Asian population', 'European population', 'South Asian population']

for y, ytitle in zip(y_labels, y_titles):

    # %% boxenplot
    ax = seaborn.boxenplot(x=x, y=y, data=m_df, order=order, **boxenplot_kws)

    annotator = Annotator(ax, pairs, data=m_df, x=x, y=y, order=order)
    annotator.configure(test='Mann-Whitney', text_format='star', **annotator_config_dic)
    annotator.apply_and_annotate()

    ax.set_xticklabels(xticklabels)
    plt.title(ytitle, fontsize=label_fontsize)
    plt.xlabel(xlabel, fontsize=label_fontsize)
    xticks_labels = [str(x) for x in (plt.xticks()[0] + 1)]
    xticks_labels[-1] = '≥' + str(xticks_labels[-1])
    plt.xticks(ticks=(plt.xticks()[0]), labels=xticks_labels, fontsize=tick_fontsize, rotation=0)
    plt.ylabel(ylabel, fontsize=label_fontsize)
    plt.yticks(fontsize=tick_fontsize)

    plt.tight_layout()
    this_af_png_path = eur_af_png_path.replace('eur_af', y)
    plt.savefig(this_af_png_path)
    plt.close()

    # %% boxenplot ms
    ax = seaborn.boxenplot(x=x, y=y, data=m_df, order=order, **boxenplot_kws)

    ylim = [0, 1.2]
    x1_annot1 = 0.

    delta_h = 0.03
    y_annot1 = ylim[1] - 0.15 * ylim[1]
    y_annot2 = y_annot1 - 0.15 * ylim[1]
    h_annot = ylim[1] * delta_h
    x2_annot1 = x1_annot1 + 1
    x2_annot2 = x1_annot1 + 2

    group1 = m_df.where(m_df.gwas_category_count == '1').dropna()[y]
    group2 = m_df.where(m_df.gwas_category_count == '3').dropna()[y]
    boxenplot_with_mannwhitneyu(group1, group2, x1_annot1, x2_annot2, y_annot1, h_annot)

    group1 = m_df.where(m_df.gwas_category_count == '1').dropna()[y]
    group2 = m_df.where(m_df.gwas_category_count == '2').dropna()[y]
    boxenplot_with_mannwhitneyu(group1, group2, x1_annot1, x2_annot1, y_annot2, h_annot)

    ax.set_xticklabels(xticklabels)
    plt.title(ytitle, fontsize=label_fontsize)
    plt.xlabel(xlabel, fontsize=label_fontsize)
    xticks_labels = [str(x) for x in (plt.xticks()[0] + 1)]
    xticks_labels[-1] = '≥' + str(xticks_labels[-1])
    plt.xticks(ticks=(plt.xticks()[0]), labels=xticks_labels, fontsize=tick_fontsize, rotation=0)
    plt.ylabel(ylabel, fontsize=label_fontsize)
    plt.yticks(fontsize=tick_fontsize)
    plt.ylim(ylim)

    plt.tight_layout()
    this_af_png_path = eur_af_png_path.replace('eur_af', y + "_custom")
    plt.savefig(this_af_png_path)
    plt.close()
