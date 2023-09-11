import scipy
from gwas2eqtl_pleiotropy import boxenplot_with_mannwhitneyu

from gwas2eqtl_pleiotropy.constants import label_fontsize, tick_fontsize, annotator_config_dic
from statannotations.Annotator import Annotator

import sys
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
snp_pp_h4 = 0.5
pleio_high_cutoff = 3
sa_url = "postgresql://postgres:postgres@0.0.0.0:5435/postgres"
count_per_rsid_gwas_ods_path = "out/v20230901/pval_5e-08/r2_0.1/kb_1000/window_1000000/75_50/cmpt_count_per_rsid.py/count_per_rsid_gwas_egene_etissue.ods"
eqtl_beta_png_path = "out/v20230901/pval_5e-08/r2_0.1/kb_1000/window_1000000/75_50/plt_x_per_gwas_count_y_beta_neglog10pval.py/eqtl_beta.png"
gwas_beta_png_path = "out/v20230901/pval_5e-08/r2_0.1/kb_1000/window_1000000/75_50/plt_x_per_gwas_count_y_beta_neglog10pval.py/gwas_beta.png"
eqtl_neglogpval_png_path = "out/v20230901/pval_5e-08/r2_0.1/kb_1000/window_1000000/75_50/plt_x_per_gwas_count_y_beta_neglog10pval.py/eqtl_neglog10pval.png"
gwas_neglogpval_png_path = "out/v20230901/pval_5e-08/r2_0.1/kb_1000/window_1000000/75_50/plt_x_per_gwas_count_y_beta_neglog10pval.py/gwas_neglog10pval.png"

#%%
help_cmd_str = "todo"
try:
    snp_pp_h4 = float(sys.argv[1])
    pleio_high_cutoff = int(sys.argv[2])
    sa_url = sys.argv[3]
    count_per_rsid_gwas_ods_path = sys.argv[4]
    eqtl_beta_png_path = sys.argv[5]
    gwas_beta_png_path = sys.argv[6]
    eqtl_neglogpval_png_path = sys.argv[7]
    gwas_neglogpval_png_path = sys.argv[8]
    if len(sys.argv) > 9:
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

outdir_path = os.path.dirname(gwas_beta_png_path)
pathlib.Path(outdir_path).mkdir(parents=True, exist_ok=True)

#%%
sql = 'select * from colocpleio where snp_pp_h4>={}'.format(snp_pp_h4)
engine = sqlalchemy.create_engine(sa_url)
with engine.begin() as conn:
    h4_df = pandas.read_sql(sqlalchemy.text(sql), con=conn).drop_duplicates()

#%%
count_per_rsid_gwas_df = pandas.read_excel(count_per_rsid_gwas_ods_path, engine='odf')
# max_gwas_category_count = count_per_rsid_gwas_df['gwas_category_count'].max()

#%%
m_df = h4_df.merge(count_per_rsid_gwas_df, on=['chrom', 'pos38', 'rsid', 'ref', 'alt'])

#%%
m_gwas_df = m_df[['chrom', 'pos38', 'rsid', 'ref', 'alt', 'gwas_beta', 'gwas_pval', 'gwas_id', 'gwas_category_count']].drop_duplicates()
m_eqtl_df = m_df[['chrom', 'pos38', 'rsid', 'ref', 'alt', 'eqtl_beta', 'eqtl_pval', 'eqtl_gene_id', 'eqtl_id', 'gwas_category_count']].drop_duplicates()

#%%
m_gwas_df = m_gwas_df.loc[m_gwas_df['gwas_pval']!=0]  # remove pval with zeros
m_gwas_df[['gwas_beta_abs']] = m_gwas_df[['gwas_beta']].abs()
m_gwas_df[['gwas_neglog10pval']] = -numpy.log10(m_gwas_df[['gwas_pval']])
m_gwas_df = m_gwas_df[['gwas_category_count', 'gwas_beta_abs', 'gwas_neglog10pval']]
m_gwas_df.loc[m_gwas_df['gwas_category_count'] >= pleio_high_cutoff, "gwas_category_count"] = pleio_high_cutoff

m_eqtl_df[['eqtl_beta_abs']] = m_eqtl_df[['eqtl_beta']].abs()
m_eqtl_df[['eqtl_neglog10pval']] = -numpy.log10(m_eqtl_df[['eqtl_pval']])
m_eqtl_df = m_eqtl_df[['gwas_category_count', 'eqtl_beta_abs', 'eqtl_neglog10pval']]
m_eqtl_df.loc[m_eqtl_df['gwas_category_count'] >= pleio_high_cutoff, "gwas_category_count"] = pleio_high_cutoff


#%%
order = [str(x) for x in range(1, pleio_high_cutoff + 1)]
xticklabels = order.copy()
# xticklabels[-1] = '≥{}'.format(order[-1])
box_pairs = [(1, i) for i in range(1, pleio_high_cutoff + 1)]
x = 'gwas_category_count'
xlabel = "Trait category count"

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# boxplot and mann-whitney
y = "eqtl_beta_abs"
title = "eQTL effect"
ylabel = "Absolute beta mean"

#%%
pairs = [(str(1), str(i)) for i in range(2, pleio_high_cutoff + 1)]
m_eqtl_df[x] = m_eqtl_df[x].astype(str)

ax = seaborn.barplot(x=x, y=y, data=m_eqtl_df, order=order, estimator=numpy.mean, palette="rocket_r")

annotator = Annotator(ax, pairs, data=m_eqtl_df, x=x, y=y, order=order)
annotator.configure(test='Mann-Whitney', text_format='star', **annotator_config_dic)
annotator.apply_and_annotate()

plt.title(title, fontsize=label_fontsize)
plt.xlabel(xlabel, fontsize=label_fontsize)
xticks_labels = [str(x) for x in (plt.xticks()[0] + 1)]
xticks_labels[-1] = '≥' + str(xticks_labels[-1])
plt.xticks(ticks=(plt.xticks()[0]), labels=xticks_labels, fontsize=tick_fontsize, rotation=0)
plt.ylabel(ylabel, fontsize=label_fontsize)
plt.yticks(fontsize=tick_fontsize)
plt.ylim([0.3, 0.7])

plt.tight_layout()
plt.savefig(eqtl_beta_png_path)
plt.close()

#%% boxenplot
ax = seaborn.boxenplot(x=x, y=y, data=m_eqtl_df, showfliers=False, palette="rocket_r")

annotator = Annotator(ax, pairs, data=m_eqtl_df, x=x, y=y, order=order)
annotator.configure(test='Mann-Whitney', text_format='star', **annotator_config_dic)
annotator.apply_and_annotate()

plt.title(title, fontsize=label_fontsize)
plt.xlabel(xlabel, fontsize=label_fontsize)
xticks_labels = [str(x) for x in (plt.xticks()[0] + 1)]
xticks_labels[-1] = '≥' + str(xticks_labels[-1])
plt.xticks(ticks=(plt.xticks()[0]), labels=xticks_labels, fontsize=tick_fontsize, rotation=0)
plt.ylabel(ylabel, fontsize=label_fontsize)
plt.yticks(fontsize=tick_fontsize)
plt.ylim([0, 1.5])

plt.tight_layout()
png_path = os.path.join(outdir_path, "eqtl_effect_boxenplot.png")
plt.savefig(png_path)
plt.close()

#%% boxenplot
ax = seaborn.boxenplot(x=x, y=y, data=m_eqtl_df, showfliers=False, palette="rocket_r")

group1 = m_eqtl_df.where(m_eqtl_df.gwas_category_count == '1').dropna()['eqtl_beta_abs']
group2 = m_eqtl_df.where(m_eqtl_df.gwas_category_count == '2').dropna()['eqtl_beta_abs']
x1 = 0; x2 = 1; y = 0.6; h = 0.05;
boxenplot_with_mannwhitneyu(group1, group2, x1, x2, y, h)

x1=0; x2=2; y=0.8; h=0.1
group1 = m_eqtl_df.where(m_eqtl_df.gwas_category_count == '1').dropna()['eqtl_beta_abs']
group2 = m_eqtl_df.where(m_eqtl_df.gwas_category_count == '3').dropna()['eqtl_beta_abs']
boxenplot_with_mannwhitneyu(group1, group2, x1, x2, y, h)

plt.title(title, fontsize=label_fontsize)
plt.xlabel(xlabel, fontsize=label_fontsize)
xticks_labels = [str(x) for x in (plt.xticks()[0] + 1)]
xticks_labels[-1] = '≥' + str(xticks_labels[-1])
plt.xticks(ticks=(plt.xticks()[0]), labels=xticks_labels, fontsize=tick_fontsize, rotation=0)
plt.ylabel(ylabel, fontsize=label_fontsize)
plt.yticks(fontsize=tick_fontsize)
plt.ylim([0, 1])

plt.tight_layout()
png_path = os.path.join(outdir_path, "eqtl_effect_ms_boxenplot.png")
plt.savefig(png_path)
plt.close()

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# boxplot and mann-whitney
y = "gwas_beta_abs"
title = "GWAS effect"
ylabel = "Absolute beta mean"

#%%
pairs = [(str(1), str(i)) for i in range(2, pleio_high_cutoff + 1)]
m_gwas_df[x] = m_gwas_df[x].astype(str)

ax = seaborn.barplot(x=x, y=y, data=m_gwas_df, order=order, estimator=numpy.mean, palette="rocket_r")

annotator = Annotator(ax, pairs, data=m_gwas_df, x=x, y=y, order=order)
annotator.configure(test='Mann-Whitney', text_format='star', **annotator_config_dic)
annotator.apply_and_annotate()

plt.title(title, fontsize=label_fontsize)
plt.xlabel(xlabel, fontsize=label_fontsize)
xticks_labels = [str(x) for x in (plt.xticks()[0] + 1)]
xticks_labels[-1] = '≥' + str(xticks_labels[-1])
plt.xticks(ticks=(plt.xticks()[0]), labels=xticks_labels, fontsize=tick_fontsize, rotation=0)
plt.ylabel(ylabel, fontsize=label_fontsize)
plt.yticks(fontsize=tick_fontsize)

plt.tight_layout()
plt.savefig(gwas_beta_png_path)
plt.close()

#%% boxenplot
ax = seaborn.boxenplot(x=x, y=y, data=m_gwas_df, showfliers=True, palette="rocket_r")

annotator = Annotator(ax, pairs, data=m_gwas_df, x=x, y=y, order=order, loc='inside')
annotator.configure(test='Mann-Whitney', text_format='star', **annotator_config_dic)
annotator.apply_and_annotate()

plt.title(title, fontsize=label_fontsize)
plt.xlabel(xlabel, fontsize=label_fontsize)
xticks_labels = [str(x) for x in (plt.xticks()[0] + 1)]
xticks_labels[-1] = '≥' + str(xticks_labels[-1])
plt.xticks(ticks=(plt.xticks()[0]), labels=xticks_labels, fontsize=tick_fontsize, rotation=0)
plt.ylabel(ylabel, fontsize=label_fontsize)
plt.yticks(fontsize=tick_fontsize)

plt.tight_layout()
png_path = os.path.join(outdir_path, "gwas_effect_boxenplot.png")
plt.savefig(png_path)
plt.close()

#%% boxenplot ms
ax = seaborn.boxenplot(x=x, y=y, data=m_gwas_df, showfliers=False, palette="rocket_r")

# annotator = Annotator(ax, pairs, data=m_gwas_df, x=x, y=y, order=order, loc='inside')
# annotator.configure(test='Mann-Whitney', text_format='star', **annotator_config_dic)
# annotator.apply_and_annotate()

group1 = m_gwas_df.where(m_gwas_df.gwas_category_count == '1').dropna()[y]
group2 = m_gwas_df.where(m_gwas_df.gwas_category_count == '2').dropna()[y]
x1 = 0; x2 = 1; annot_y = 0.4; h = 0.05;
boxenplot_with_mannwhitneyu(group1, group2, x1, x2, annot_y, h)

group1 = m_gwas_df.where(m_gwas_df.gwas_category_count == '1').dropna()[y]
group2 = m_gwas_df.where(m_gwas_df.gwas_category_count == '3').dropna()[y]
x1 = 0; x2 = 2; annot_y = 0.6; h = 0.05;
boxenplot_with_mannwhitneyu(group1, group2, x1, x2, annot_y, h)

plt.title(title, fontsize=label_fontsize)
plt.xlabel(xlabel, fontsize=label_fontsize)
xticks_labels = [str(x) for x in (plt.xticks()[0] + 1)]
xticks_labels[-1] = '≥' + str(xticks_labels[-1])
plt.xticks(ticks=(plt.xticks()[0]), labels=xticks_labels, fontsize=tick_fontsize, rotation=0)
plt.ylabel(ylabel, fontsize=label_fontsize)
plt.yticks(fontsize=tick_fontsize)
plt.ylim([0, 0.8])

plt.tight_layout()
png_path = os.path.join(outdir_path, "gwas_effect_ms_boxenplot.png")
plt.savefig(png_path)
plt.close()

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# violin and mann-whitney
y = "gwas_beta_abs"
title = "GWAS effect"
ylabel = "Absolute beta"

#%%
pairs = [(str(1), str(i)) for i in range(2, pleio_high_cutoff + 1)]
m_gwas_df[x] = m_gwas_df[x].astype(str)

ax = seaborn.violinplot(x=x, y=y, data=m_gwas_df, order=order, estimator=numpy.mean, palette="rocket_r")

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
gwas_beta_png_path = os.path.join(outdir_path, "gwas_beta_violin.png")
plt.savefig(gwas_beta_png_path)
plt.close()


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# boxplot and mann-whitney
y = "eqtl_neglog10pval"
title = "eQTL significance"
ylabel = "Neg. log10 p-val mean"

#%%
pairs = [(str(1), str(i)) for i in range(2, pleio_high_cutoff + 1)]
m_eqtl_df[x] = m_eqtl_df[x].astype(str)

ax = seaborn.barplot(x=x, y=y, data=m_eqtl_df, order=order, estimator=numpy.mean, palette="rocket_r")
annotator = Annotator(ax, pairs, data=m_eqtl_df, x=x, y=y, order=order)
annotator.configure(test='Mann-Whitney', text_format='star', **annotator_config_dic)
annotator.apply_and_annotate()

plt.title(title, fontsize=label_fontsize)
plt.xlabel(xlabel, fontsize=label_fontsize)
xticks_labels = [str(x) for x in (plt.xticks()[0] + 1)]
xticks_labels[-1] = '≥' + str(xticks_labels[-1])
plt.xticks(ticks=(plt.xticks()[0]), labels=xticks_labels, fontsize=tick_fontsize, rotation=0)
plt.ylabel(ylabel, fontsize=label_fontsize)
plt.yticks(fontsize=tick_fontsize)

plt.tight_layout()
plt.savefig(eqtl_neglogpval_png_path)
plt.close()

#%% boxenplot ms
ax = seaborn.boxenplot(x=x, y=y, data=m_eqtl_df, showfliers=True, palette="rocket_r")

group1 = m_eqtl_df.where(m_eqtl_df.gwas_category_count == '1').dropna()[y]
group2 = m_eqtl_df.where(m_eqtl_df.gwas_category_count == '2').dropna()[y]
x1 = 0; x2 = 1; annot_y = 20; h = 1;
boxenplot_with_mannwhitneyu(group1, group2, x1, x2, annot_y, h)

group1 = m_eqtl_df.where(m_eqtl_df.gwas_category_count == '1').dropna()[y]
group2 = m_eqtl_df.where(m_eqtl_df.gwas_category_count == '3').dropna()[y]
x1 = 0; x2 = 2; annot_y = 25; h = 1;
boxenplot_with_mannwhitneyu(group1, group2, x1, x2, annot_y, h)

plt.title(title, fontsize=label_fontsize)
plt.xlabel(xlabel, fontsize=label_fontsize)
xticks_labels = [str(x) for x in (plt.xticks()[0] + 1)]
xticks_labels[-1] = '≥' + str(xticks_labels[-1])
plt.xticks(ticks=(plt.xticks()[0]), labels=xticks_labels, fontsize=tick_fontsize, rotation=0)
plt.ylabel(ylabel, fontsize=label_fontsize)
plt.yticks(fontsize=tick_fontsize)
plt.ylim([0, 30])

plt.tight_layout()
png_path = os.path.join(outdir_path, "eqtl_signif_ms_boxenplot.png")
plt.savefig(png_path)
plt.close()

#%% boxenplot
ax = seaborn.boxenplot(x=x, y=y, data=m_eqtl_df, showfliers=True, palette="rocket_r")

annotator = Annotator(ax, pairs, data=m_eqtl_df, x=x, y=y, order=order)
annotator.configure(test='Mann-Whitney', text_format='star', loc='inside', **annotator_config_dic)
annotator.apply_and_annotate()

plt.title(title, fontsize=label_fontsize)
plt.xlabel(xlabel, fontsize=label_fontsize)
xticks_labels = [str(x) for x in (plt.xticks()[0] + 1)]
xticks_labels[-1] = '≥' + str(xticks_labels[-1])
plt.xticks(ticks=(plt.xticks()[0]), labels=xticks_labels, fontsize=tick_fontsize, rotation=0)
plt.ylabel(ylabel, fontsize=label_fontsize)
plt.yticks(fontsize=tick_fontsize)

plt.tight_layout()
png_path = os.path.join(outdir_path, "eqtl_signif_boxenplot.png")
plt.savefig(png_path)
plt.close()

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# boxplot and mann-whitney
y = "gwas_neglog10pval"
title = "GWAS significance"
ylabel = "Neg. log10 p-val mean"

#%% barplot
pairs = [(str(1), str(i)) for i in range(2, pleio_high_cutoff + 1)]
m_gwas_df[x] = m_gwas_df[x].astype(str)
ax = seaborn.barplot(x=x, y=y, data=m_gwas_df, order=order, estimator=numpy.mean, palette="rocket_r")

annotator = Annotator(ax, pairs, data=m_gwas_df, x=x, y=y, order=order)
annotator.configure(test='Mann-Whitney', text_format='star', **annotator_config_dic)
annotator.apply_and_annotate()

plt.title(title, fontsize=label_fontsize)
plt.xlabel(xlabel, fontsize=label_fontsize)
xticks_labels = [str(x) for x in (plt.xticks()[0] + 1)]
xticks_labels[-1] = '≥' + str(xticks_labels[-1])
plt.xticks(ticks=(plt.xticks()[0]), labels=xticks_labels, fontsize=tick_fontsize, rotation=0)
plt.ylabel(ylabel, fontsize=label_fontsize)
plt.yticks(fontsize=tick_fontsize)
plt.ylim([8, 22])

plt.tight_layout()
plt.savefig(gwas_neglogpval_png_path)
plt.close()

#%% boxenplot ms
pairs = [(str(1), str(i)) for i in range(2, pleio_high_cutoff + 1)]
m_gwas_df[x] = m_gwas_df[x].astype(str)
ax = seaborn.boxenplot(x=x, y=y, data=m_gwas_df, showfliers=False, palette="rocket_r")

group1 = m_gwas_df.where(m_gwas_df.gwas_category_count == '1').dropna()[y]
group2 = m_gwas_df.where(m_gwas_df.gwas_category_count == '2').dropna()[y]
x1 = 0; x2 = 1; annot_y = 20; h = 1;
boxenplot_with_mannwhitneyu(group1, group2, x1, x2, annot_y, h)

group1 = m_gwas_df.where(m_gwas_df.gwas_category_count == '1').dropna()[y]
group2 = m_gwas_df.where(m_gwas_df.gwas_category_count == '3').dropna()[y]
x1 = 0; x2 = 2; annot_y = 25; h = 1;
boxenplot_with_mannwhitneyu(group1, group2, x1, x2, annot_y, h)

plt.title(title, fontsize=label_fontsize)
plt.xlabel(xlabel, fontsize=label_fontsize)
xticks_labels = [str(x) for x in (plt.xticks()[0] + 1)]
xticks_labels[-1] = '≥' + str(xticks_labels[-1])
plt.xticks(ticks=(plt.xticks()[0]), labels=xticks_labels, fontsize=tick_fontsize, rotation=0)
plt.ylabel(ylabel, fontsize=label_fontsize)
plt.yticks(fontsize=tick_fontsize)
plt.ylim([0, 30])

plt.tight_layout()
png_path = os.path.join(outdir_path, "gwas_signif_ms_boxenplot.png")
plt.savefig(png_path)
plt.close()

#%% boxenplot
pairs = [(str(1), str(i)) for i in range(2, pleio_high_cutoff + 1)]
m_gwas_df[x] = m_gwas_df[x].astype(str)
ax = seaborn.boxenplot(x=x, y=y, data=m_gwas_df, showfliers=True, palette="rocket_r")

annotator = Annotator(ax, pairs, data=m_gwas_df, x=x, y=y, order=order)
annotator.configure(test='Mann-Whitney', text_format='star', **annotator_config_dic)
annotator.apply_and_annotate()

plt.title(title, fontsize=label_fontsize)
plt.xlabel(xlabel, fontsize=label_fontsize)
xticks_labels = [str(x) for x in (plt.xticks()[0] + 1)]
xticks_labels[-1] = '≥' + str(xticks_labels[-1])
plt.xticks(ticks=(plt.xticks()[0]), labels=xticks_labels, fontsize=tick_fontsize, rotation=0)
plt.ylabel(ylabel, fontsize=label_fontsize)
plt.yticks(fontsize=tick_fontsize)

plt.tight_layout()
png_path = os.path.join(outdir_path, "gwas_signif_boxenplot.png")
plt.savefig(png_path)
plt.close()