from various import boxenplot_with_mannwhitneyu
from constants import label_fontsize, tick_fontsize, annotator_config_dic, boxenplot_kws, \
    boxenplot_line_kws
from statannotations.Annotator import Annotator
from constants import seaborn_theme_dic

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

seaborn.set_theme(**seaborn_theme_dic)

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
    sample_size_min = int(sys.argv[9])
    sample_size_max = int(sys.argv[10])
    eur_af_min = float(sys.argv[11])
    eur_af_max = float(sys.argv[12])
    if len(sys.argv) > 13:
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

h4_df['sample_size'] = h4_df['ncontrol'] + h4_df['ncase']
#%%
count_per_rsid_gwas_df = pandas.read_excel(count_per_rsid_gwas_ods_path, engine='odf')

#%%
m_df = h4_df.merge(count_per_rsid_gwas_df, on=['chrom', 'pos38', 'rsid', 'ref', 'alt'])
# import pdb; pdb.set_trace()

#%% set sample size
m_df = m_df.query("sample_size >= {} & sample_size <= {}".format(sample_size_min, sample_size_max))
m_df = m_df.query("eur_af >= {} & eur_af <= {}".format(eur_af_min, eur_af_max))

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

#%% eqtl effect %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# boxplot and mann-whitney
y = "eqtl_beta_abs"
title = "eQTL effect"
ylabel = "Absolute beta"

#%%
pairs = [(str(1), str(i)) for i in range(2, pleio_high_cutoff + 1)]
m_eqtl_df[x] = m_eqtl_df[x].astype(str)

#%% boxenplot
ax = seaborn.boxplot(x=x, y=y, data=m_eqtl_df, **boxenplot_kws)

ylim = [0, 1.6]
x1_annot1 = 0.
delta_h = 0.03

y_annot1 = ylim[1] - 0.15 * ylim[1]
y_annot2 = y_annot1 - 0.15 * ylim[1]
h_annot = ylim[1] * delta_h
x2_annot1 = x1_annot1 + 1
x2_annot2 = x1_annot1 + 2

group1 = m_eqtl_df.where(m_eqtl_df.gwas_category_count == '1').dropna()[y]
group2 = m_eqtl_df.where(m_eqtl_df.gwas_category_count == '3').dropna()[y]
boxenplot_with_mannwhitneyu(group1, group2, x1_annot1, x2_annot2, y_annot1, h_annot)

group1 = m_eqtl_df.where(m_eqtl_df.gwas_category_count == '1').dropna()[y]
group2 = m_eqtl_df.where(m_eqtl_df.gwas_category_count == '2').dropna()[y]
boxenplot_with_mannwhitneyu(group1, group2, x1_annot1, x2_annot1, y_annot2, h_annot)

# plt.title(title, fontsize=label_fontsize, pad=100)
plt.title(title, fontsize=label_fontsize)
plt.xlabel(xlabel, fontsize=label_fontsize)
# xticks_labels = [str(x) for x in (plt.xticks()[0] + 1)]
xticks_labels = [str(x+1) for x in plt.xticks()[0]]
xticks_labels[-1] = '≥' + str(xticks_labels[-1])
plt.xticks(ticks=(plt.xticks()[0]), labels=xticks_labels, fontsize=tick_fontsize, rotation=0)
plt.ylabel(ylabel, fontsize=label_fontsize)
plt.yticks(fontsize=tick_fontsize)
plt.ylim(ylim)

# Calculate the medians directly from the data
medians = m_eqtl_df.groupby(x)[y].median()

# Annotate each median on the boxplot
for pos, median in enumerate(medians):
    x_loc = pos  # x location corresponds to the position of the box
    ax.annotate(f'{median:.2f}', xy=(x_loc, median), xytext=(0, 5),
                textcoords="offset points", ha='center', va='bottom', color='black', fontsize=tick_fontsize)

plt.tight_layout()
# plt.subplots_adjust(top=0.7)
png_path = os.path.join(eqtl_beta_png_path)
plt.savefig(png_path)
plt.close()

#%% gwas effect %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# boxplot and mann-whitney
y = "gwas_beta_abs"
title = "GWAS effect"
ylabel = "Absolute beta"
pairs = [(str(1), str(i)) for i in range(2, pleio_high_cutoff + 1)]
m_gwas_df[x] = m_gwas_df[x].astype(str)

#%% boxenplot custom
ax = seaborn.boxplot(x=x, y=y, data=m_gwas_df, **boxenplot_kws)

ylim = [-0.05, 0.45]
x1_annot1 = 0.
delta_h = 0.03

y_annot1 = ylim[1] - 0.15 * ylim[1]
y_annot2 = y_annot1 - 0.15 * ylim[1]
h_annot = ylim[1] * delta_h
x2_annot1 = x1_annot1 + 1
x2_annot2 = x1_annot1 + 2

group1 = m_gwas_df.where(m_gwas_df.gwas_category_count == '1').dropna()[y]
group2 = m_gwas_df.where(m_gwas_df.gwas_category_count == '3').dropna()[y]
boxenplot_with_mannwhitneyu(group1, group2, x1_annot1, x2_annot2, y_annot1, h_annot)

group1 = m_gwas_df.where(m_gwas_df.gwas_category_count == '1').dropna()[y]
group2 = m_gwas_df.where(m_gwas_df.gwas_category_count == '2').dropna()[y]
boxenplot_with_mannwhitneyu(group1, group2, x1_annot1, x2_annot1, y_annot2, h_annot)

# plt.title(title, fontsize=label_fontsize, pad=100)
plt.title(title, fontsize=label_fontsize)
plt.xlabel(xlabel, fontsize=label_fontsize)
# xticks_labels = [str(x) for x in (plt.xticks()[0] + 1)]
xticks_labels = [str(x+1) for x in plt.xticks()[0]]
xticks_labels[-1] = '≥' + str(xticks_labels[-1])
plt.xticks(ticks=(plt.xticks()[0]), labels=xticks_labels, fontsize=tick_fontsize, rotation=0)
plt.ylabel(ylabel, fontsize=label_fontsize)
plt.yticks(fontsize=tick_fontsize)
plt.ylim(ylim)

# Calculate the medians directly from the data
medians = m_gwas_df.groupby(x)[y].median()

# Annotate each median on the boxplot
for pos, median in enumerate(medians):
    x_loc = pos  # x location corresponds to the position of the box
    ax.annotate(f'{median:.2f}', xy=(x_loc, median), xytext=(0, 5),
                textcoords="offset points", ha='center', va='bottom', color='black', fontsize=tick_fontsize)

plt.tight_layout()
# plt.subplots_adjust(top=0.7)
png_path = os.path.join(gwas_beta_png_path)
plt.savefig(png_path)
plt.close()

#%% eqtl signif %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# boxplot and mann-whitney
y = "eqtl_neglog10pval"
title = "eQTL significance"
ylabel = "Neg. log10 p-val"

pairs = [(str(1), str(i)) for i in range(2, pleio_high_cutoff + 1)]
m_eqtl_df[x] = m_eqtl_df[x].astype(str)

#%% boxenplot ms
ax = seaborn.boxplot(x=x, y=y, data=m_eqtl_df, **boxenplot_kws)

ylim = [0, 40]
x1_annot1 = 0.
delta_h = 0.03

y_annot1 = ylim[1] - 0.15 * ylim[1]
y_annot2 = y_annot1 - 0.15 * ylim[1]
h_annot = ylim[1] * delta_h
x2_annot1 = x1_annot1 + 1
x2_annot2 = x1_annot1 + 2

group1 = m_eqtl_df.where(m_eqtl_df.gwas_category_count == '1').dropna()[y]
group2 = m_eqtl_df.where(m_eqtl_df.gwas_category_count == '3').dropna()[y]
boxenplot_with_mannwhitneyu(group1, group2, x1_annot1, x2_annot2, y_annot1, h_annot)

group1 = m_eqtl_df.where(m_eqtl_df.gwas_category_count == '1').dropna()[y]
group2 = m_eqtl_df.where(m_eqtl_df.gwas_category_count == '2').dropna()[y]
boxenplot_with_mannwhitneyu(group1, group2, x1_annot1, x2_annot1, y_annot2, h_annot)

# plt.title(title, fontsize=label_fontsize, pad=100)
plt.title(title, fontsize=label_fontsize)
plt.xlabel(xlabel, fontsize=label_fontsize)
# xticks_labels = [str(x) for x in (plt.xticks()[0] + 1)]
xticks_labels = [str(x+1) for x in plt.xticks()[0]]
xticks_labels[-1] = '≥' + str(xticks_labels[-1])
plt.xticks(ticks=(plt.xticks()[0]), labels=xticks_labels, fontsize=tick_fontsize, rotation=0)
plt.ylabel(ylabel, fontsize=label_fontsize)
plt.yticks(fontsize=tick_fontsize)
plt.ylim(ylim)

# Calculate the medians directly from the data
medians = m_eqtl_df.groupby(x)[y].median()

# Annotate each median on the boxplot
for pos, median in enumerate(medians):
    x_loc = pos  # x location corresponds to the position of the box
    ax.annotate(f'{median:.2f}', xy=(x_loc, median), xytext=(0, 5),
                textcoords="offset points", ha='center', va='bottom', color='black', fontsize=tick_fontsize)

plt.tight_layout()
# plt.subplots_adjust(top=0.7)
png_path = os.path.join(eqtl_neglogpval_png_path)
plt.savefig(png_path)
plt.close()

#%% gwas signif %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# boxplot and mann-whitney
y = "gwas_neglog10pval"
title = "GWAS significance"
ylabel = "Neg. log10 p-val"

pairs = [(str(1), str(i)) for i in range(2, pleio_high_cutoff + 1)]
m_gwas_df[x] = m_gwas_df[x].astype(str)

#%% boxenplot ms
ax = seaborn.boxplot(x=x, y=y, data=m_gwas_df, **boxenplot_kws)

ylim = [5, 35]
x1_annot1 = 0.
delta_h = 0.03

y_annot1 = ylim[1] - 0.15 * ylim[1]
y_annot2 = y_annot1 - 0.15 * ylim[1]
h_annot = ylim[1] * delta_h
x2_annot1 = x1_annot1 + 1
x2_annot2 = x1_annot1 + 2

group1 = m_gwas_df.where(m_gwas_df.gwas_category_count == '1').dropna()[y]
group2 = m_gwas_df.where(m_gwas_df.gwas_category_count == '3').dropna()[y]
boxenplot_with_mannwhitneyu(group1, group2, x1_annot1, x2_annot2, y_annot1, h_annot)

group1 = m_gwas_df.where(m_gwas_df.gwas_category_count == '1').dropna()[y]
group2 = m_gwas_df.where(m_gwas_df.gwas_category_count == '2').dropna()[y]
boxenplot_with_mannwhitneyu(group1, group2, x1_annot1, x2_annot1, y_annot2, h_annot)

plt.title(title, fontsize=label_fontsize)
plt.xlabel(xlabel, fontsize=label_fontsize)
# xticks_labels = [str(x) for x in (plt.xticks()[0] + 1)]
xticks_labels = [str(x+1) for x in plt.xticks()[0]]
xticks_labels[-1] = '≥' + str(xticks_labels[-1])
plt.xticks(ticks=(plt.xticks()[0]), labels=xticks_labels, fontsize=tick_fontsize, rotation=0)
plt.ylabel(ylabel, fontsize=label_fontsize)
plt.yticks(fontsize=tick_fontsize)
plt.ylim(ylim)

# Calculate the medians directly from the data
medians = m_gwas_df.groupby(x)[y].median()

# Annotate each median on the boxplot
for pos, median in enumerate(medians):
    x_loc = pos  # x location corresponds to the position of the box
    ax.annotate(f'{median:.2f}', xy=(x_loc, median), xytext=(0, 5),
                textcoords="offset points", ha='center', va='bottom', color='black', fontsize=tick_fontsize)

plt.tight_layout()
# plt.subplots_adjust(top=0.7)
png_path = os.path.join(gwas_neglogpval_png_path)
plt.savefig(png_path)
plt.close()
