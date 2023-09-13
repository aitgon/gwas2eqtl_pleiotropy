import os
import pathlib

import numpy
import pandas
import seaborn
import sys

from gwas2eqtl_pleiotropy import boxenplot_with_mannwhitneyu

from gwas2eqtl_pleiotropy.constants import label_fontsize, tick_fontsize, boxplot_kwargs, annotator_config_dic, \
    boxenplot_kws, boxenplot_line_kws
from matplotlib import pyplot as plt
from statannotations.Annotator import Annotator
from gwas2eqtl_pleiotropy.constants import seaborn_theme_dic

#%%

seaborn.set_theme(**seaborn_theme_dic)

#%%
help_cmd_str = "todo"
try:
    pleio_high_cutoff = int(sys.argv[1])
    remap_nr_pleio_1_flank_10_hg38_bed = sys.argv[2]
    tf_flank_10_png = sys.argv[3]
    if len(sys.argv) > 4:
        print("""Two many arguments!
        {}""".format(help_cmd_str))
        sys.exit(1)
except IndexError:
    print("""Argument missing!
    {}""".format(help_cmd_str))
    sys.exit(1)

i_path_lst = []

indir_path = os.path.dirname(remap_nr_pleio_1_flank_10_hg38_bed)
outdir_path = os.path.dirname(tf_flank_10_png)
pathlib.Path(outdir_path).mkdir(exist_ok=True, parents=True)

########################################################################################################################
flank = 10

cat_df = pandas.DataFrame({'gwas_category_count': [], 'rsid': [], 'tf': []})
# for pleio in range(1, max_gwas_category_count+1):
for pleio in range(1, pleio_high_cutoff + 1):
    pleio_path = remap_nr_pleio_1_flank_10_hg38_bed.replace('pleio_1', 'pleio_{}'.format(pleio))
    if not os.path.isfile(pleio_path):
        break
    if os.stat(pleio_path).st_size == 0:
        continue
    df = pandas.read_csv(pleio_path, sep="\t", header=None)
    df['tf'] = df[8].str.split(':', expand=True)[0]
    df = df[[3, 'tf']].drop_duplicates()
    df.columns = ['rsid', 'tf']
    df = df.groupby('rsid').count().reset_index()
    df['gwas_category_count'] = pleio
    cat_df = pandas.concat([cat_df, df], axis=0)

cat_df.loc[cat_df['gwas_category_count'] >= pleio_high_cutoff, "gwas_category_count"] = pleio_high_cutoff

describe_tsv_path = os.path.join(os.path.dirname(tf_flank_10_png), "describe.tsv")
cat_df.groupby('gwas_category_count')['tf'].apply(lambda x: x.describe()).to_csv(describe_tsv_path, sep="\t")


#%%
order = [str(int(x)) for x in cat_df['gwas_category_count'].unique()]
xticklabels = order.copy()
title = "TFs per variant"
xlabel = "Trait category count"
ylabel = "TF count"
y = "tf"
x = "gwas_category_count"

#%%
pairs = [(str(1), str(int(i))) for i in sorted(cat_df['gwas_category_count'].unique())]
cat_df[x] = cat_df[x].astype(int).astype(str)

#%% boxenplot
ax = seaborn.boxenplot(x=x, y=y, data=cat_df, order=order, **boxenplot_kws, line_kws=boxenplot_line_kws)

annotator = Annotator(ax, pairs, data=cat_df, x=x, y=y, order=order)
annotator.configure(test='Mann-Whitney', text_format='star', **annotator_config_dic)
annotator.apply_and_annotate()

plt.title(title, fontsize=label_fontsize)
plt.xlabel(xlabel, fontsize=label_fontsize)
# plt.xticks(fontsize=tick_fontsize, rotation=0)
xticks_labels = [str(x) for x in (plt.xticks()[0] + 1)]
xticks_labels[-1] = '≥' + str(xticks_labels[-1])
plt.xticks(ticks=(plt.xticks()[0]), labels=xticks_labels, fontsize=tick_fontsize, rotation=0)
plt.ylabel(ylabel, fontsize=label_fontsize)
plt.yticks(fontsize=tick_fontsize)
ax.set_xticklabels(xticklabels)

plt.tight_layout()
plt.savefig(tf_flank_10_png)
plt.close()

#%% boxenplot ms
ax = seaborn.boxenplot(x=x, y=y, data=cat_df, order=order, **boxenplot_kws, line_kws=boxenplot_line_kws)

ylim = [0, 100]
x1_annot1 = 0.
delta_h = 0.03

y_annot1 = ylim[1] - 0.15 * ylim[1]
y_annot2 = y_annot1 - 0.15 * ylim[1]
h_annot = ylim[1] * delta_h
x2_annot1 = x1_annot1 + 1
x2_annot2 = x1_annot1 + 2

group1 = cat_df.where(cat_df.gwas_category_count == '1').dropna()[y]
group2 = cat_df.where(cat_df.gwas_category_count == '3').dropna()[y]
boxenplot_with_mannwhitneyu(group1, group2, x1_annot1, x2_annot2, y_annot1, h_annot)

group1 = cat_df.where(cat_df.gwas_category_count == '1').dropna()[y]
group2 = cat_df.where(cat_df.gwas_category_count == '2').dropna()[y]
boxenplot_with_mannwhitneyu(group1, group2, x1_annot1, x2_annot1, y_annot2, h_annot)

plt.title(title, fontsize=label_fontsize)
plt.xlabel(xlabel, fontsize=label_fontsize)
xticks_labels = [str(x) for x in (plt.xticks()[0] + 1)]
xticks_labels[-1] = '≥' + str(xticks_labels[-1])
plt.xticks(ticks=(plt.xticks()[0]), labels=xticks_labels, fontsize=tick_fontsize, rotation=0)
plt.ylabel(ylabel, fontsize=label_fontsize)
plt.yticks(fontsize=tick_fontsize)
ax.set_xticklabels(xticklabels)
plt.ylim(ylim)

plt.tight_layout()
tf_flank_10_png = os.path.join(outdir_path, "boxenplot_custom.png")
plt.savefig(tf_flank_10_png)
plt.close()
