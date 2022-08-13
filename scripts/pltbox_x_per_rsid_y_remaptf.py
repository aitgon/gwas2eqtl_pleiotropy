import os
import pandas
import seaborn
import sys

from eqtl2gwas_pleiotropy.constants import label_fontsize, tick_fontsize, boxplot_kwargs, annotator_config_dic
from matplotlib import pyplot as plt
from statannotations.Annotator import Annotator
from eqtl2gwas_pleiotropy.constants import seaborn_theme_dic

#%%

seaborn.set_theme(**seaborn_theme_dic)

#%%
help_cmd_str = "todo"
try:
    remap_nr_pleio_1_flank_0_hg38_bed = sys.argv[1]
    remap_nr_pleio_1_flank_50_hg38_bed = sys.argv[2]
    upper_var_gwas_cat_count = int(sys.argv[3])
    tf_flank_0_png = sys.argv[4]
    tf_flank_50_png = sys.argv[5]
    if len(sys.argv) > 6:
        print("""Two many arguments!
        {}""".format(help_cmd_str))
        sys.exit(1)
except IndexError:
    print("""Argument missing!
    {}""".format(help_cmd_str))
    sys.exit(1)

i_path_lst = []

indir_path = os.path.dirname(remap_nr_pleio_1_flank_0_hg38_bed)


########################################################################################################################
flank = 0

cat_df = pandas.DataFrame({'gwas_category_count': [], 'rsid': [], 'tf': []})
for pleio in range(1, upper_var_gwas_cat_count+1):
    pleio_path = os.path.join(indir_path, "remap_nr_variant_pleio_{}_flank_{}_hg38.bed".format(pleio, flank))
    df = pandas.read_csv(pleio_path, sep="\t", header=None)
    df['tf'] = df[9].str.split(':', expand=True)[0]
    df = df[[3, 'tf']].drop_duplicates()
    df.columns = ['rsid', 'tf']
    df = df.groupby('rsid').count().reset_index()
    df['gwas_category_count'] = pleio
    cat_df = pandas.concat([cat_df, df], axis=0)

#%%
order = [str(x) for x in range(1, upper_var_gwas_cat_count+1)]
xticklabels = order.copy()
xticklabels[-1] = '≥{}'.format(order[-1])
title = "TFs per variant"
xlabel = "GWAS category count"
ylabel = "TF count"
y = "tf"
x = "gwas_category_count"

#%%
pairs = [(str(1), str(i)) for i in range(2, upper_var_gwas_cat_count + 1)]
cat_df[x] = cat_df[x].astype(int).astype(str)
ax = seaborn.boxplot(x=x, y=y, data=cat_df, order=order, **boxplot_kwargs)
annotator = Annotator(ax, pairs, data=cat_df, x=x, y=y, order=order)
annotator.configure(test='Mann-Whitney', text_format='star', **annotator_config_dic)
annotator.apply_and_annotate()

# box_pairs = [(1, i) for i in range(2, upper_var_gwas_cat_count+1) ]
# ax = seaborn.violinplot(x="gwas_category_count", y=y, data=cat_df, order=order, palette="rocket_r")
# test_results = add_stat_annotation(ax, data=cat_df, x="gwas_category_count", y=y, order=order,
#                                    box_pairs=box_pairs,
#                                    test='Mann-Whitney', text_format='star',
#                                    loc='inside', verbose=2)

plt.title(title, fontsize=label_fontsize)
plt.xlabel(xlabel, fontsize=label_fontsize)
plt.xticks(fontsize=tick_fontsize, rotation=0)
plt.ylabel(ylabel, fontsize=label_fontsize)
plt.yticks(fontsize=tick_fontsize)
ax.set_xticklabels(xticklabels)

plt.tight_layout()
plt.savefig(tf_flank_0_png)
plt.close()

########################################################################################################################
flank = 50

cat_df = pandas.DataFrame({'gwas_category_count': [], 'rsid': [], 'tf': []})
for pleio in range(1, upper_var_gwas_cat_count+1):
    pleio_path = os.path.join(indir_path, "remap_nr_variant_pleio_{}_flank_{}_hg38.bed".format(pleio, flank))
    df = pandas.read_csv(pleio_path, sep="\t", header=None)
    df['tf'] = df[9].str.split(':', expand=True)[0]
    df = df[[3, 'tf']].drop_duplicates()
    df.columns = ['rsid', 'tf']
    df = df.groupby('rsid').count().reset_index()
    df['gwas_category_count'] = pleio
    cat_df = pandas.concat([cat_df, df], axis=0)

#%%
# box_pairs = [(1, i) for i in range(2, upper_var_gwas_cat_count+1) ]
pairs = [(str(1), str(i)) for i in range(2, upper_var_gwas_cat_count + 1)]
cat_df[x] = cat_df[x].astype(int).astype(str)
ax = seaborn.boxplot(x=x, y=y, data=cat_df, order=order, **boxplot_kwargs)
annotator = Annotator(ax, pairs, data=cat_df, x=x, y=y, order=order)
annotator.configure(test='Mann-Whitney', text_format='star', **annotator_config_dic)
annotator.apply_and_annotate()

# ax = seaborn.violinplot(x="gwas_category_count", y=y, data=cat_df, order=order, palette="rocket_r")
# test_results = add_stat_annotation(ax, data=cat_df, x="gwas_category_count", y=y, order=order,
#                                    box_pairs=pairs,
#                                    test='Mann-Whitney', text_format='star',
#                                    loc='inside', verbose=2)

plt.title(title, fontsize=label_fontsize)
plt.xlabel(xlabel, fontsize=label_fontsize)
plt.xticks(fontsize=tick_fontsize, rotation=0)
plt.ylabel(ylabel, fontsize=label_fontsize)
plt.yticks(fontsize=tick_fontsize)
ax.set_xticklabels(xticklabels)

plt.tight_layout()
plt.savefig(tf_flank_50_png)
plt.close()
