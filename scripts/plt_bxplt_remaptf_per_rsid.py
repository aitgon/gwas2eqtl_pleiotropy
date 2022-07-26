import os
import sys

import numpy
import seaborn
import pandas
from matplotlib import pyplot as plt
from statannot import add_stat_annotation

#%%
from eqtl2gwas_pleiotropy.constants import label_fontsize, tick_fontsize

# upper_var_gwas_cat_count = 5

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
# i_path_lst.append("/home/gonzalez/Repositories/eqtl2gwas_pleiotropy/out/remap/pleio1.bed")
# i_path_lst.append("/home/gonzalez/Repositories/eqtl2gwas_pleiotropy/out/remap/pleio2.bed")
# i_path_lst.append("/home/gonzalez/Repositories/eqtl2gwas_pleiotropy/out/remap/pleio3.bed")
# i_path_lst.append("/home/gonzalez/Repositories/eqtl2gwas_pleiotropy/out/remap/pleio4.bed")
# i_path_lst.append("/home/gonzalez/Repositories/eqtl2gwas_pleiotropy/out/remap/pleio5.bed")

# i_path_lst.append("/home/gonzalez/Repositories/eqtl2gwas_pleiotropy/out/remap/pleio1_flank_50.bed")
# i_path_lst.append("/home/gonzalez/Repositories/eqtl2gwas_pleiotropy/out/remap/pleio2_flank_50.bed")
# i_path_lst.append("/home/gonzalez/Repositories/eqtl2gwas_pleiotropy/out/remap/pleio3_flank_50.bed")
# i_path_lst.append("/home/gonzalez/Repositories/eqtl2gwas_pleiotropy/out/remap/pleio4_flank_50.bed")
# i_path_lst.append("/home/gonzalez/Repositories/eqtl2gwas_pleiotropy/out/remap/pleio5_flank_50.bed")
#
# # vlnplt_png_path = "/home/gonzalez/Repositories/eqtl2gwas_pleiotropy/out/remap/bxplt.png"
# vlnplt_png_path = "/home/gonzalez/Repositories/eqtl2gwas_pleiotropy/out/remap/bxplt_flank_50.png"

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
order = [*range(1, upper_var_gwas_cat_count+1)]
seaborn.set_theme(style="whitegrid")
xticklabels = order.copy()
xticklabels[-1] = '≥{}'.format(order[-1])
title = "ReMap TFs per var."
xlabel = "GWAS category count"
ylabel = "TF count"
y = "tf"

#%%
box_pairs = [(1, i) for i in range(2, upper_var_gwas_cat_count+1) ]
ax = seaborn.violinplot(x="gwas_category_count", y=y, data=cat_df, order=order, palette="rocket_r")
test_results = add_stat_annotation(ax, data=cat_df, x="gwas_category_count", y=y, order=order,
                                   box_pairs=box_pairs,
                                   test='Mann-Whitney', text_format='star',
                                   loc='inside', verbose=2)

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
box_pairs = [(1, i) for i in range(2, upper_var_gwas_cat_count+1) ]
ax = seaborn.violinplot(x="gwas_category_count", y=y, data=cat_df, order=order, palette="rocket_r")
test_results = add_stat_annotation(ax, data=cat_df, x="gwas_category_count", y=y, order=order,
                                   box_pairs=box_pairs,
                                   test='Mann-Whitney', text_format='star',
                                   loc='inside', verbose=2)

plt.title(title, fontsize=label_fontsize)
plt.xlabel(xlabel, fontsize=label_fontsize)
plt.xticks(fontsize=tick_fontsize, rotation=0)
plt.ylabel(ylabel, fontsize=label_fontsize)
plt.yticks(fontsize=tick_fontsize)
ax.set_xticklabels(xticklabels)

plt.tight_layout()
plt.savefig(tf_flank_50_png)
plt.close()
