"""This script takes variants in each pleiotropic group,
removes redundant variants in pleiotropic regions,
and computes the percentage of variants belonging to each category"""

import os
import pandas
import pathlib
import shlex
import subprocess
import seaborn

from gwas2eqtl_pleiotropy.Logger import Logger
from gwas2eqtl_pleiotropy.PathManager import PathManager
from matplotlib import pyplot as plt
from gwas2eqtl_pleiotropy.constants import seaborn_theme_dic
from gwas2eqtl_pleiotropy.constants import dpi, tick_fontsize, label_fontsize

plt.rcParams["figure.figsize"] = (16, 6)
seaborn.set_theme(**seaborn_theme_dic)

#%% outdir path
if not '__file__' in locals():
    __file__ = "plt_bar_nonredundant_variant_category_percentage.py"
outdir_path = os.path.join(PathManager.get_project_path(), "out", os.path.basename(__file__))
pathlib.Path(outdir_path).mkdir(parents=True, exist_ok=True)

#%% input region window
region_window_bed_path = os.path.join(PathManager.get_outdir_path(), "cmpt_pleiotropic_regions.py/region_window_100000.bed")

#%% input h4_annotated
h4_annotated_tsv_path = os.path.join(PathManager.get_outdir_path(), "annotate_db.py/h4_annotated.tsv")
h4_annotated_df = pandas.read_csv(h4_annotated_tsv_path, sep="\t", header=0)
gwas_cat_lst = sorted(h4_annotated_df['gwas_subcategory'].unique().tolist())

x_lst = [x[0:6] for x in gwas_cat_lst]

out_lst_2d_lst = []
for pleiotropy in range(1, 6):

    #%% input variant pleio
    variant_pleio_flank_0_bed_path = os.path.join(PathManager.get_outdir_path(), "cmpt_count_per_rsid.py/variant_pleio_{}_flank_0.bed".format(pleiotropy))

    #%%
    variant_region_bed = "variant_region.bed"

    #%% intersect
    cmd = "bedtools intersect -a {} -b {} -wb"
    cmd = cmd.format(variant_pleio_flank_0_bed_path, region_window_bed_path)
    Logger.info(cmd)
    with open(variant_region_bed, 'w') as fout:
        result = subprocess.run(shlex.split(cmd), stdout=fout)

    #%%
    df_variant_region = pandas.read_csv(variant_region_bed, header=None, sep="\t")[range(9)]

    #%%
    df_variant_region.columns = ['chrom', 'start', 'end', 'rsid', 'gwas_subcategory_count', 'gwas_subcategory_lst', 'region_chrom', 'region_start', 'region_end']

    #%%
    df_variant_region = df_variant_region.drop_duplicates(['region_chrom', 'region_start', 'region_end'])
    loci_variant_count = df_variant_region.shape[0]

    #%%
    df_long = df_variant_region.copy()
    df_long["gwas_subcategory_lst"] = df_long["gwas_subcategory_lst"].str.split(",")
    df_long = df_long.explode("gwas_subcategory_lst")

    y_lst = []
    for category_i, category in enumerate(gwas_cat_lst):
        cat_variant_count = df_long.loc[(df_long['gwas_subcategory_lst'] == category)].shape[0]
        cat_variant_perc = round(cat_variant_count / loci_variant_count * 100)
        out_lst_2d_lst.append(["Pleio. {}".format(pleiotropy), x_lst[category_i], cat_variant_perc])

variant_cat_perc_df = pandas.DataFrame(out_lst_2d_lst, columns=['Pleiotropy', 'Category', 'variant_cat_perc'])

ax = seaborn.barplot(x="Category", y="variant_cat_perc", hue="Pleiotropy", data=variant_cat_perc_df)

plt.grid(axis='y')
plt.legend(fontsize=16) # using a size in points
plt.title("Variant category and pleiotropy", fontsize=label_fontsize)
plt.xlabel("Category", fontsize=label_fontsize)
plt.xticks(fontsize=tick_fontsize, rotation=0)
plt.xticks(fontsize=tick_fontsize, rotation=90)
plt.ylabel("Variant category [%]", fontsize=label_fontsize)
plt.yticks(fontsize=tick_fontsize)

plt.tight_layout()
png_path = os.path.join(outdir_path, "non_redundant_variant_category_pleio2.png")
plt.savefig(png_path, dpi=dpi)
plt.close()
