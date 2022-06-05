"""This script will compute a tissue enrichment in a given GWAS"""
import os
import pathlib
import pandas
import seaborn as sns
from statannot import add_stat_annotation

from eqtl2gwas_pleiotropy.PathManager import PathManager
from eqtl2gwas_pleiotropy.URL import URL
from matplotlib import pyplot as plt
from eqtl2gwas_pleiotropy.constants import label_fontsize, tick_fontsize, dpi


#%% download gencode annotation
url = "https://remap.univ-amu.fr/storage/remap2022/hg38/MACS2/remap2022_crm_macs2_hg38_v1_0.bed.gz"
remap_crm_path = URL(url).download()

#%% outdir path
if not '__file__' in locals():
    __file__ = "plt_overlap_crm_remap.py"
outdir_path = os.path.join(PathManager.get_project_path(), "out", os.path.basename(__file__))
pathlib.Path(outdir_path).mkdir(parents=True, exist_ok=True)

#%% input dir cmpt_count_per_rsid
cmpt_count_per_rsid_indir_path = os.path.join(PathManager.get_outdir_path(), "cmpt_count_per_rsid.py")

#%% input dir cmpt_overlap_crm_remap
cmpt_overlap_crm_remap_indir_path = os.path.join(PathManager.get_outdir_path(), "cmpt_overlap_crm_remap.py")

# #%% bedtools intersect
# for count_pleio in range(1, 6):
#     bed_path = os.path.join(indir_path, "variant_pleio_{}_flank_0.bed".format(count_pleio))
#     intersect_bed_path = os.path.join(outdir_path, "remap_crm_" + os.path.basename(bed_path))
#     bedtool_intersect = "bedtools intersect -sorted -a {bed_path} -b {remap_crm_path} -wb"
#     bedtool_intersect = bedtool_intersect.format(**{'bed_path': bed_path, 'remap_crm_path': remap_crm_path, 'output_bed': intersect_bed_path})
#     Logger.info(bedtool_intersect)
#     with open(intersect_bed_path, 'w') as fout:
#         result = subprocess.run(shlex.split(bedtool_intersect), stdout=fout)

#%% barplot with proportion of variants in CRMs
x_lst = []
y_lst = []
for count_pleio in range(1, 6):
    bed_path = os.path.join(cmpt_count_per_rsid_indir_path, "variant_pleio_{}_flank_0.bed".format(count_pleio))
    intersect_bed_path = os.path.join(cmpt_overlap_crm_remap_indir_path, "remap_crm_" + os.path.basename(bed_path))
    count_variants = sum(1 for line in open(bed_path))
    count_intersection = sum(1 for line in open(intersect_bed_path))
    x_lst.append(count_pleio)
    y_lst.append(count_intersection/count_variants)

#%% barplot: Proportion varians in CRMs
plt.bar(x_lst, y_lst)
plt.grid(True)
plt.title("Variants and TF CRMs", fontsize=label_fontsize)
plt.xlabel("GWAS category count", fontsize=label_fontsize)
plt.xticks(fontsize=tick_fontsize)
plt.xticks(fontsize=tick_fontsize, rotation=0)
plt.ylabel("Ratio CRMs/variants", fontsize=label_fontsize)
plt.yticks(fontsize=tick_fontsize)

plt.tight_layout()
png_path = os.path.join(outdir_path, "remap_crm_bar_plot_ratio_variants.png")
plt.savefig(png_path, dpi=dpi)
plt.clf()
plt.close()

#%% Correlation CRMs complexity and pleiotropy
crm_out_df = pandas.DataFrame({'gwas_category_count': [], 'crm_tf_count': []})
for count_pleio in range(1, 6):
    bed_path = os.path.join(cmpt_count_per_rsid_indir_path, "variant_pleio_{}_flank_0.bed".format(count_pleio))
    intersect_bed_path = os.path.join(cmpt_overlap_crm_remap_indir_path, "remap_crm_" + os.path.basename(bed_path))
    crm_pleio_df = pandas.read_csv(intersect_bed_path, sep="\t", header=None)
    crm_pleio_df.rename({4: 'gwas_category_count', 10: 'crm_tf_count'}, axis=1, inplace=True)
    crm_pleio_df.loc[crm_pleio_df['crm_tf_count'].isna(), 'crm_tf_count'] = 0  # replace crm na with 0
    crm_out_df = pandas.concat([crm_out_df, crm_pleio_df[['gwas_category_count', 'crm_tf_count']]], axis=0)
    # count_variants = sum(1 for line in open(bed_path))
    # count_intersection = sum(1 for line in open(intersect_bed_path))
    # x_lst.append(count_pleio)
    # y_lst.append(count_intersection/count_variants)

crm_out_df['gwas_category_count'] = crm_out_df['gwas_category_count'].astype('int')
crm_out_df['crm_tf_count'] = crm_out_df['crm_tf_count'].astype('int')

sns.set_theme(style="whitegrid")
order = [1, 2, 3, 4, 5]
ax = sns.boxplot(x="gwas_category_count", y="crm_tf_count", data=crm_out_df, order=order)
test_results = add_stat_annotation(ax, data=crm_out_df, x="gwas_category_count", y="crm_tf_count", order=order, box_pairs=[(1, 2), (1, 3), (1, 4), (1, 5)], test='Mann-Whitney', text_format='star',loc='inside', verbose=2)

plt.xlabel("GWAS category count", fontsize=label_fontsize)
plt.xticks(fontsize=tick_fontsize, rotation=0)
plt.ylabel("CRM TF count", fontsize=label_fontsize)
plt.yticks(fontsize=tick_fontsize)

plt.tight_layout()
fig = ax.get_figure()
png_path = os.path.join(outdir_path, "boxplot_gwas_crm_tf.png")
fig.savefig(png_path, dpi=dpi)
plt.clf()
plt.close()

