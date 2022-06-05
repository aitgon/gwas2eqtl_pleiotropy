"""This script will compute a tissue enrichment in a given GWAS"""
import os
import pathlib
import shlex
import subprocess

from eqtl2gwas_pleiotropy.Logger import Logger
from eqtl2gwas_pleiotropy.PathManager import PathManager
from eqtl2gwas_pleiotropy.URL import URL
from matplotlib import pyplot as plt
from eqtl2gwas_pleiotropy.constants import label_fontsize, tick_fontsize, dpi


#%% download gencode annotation
url = "https://remap.univ-amu.fr/storage/remap2022/hg38/MACS2/remap2022_crm_macs2_hg38_v1_0.bed.gz"
remap_crm_path = URL(url).download()

#%% outdir path
if not '__file__' in locals():
    __file__ = "cmpt_overlap_crm_remap.py"
outdir_path = os.path.join(PathManager.get_project_path(), "out", os.path.basename(__file__))
pathlib.Path(outdir_path).mkdir(parents=True, exist_ok=True)

#%% input dir cmpt_count_per_rsid
indir_path = os.path.join(PathManager.get_outdir_path(), "cmpt_count_per_rsid.py")

#%% bedtools intersect
for count_pleio in range(1, 6):
    bed_path = os.path.join(indir_path, "variant_pleio_{}_flank_0.bed".format(count_pleio))
    intersect_bed_path = os.path.join(outdir_path, "remap_crm_" + os.path.basename(bed_path))
    bedtool_intersect = "bedtools intersect -sorted -a {bed_path} -b {remap_crm_path} -wb"
    bedtool_intersect = bedtool_intersect.format(**{'bed_path': bed_path, 'remap_crm_path': remap_crm_path, 'output_bed': intersect_bed_path})
    Logger.info(bedtool_intersect)
    with open(intersect_bed_path, 'w') as fout:
        result = subprocess.run(shlex.split(bedtool_intersect), stdout=fout)

#%% barplot with proportion of variants in CRMs
x_lst = []
y_lst = []
for count_pleio in range(1, 6):
    bed_path = os.path.join(indir_path, "variant_pleio_{}_flank_0.bed".format(count_pleio))
    intersect_bed_path = os.path.join(outdir_path, "remap_crm_" + os.path.basename(bed_path))
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


