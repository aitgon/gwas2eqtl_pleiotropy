import os
import pathlib
import shlex
import subprocess
import sys

from eqtl2gwas_pleiotropy.Logger import Logger
from eqtl2gwas_pleiotropy.constants import public_data_dir

#%%
help_cmd_str = "todo"
try:
    count_per_rsid_gwas_tsv_path = sys.argv[1]
    variant_pleio_1_flank_0_hg38_bed = sys.argv[2]
    variant_pleio_1_flank_50_hg38_bed = sys.argv[3]
    upper_var_gwas_cat_count = int(sys.argv[4])
    remap_nr_pleio_1_flank_0_hg38_bed = sys.argv[5]
    remap_nr_pleio_1_flank_50_hg38_bed = sys.argv[6]
    if len(sys.argv) > 7:
        print("""Two many arguments!
        {}""".format(help_cmd_str))
        sys.exit(1)
except IndexError:
    print("""Argument missing!
    {}""".format(help_cmd_str))
    sys.exit(1)

# #%% outdir path
# if not '__file__' in locals():
#     __file__ = "intrsct_remap_nr_gwas_categories.py"

outdir_path = os.path.join(os.path.dirname(remap_nr_pleio_1_flank_0_hg38_bed))
pathlib.Path(outdir_path).mkdir(parents=True, exist_ok=True)
outdir_path = os.path.join(os.path.dirname(remap_nr_pleio_1_flank_50_hg38_bed))
pathlib.Path(outdir_path).mkdir(parents=True, exist_ok=True)

#%% input dir cmpt_count_per_rsid
indir_path = os.path.dirname(count_per_rsid_gwas_tsv_path)

remap_crm_path = os.path.join(public_data_dir, "remap.univ-amu.fr/storage/remap2022/hg38/MACS2/remap2022_nr_macs2_hg38_v1_0.bed.gz")


# def f(flank, count_pleio):
#     bed_path = os.path.join(indir_path, "variant_pleio_{}_flank_{}_hg38.bed".format(count_pleio, flank))
#     intersect_bed_path = os.path.join(outdir_path, "remap_crm_" + os.path.basename(bed_path))
#     cmd_stf = "bedtools intersect -sorted -a {bed_path} -b {remap_crm_path} -wb"
#     cmd = cmd_stf.format(**{'bed_path': bed_path, 'remap_crm_path': remap_crm_path, 'output_bed': intersect_bed_path})
#     Logger.info(cmd)
#     with open(intersect_bed_path, 'w') as fout:
#         result = subprocess.run(shlex.split(cmd), stdout=fout)

#%% bedtools intersect
flank = 0
for count_pleio in range(1, upper_var_gwas_cat_count+1):
    bed_path = os.path.join(indir_path, "variant_pleio_{}_flank_{}_hg38.bed".format(count_pleio, flank))
    intersect_basename = os.path.basename(remap_nr_pleio_1_flank_0_hg38_bed).replace('pleio_1', 'pleio_' + str(count_pleio))
    intersect_bed_path = os.path.join(os.path.dirname(remap_nr_pleio_1_flank_0_hg38_bed), intersect_basename)
    cmd_stf = "bedtools intersect -sorted -a {bed_path} -b {remap_crm_path} -wb"
    cmd = cmd_stf.format(**{'bed_path': bed_path, 'remap_crm_path': remap_crm_path, 'output_bed': intersect_bed_path})
    Logger.info(cmd)
    with open(intersect_bed_path, 'w') as fout:
        result = subprocess.run(shlex.split(cmd), stdout=fout)

#%% bedtools intersect
flank = 50
for count_pleio in range(1, upper_var_gwas_cat_count+1):
    bed_path = os.path.join(indir_path, "variant_pleio_{}_flank_{}_hg38.bed".format(count_pleio, flank))
    intersect_basemane = os.path.basename(remap_nr_pleio_1_flank_50_hg38_bed).replace('pleio_1', 'pleio_' + str(count_pleio))
    intersect_bed_path = os.path.join(os.path.dirname(remap_nr_pleio_1_flank_50_hg38_bed), intersect_basemane)
    cmd_stf = "bedtools intersect -sorted -a {bed_path} -b {remap_crm_path} -wb"
    cmd = cmd_stf.format(**{'bed_path': bed_path, 'remap_crm_path': remap_crm_path, 'output_bed': intersect_bed_path})
    Logger.info(cmd)
    with open(intersect_bed_path, 'w') as fout:
        result = subprocess.run(shlex.split(cmd), stdout=fout)
