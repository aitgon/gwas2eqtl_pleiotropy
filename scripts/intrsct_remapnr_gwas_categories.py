import os
import pathlib
import shlex
import subprocess
import sys
from multiprocessing import Pool

from gwas2eqtl_pleiotropy.Logger import Logger

#%%
help_cmd_str = "todo"
try:
    threads = int(sys.argv[1])
    max_gwas_category_count = int(sys.argv[2])
    remap_nr_path = sys.argv[3]
    eqtl_pleio_1_flank_10_hg38_bed = sys.argv[4]
    remap_nr_pleio_1_flank_10_hg38_bed = sys.argv[5]
    if len(sys.argv) > 6:
        print("""Two many arguments!
        {}""".format(help_cmd_str))
        sys.exit(1)
except IndexError:
    print("""Argument missing!
    {}""".format(help_cmd_str))
    sys.exit(1)

outdir_path = os.path.join(os.path.dirname(remap_nr_pleio_1_flank_10_hg38_bed))
pathlib.Path(outdir_path).mkdir(parents=True, exist_ok=True)

#%% input dir cmpt_count_per_rsid
indir_path = os.path.dirname(eqtl_pleio_1_flank_10_hg38_bed)

#%% bedtools intersect
flank = 10

for gwas_category_count_i in range(max_gwas_category_count + 1):
    bed_path = os.path.join(indir_path, "eqtl_pleio_{}_flank_{}_hg38.bed".format(gwas_category_count_i, flank))
    if not os.path.isfile(bed_path):
        break

# max_gwas_category_count = max_gwas_category_count - 1
# print(max_gwas_category_count)
def intrsct_remapnr(count_pleio):
    bed_path = os.path.join(indir_path, "eqtl_pleio_{}_flank_{}_hg38.bed".format(count_pleio, flank))
    intersect_basename = os.path.basename(remap_nr_pleio_1_flank_10_hg38_bed).replace('pleio_1', 'pleio_' + str(count_pleio))
    intersect_bed_path = os.path.join(os.path.dirname(remap_nr_pleio_1_flank_10_hg38_bed), intersect_basename)
    cmd_stf = "bedtools intersect -sorted -a {bed_path} -b {remap_nr_path} -wb"
    cmd = cmd_stf.format(**{'bed_path': bed_path, 'remap_nr_path': remap_nr_path, 'output_bed': intersect_bed_path})
    Logger.info(cmd)
    with open(intersect_bed_path, 'w') as fout:
        result = subprocess.run(shlex.split(cmd), stdout=fout)

with Pool(processes=threads) as p:
    p.map(intrsct_remapnr, range(1, max_gwas_category_count + 1))
