"""For each rsid pleiotropy group, this script compute the vep consequence"""
from eqtl2gwas_pleiotropy.PathManager import PathManager
from eqtl2gwas_pleiotropy.constants import region_bin, label_fontsize, tick_fontsize
from matplotlib import pyplot as plt

import sys
import math
import numpy
import os
import pandas
import pathlib
import shlex
import subprocess


#%%
help_cmd_str = "todo"
try:
    count_per_rsid_gwas_tsv_path = sys.argv[1]
    cmd_path = sys.argv[2]
    if len(sys.argv) > 3:
        print("""Two many arguments!
        {}""".format(help_cmd_str))
        sys.exit(1)
except IndexError:
    print("""Argument missing!
    {}""".format(help_cmd_str))
    sys.exit(1)

outdir_path = os.path.dirname(cmd_path)
pathlib.Path(outdir_path).mkdir(parents=True, exist_ok=True)

#%%
# count_per_rsid_gwas_tsv_path = os.path.join(PathManager.get_outdir_path(), "cmpt_count_per_rsid.py", "count_per_rsid_gwas.tsv")
df = pandas.read_csv(count_per_rsid_gwas_tsv_path, sep="\t")

#%%
cmd_lst = []
for pleio in df['gwas_category_count'].unique():
    rsid_path = os.path.join(outdir_path, "rsid_pleio{}.txt".format(pleio))
    df.loc[df['gwas_category_count']==pleio, 'rsid'].to_csv(rsid_path, header=False, index=False)
    vep_path = os.path.join(outdir_path, "vep_pleio{}.tsv".format(pleio))
    cmd_str = "vep -i " + rsid_path + " --cache --output_file " + vep_path
    cmd_lst.append(cmd_str)
    # cmd_str_path = os.path.join(outdir_path, "_pleio{}.tsv".format(pleio))
    # stdout = subprocess.run(cmd_str, shell=True, capture_output=True)

with open(cmd_path, "w") as fout:
    for item in cmd_lst:
        fout.write(item + "\n")


