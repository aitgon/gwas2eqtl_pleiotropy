"""Variant effect predictor"""
import sqlalchemy

from gwas2eqtl_pleiotropy.Logger import Logger
from gwas2eqtl_pleiotropy.PathManager import PathManager
from gwas2eqtl_pleiotropy.constants import region_bin, label_fontsize, tick_fontsize
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
    snp_pp_h4 = float(sys.argv[1])
    url = sys.argv[2]
    count_per_rsid_gwas_ods_path = sys.argv[3]
    vep_cache_info = sys.argv[4]
    vep_input_path = sys.argv[5]
    vep_output_path = sys.argv[6]
    if len(sys.argv) > 7:
        print("""Two many arguments!
        {}""".format(help_cmd_str))
        sys.exit(1)
except IndexError:
    print("""Argument missing!
    {}""".format(help_cmd_str))
    sys.exit(1)

outdir_path = os.path.dirname(vep_input_path)
pathlib.Path(outdir_path).mkdir(parents=True, exist_ok=True)

#%%
sql = 'select * from colocpleio where snp_pp_h4>={}'.format(snp_pp_h4)
# columns = ['rsid', 'eqtl_beta', 'eqtl_gene_id', 'gwas_id', 'eqtl_id']
# h4_df = pandas.read_sql(sql, con=url).drop_duplicates()
engine = sqlalchemy.create_engine(url)
with engine.begin() as conn:
    h4_df = pandas.read_sql(sqlalchemy.text(sql), con=conn).drop_duplicates()

#%%
# gwas_cat_df = pandas.read_csv(count_per_rsid_gwas_ods_path, sep="\t")
gwas_cat_df = pandas.read_excel(count_per_rsid_gwas_ods_path, engine='odf')

#%%
vep_df = h4_df[['chrom', 'pos38', 'rsid', 'ref', 'alt']].merge(gwas_cat_df[['chrom', 'pos38', 'rsid', 'gwas_category_count']], on=['chrom', 'pos38', 'rsid']).drop_duplicates()

vep_df.rename({'pos38': 'start'}, axis=1, inplace=True)
vep_df['end'] = vep_df['start']

# A deletion is indicated by the exact nucleotide coordinates.
# https://www.ensembl.org/info/docs/tools/vep/vep_formats.html
vep_df['alt-ref'] = vep_df['alt'].apply(len) - vep_df['ref'].apply(len)
vep_df.loc[vep_df['alt-ref'] < 0, 'end'] = vep_df.loc[vep_df['alt-ref'] < 0, 'start'] - vep_df.loc[vep_df['alt-ref'] < 0, 'alt-ref']

# An insertion (of any size) is indicated by start coordinate = end coordinate + 1
# Looks like not necessary
# https://www.ensembl.org/info/docs/tools/vep/vep_formats.html
# vep_df.loc[vep_df['alt-ref'] > 0, 'start'] = vep_df['end'] + 1

#%%
vep_df['alleles'] = vep_df['ref'] + '/' + vep_df['alt']
vep_df['strand'] = '+'

#%% prepare vep input
out_col_lst = ['chrom', 'start', 'end', 'alleles', 'strand', 'rsid', 'gwas_category_count']
vep_df = vep_df[out_col_lst]
vep_df.sort_values(vep_df.columns.tolist(), inplace=True)

#%% write vep input
vep_df.to_csv(vep_input_path, sep="\t", index=False, header=False)

#%% run vep command
dir_cache = "/".join(vep_cache_info.split('/')[:-3])
cmd_str = "vep --offline --dir_cache {dir_cache} --cache_version 108 --plugin TSSDistance --force_overwrite -i {vep_input} " \
          "--output_file {vep_output}".format(vep_input=vep_input_path, dir_cache=dir_cache, vep_output=vep_output_path)
Logger.info(cmd_str)
output = subprocess.run(shlex.split(cmd_str), capture_output=True)
