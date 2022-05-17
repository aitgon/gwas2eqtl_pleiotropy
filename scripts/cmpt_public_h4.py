from eqtl2gwas_pleiotropy.EBIeQTLinfo import EBIeQTLinfo
from eqtl2gwas_pleiotropy.OpenGWASinfo import OpenGWASinfo
from eqtl2gwas_pleiotropy.PathManager import PathManager
from eqtl2gwas_pleiotropy.PyTabix import bgzip, tabix_index
from eqtl2gwas_pleiotropy.constants import coloc_raw_tsv_path, h4_cutoff

import os
import pandas
import pathlib
import sys


#%% Parameters
if not '__file__' in locals():
    __file__ = "cmpt_public_h4.py"
if not os.path.isfile(coloc_raw_tsv_path):
    print("input file does not exit")
    sys.exit(1)

outdir_path = os.path.join(PathManager.get_outdir_path(), os.path.basename(__file__))
pathlib.Path(outdir_path).mkdir(parents=True, exist_ok=True)

#%% Input file
coloc_all_tsv_path = os.path.join(PathManager.get_outdir_path(), "cmpt_public_all.py", "coloc_all.tsv")

#%%
coloc_all_df = pandas.read_csv(coloc_all_tsv_path, sep="\t")

# %%
coloc_h4_df = coloc_all_df.loc[coloc_all_df['PP.H4.abf'] >= h4_cutoff, ]
coloc_h4_tsv_path = os.path.join(outdir_path, "coloc_h4_all.tsv")
coloc_h4_df.to_csv(coloc_h4_tsv_path, sep="\t", index=False)

#%%
bgzip(coloc_h4_tsv_path)
tabix_index(coloc_h4_tsv_path + ".gz", chrom=1, start=2, end=2, comment="#")
coloc_h4_df.to_csv(coloc_h4_tsv_path, index=False, sep="\t")
