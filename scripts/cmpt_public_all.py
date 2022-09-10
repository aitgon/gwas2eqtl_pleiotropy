from gwas2eqtl_pleiotropy.EBIeQTLinfo import EBIeQTLinfo
from gwas2eqtl_pleiotropy.OpenGWASinfo import OpenGWASinfo
from gwas2eqtl_pleiotropy.PathManager import PathManager
from gwas2eqtl_pleiotropy.PyTabix import bgzip, tabix_index
from gwas2eqtl_pleiotropy.constants import coloc_raw_tsv_path, h4_cutoff

import os
import pandas
import pathlib
import sys


#%% Parameters
if not '__file__' in locals():
    __file__ = "cmpt_public_all.py"
if not os.path.isfile(coloc_raw_tsv_path):
    print("input file does not exit")
    sys.exit(1)

outdir_path = os.path.join(PathManager.get_project_path(), "out", os.path.basename(__file__))
pathlib.Path(outdir_path).mkdir(parents=True, exist_ok=True)

#%% Download eQTL annotations
eqtl_info_df = EBIeQTLinfo().df
eqtl_info_df.rename({'identifier': "eqtl_identifier", 'tissue_label': 'etissue_subcategory'}, axis=1, inplace=True)

#%% Download OpenGWAS annotations
open_gwas_df = OpenGWASinfo().df
open_gwas_df.rename({'subcategory': 'gwas_subcategory'}, axis=1, inplace=True)
open_gwas_df = open_gwas_df[['gwas_identifier', 'gwas_subcategory']].drop_duplicates()
# remove subcat with nan
open_gwas_df = open_gwas_df.loc[~open_gwas_df["gwas_subcategory"].isna(), ]

#%%
coloc_df = pandas.read_csv(coloc_raw_tsv_path, sep="\t")

#%%
coloc_df = coloc_df.merge(open_gwas_df, on='gwas_identifier')
coloc_df = coloc_df.merge(eqtl_info_df[['eqtl_identifier', 'etissue_subcategory']].drop_duplicates(), on='eqtl_identifier')

#%%
coloc_columns = ['chrom', 'pos', 'rsid', 'ref', 'alt',
                 'gwas_trait_name', 'gwas_subcategory', 'gwas_identifier', 'gwas_beta', 'gwas_pvalue',
                 'eqtl_identifier', 'etissue_subcategory', 'egene_symbol', 'egene', 'eqtl_beta', 'eqtl_pvalue',
                 'pp_h4', 'PP.H4.abf', 'nsnps', 'PP.H3.abf', 'PP.H2.abf', 'PP.H1.abf', 'PP.H0.abf', 'leadeqtl_pos', 'leadeqtl_egene']
coloc_df = coloc_df[coloc_columns]
coloc_df.rename({'chrom': '#chrom'}, axis=1, inplace=True)

#%%
coloc_df.sort_values(by=['#chrom', 'pos', "rsid"], inplace=True)
coloc_all_tsv_path = os.path.join(outdir_path, "coloc_all.tsv")
coloc_df.to_csv(coloc_all_tsv_path, index=False, sep="\t")

#%%
bgzip(coloc_all_tsv_path)
tabix_index(coloc_all_tsv_path + ".gz", chrom=1, start=2, end=2, comment="#")
coloc_df.to_csv(coloc_all_tsv_path, index=False, sep="\t")
