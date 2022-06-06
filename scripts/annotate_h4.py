from eqtl2gwas_pleiotropy.EBIeQTLinfo import EBIeQTLinfo
from eqtl2gwas_pleiotropy.OpenGWASinfo import OpenGWASinfo
from eqtl2gwas_pleiotropy.PathManager import PathManager
from eqtl2gwas_pleiotropy.constants import coloc_h4_tsv_path

import os
import pandas
import pathlib
import sys


#%% Parameters
# if not '__file__' in locals():
#     __file__ = "annotate_h4.py"
if not os.path.isfile(coloc_h4_tsv_path):
    print("input file does not exit")
    sys.exit(1)

outdir_path = os.path.join(PathManager.get_outdir_path(), os.path.basename(__file__))
pathlib.Path(outdir_path).mkdir(parents=True, exist_ok=True)

h4_annotated_tsv_path = os.path.join(outdir_path, 'h4_annotated.tsv')
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
coloc_h4_df = pandas.read_csv(coloc_h4_tsv_path, sep="\t")

#%%
coloc_h4_df = coloc_h4_df.merge(open_gwas_df, on='gwas_identifier')
coloc_h4_df = coloc_h4_df.merge(eqtl_info_df[['eqtl_identifier', 'etissue_subcategory']].drop_duplicates(), on='eqtl_identifier')

#%%
coloc_h4_df.sort_values(by=coloc_h4_df.columns.tolist(), inplace=True)
coloc_h4_df.to_csv(h4_annotated_tsv_path, sep="\t", index=False)
