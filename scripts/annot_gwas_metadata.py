import numpy
from gwas2eqtl_pleiotropy.Logger import Logger
from gwas2eqtl_pleiotropy.URL import URL

import os
import pandas

import pathlib
import sys

# gwas_metadata_ods = "config/gwas418.ods"
# gwas_annot_tsv_path = "t.tsv"

#%%
help_cmd_str = "todo"
try:
    gwas_metadata_ods = sys.argv[1]
    gwas_annot_ods_path = sys.argv[2]
    if len(sys.argv) > 3:
        print("""Two many arguments!
        {}""".format(help_cmd_str))
        sys.exit(1)
except IndexError:
    print("""Argument missing!
    {}""".format(help_cmd_str))
    sys.exit(1)

# #%% Outdir
# if not '__file__' in locals():
#     __file__ = "dwnld_gwas_info.py"
outdir_path = os.path.dirname(gwas_annot_ods_path)
pathlib.Path(outdir_path).mkdir(parents=True, exist_ok=True)

#%%
metadata_df = pandas.read_excel(gwas_metadata_ods, engine='odf')

#%%
# import pdb; pdb.set_trace()
metadata_df.drop(['Unnamed: 4'], axis=1, inplace=True)

#%%
gwas_id_df = metadata_df[['id', 'trait', 'icd10_disease_level1', 'icd10_code_level1']].drop_duplicates()
gwas_id_df.rename({'id': 'gwas_id'}, axis=1, inplace=True)

#%%
gwas_class_df = metadata_df[["icd10_disease_level1.1", "icd10_code_level1.1", "count", "icd10_disease_level2", "icd10_code_level2", "class_pleiotropy"]].dropna(axis=0).drop_duplicates()

#%%
metadata2_df = gwas_id_df.merge(gwas_class_df, left_on=['icd10_disease_level1', 'icd10_code_level1'], right_on=['icd10_disease_level1.1', 'icd10_code_level1.1'], how='left')

if metadata2_df.loc[metadata2_df.isna().any(axis=1)].shape[0] > 0:
    Logger.error("Mapping error between gwas and classes")
    Logger.error(metadata2_df.loc[metadata2_df.isna().any(axis=1)])
    sys.exit(1)

#%%
url = "http://gwas-api.mrcieu.ac.uk/gwasinfo"
Logger.info("Downloading {}".format(url))
gwassinfo_json_path = URL(url).download()
metadata_df.rename({'manual_class': 'class'}, axis=1, inplace=True)

#%% no select
metadata_mrcieu_df = pandas.read_json(gwassinfo_json_path).T
metadata_mrcieu_df = metadata_mrcieu_df[['id', 'sample_size', 'ncontrol', 'ncase', 'consortium', 'pmid', 'year', 'author', 'nsnp']]
gwasinfo_tsv_path = os.path.join(outdir_path, "gwasinfo.tsv")
metadata_mrcieu_df.to_csv(gwasinfo_tsv_path, sep="\t", header=True, index=False)
metadata_mrcieu_df.rename({'id': 'gwas_id'}, axis=1, inplace=True)

#%%
gwas_annot_df = metadata2_df.merge(metadata_mrcieu_df, on='gwas_id')
# gwas_annot_df.to_csv(gwas_annot_tsv_path, sep="\t", header=True, index=False)
with pandas.ExcelWriter(gwas_annot_ods_path) as fout:
    gwas_annot_df.to_excel(fout, index=False)
