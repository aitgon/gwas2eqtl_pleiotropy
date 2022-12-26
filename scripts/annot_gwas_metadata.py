import requests

from gwas2eqtl_pleiotropy.Logger import Logger
from gwas2eqtl_pleiotropy.URL import URL
import urllib.request
import os
import pandas

import pathlib
import sys

from gwas2eqtl_pleiotropy.constants import public_data_dir

#%%
help_cmd_str = "todo"
try:
    config_ods = sys.argv[1]
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
outdir_path = os.path.dirname(gwas_annot_ods_path)
pathlib.Path(outdir_path).mkdir(parents=True, exist_ok=True)

#%%
config_df = pandas.read_excel(config_ods, engine='odf')

#%%
# metadata_df.drop(['Unnamed: 4'], axis=1, inplace=True)
#
# #%%
# gwas_id_df = metadata_df[['id', 'trait', 'icd10_disease_level1', 'icd10_code_level1']].drop_duplicates()
# gwas_id_df.rename({'id': 'gwas_id'}, axis=1, inplace=True)
#
# #%%
# gwas_class_df = metadata_df[["icd10_disease_level1.1", "icd10_code_level1.1", "count", "icd10_disease_level2", "icd10_code_level2", "class_pleiotropy"]].dropna(axis=0).drop_duplicates()
#
# #%%
# metadata2_df = gwas_id_df.merge(gwas_class_df, left_on=['icd10_disease_level1', 'icd10_code_level1'], right_on=['icd10_disease_level1.1', 'icd10_code_level1.1'], how='left')
#
# if metadata2_df.loc[metadata2_df.isna().any(axis=1)].shape[0] > 0:
#     Logger.error("Mapping error between gwas and classes")
#     Logger.error(metadata2_df.loc[metadata2_df.isna().any(axis=1)])
#     sys.exit(1)

#%%
url = "http://gwas-api.mrcieu.ac.uk/gwasinfo"
json_path = os.path.join(public_data_dir, url.replace('http://', ''))
pathlib.Path(os.path.dirname(json_path)).mkdir(exist_ok=True, parents=True)

# URL = "https://td-cdn.pw/api.php?download=tikdown.org-42500282235.mp4"
# FILE_TO_SAVE_AS = "myvideo.mp4" # the name you want to save file as
# resp = requests.get(url) # making requests to server
if not os.path.isfile(json_path):
    with open(json_path, "wb") as f: # opening a file handler to create new file
        f.write((requests.get(url)).content)  # writing content to file

# Logger.info("Downloading {}".format(url))
# gwassinfo_json_path = URL(url).download()
# metadata_df.rename({'manual_class': 'class'}, axis=1, inplace=True)

#%%
mrcieu_df = pandas.read_json(json_path).T
mrcieu_df = mrcieu_df[['id', 'sample_size', 'ncontrol', 'ncase', 'consortium', 'pmid', 'year', 'author', 'nsnp']]
gwasinfo_tsv_path = os.path.join(outdir_path, "gwasinfo.tsv")
mrcieu_df.to_csv(gwasinfo_tsv_path, sep="\t", header=True, index=False)
mrcieu_df.rename({'id': 'gwas_id'}, axis=1, inplace=True)

#%%
config_df.rename({'id': 'gwas_id'}, axis=1, inplace=True)
gwas_annot_df = config_df.merge(mrcieu_df, on='gwas_id')
with pandas.ExcelWriter(gwas_annot_ods_path) as fout:
    gwas_annot_df.to_excel(fout, index=False)
