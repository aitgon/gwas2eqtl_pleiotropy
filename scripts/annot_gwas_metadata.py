from gwas2eqtl_pleiotropy.Logger import Logger
from gwas2eqtl_pleiotropy.URL import URL

import os
import pandas
import pathlib
import sys


#%%
help_cmd_str = "todo"
try:
    gwas_nonannot_ods = sys.argv[1]
    gwas_annot_tsv_path = sys.argv[2]
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
outdir_path = os.path.dirname(gwas_annot_tsv_path)
pathlib.Path(outdir_path).mkdir(parents=True, exist_ok=True)

#%%
gwas_nonannot_df = pandas.read_excel(gwas_nonannot_ods, engine='odf')
gwas_nonannot_df.drop(['Unnamed: 4', 'Unnamed: 5', 'manual_category.1', 'count'], axis=1, inplace=True)
gwas_nonannot_df = gwas_nonannot_df[['id', 'trait', 'manual_category']]

#%%
url = "http://gwas-api.mrcieu.ac.uk/gwasinfo"
Logger.info("Downloading {}".format(url))
gwassinfo_json_path = URL(url).download()
gwas_nonannot_df.rename({'manual_category': 'category'}, axis=1, inplace=True)

#%% no select
gwasinfo_df = pandas.read_json(gwassinfo_json_path).T
gwasinfo_df = gwasinfo_df[['id', 'sample_size', 'ncontrol', 'ncase', 'consortium', 'pmid', 'year', 'author', 'nsnp']]
gwasinfo_df.to_csv(gwas_annot_tsv_path, sep="\t", header=True, index=False)
gwasinfo_tsv_path = os.path.join(outdir_path, "gwasinfo.tsv")
gwasinfo_df.to_csv(gwasinfo_tsv_path, sep="\t", header=True, index=False)

#%%
gwas_annot_df = gwas_nonannot_df.merge(gwasinfo_df, on='id')
gwas_annot_df.to_csv(gwas_annot_tsv_path, sep="\t", header=True, index=False)
