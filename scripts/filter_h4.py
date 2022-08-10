from eqtl2gwas_pleiotropy.PathManager import PathManager
from eqtl2gwas_pleiotropy.constants import coloc_raw_tsv_path, h4_cutoff, coloc_h4_tsv_path

import os
import pandas
import pathlib
import sys

# #%% output
# if not '__file__' in locals():
#     __file__ = "filter_h4.py"
# outdir_path = os.path.join(PathManager.get_project_path(), "out", os.path.basename(__file__))
# pathlib.Path(outdir_path).mkdir(parents=True, exist_ok=True)

#%%
help_cmd_str = "todo"
try:
    coloc_raw_tsv_path = sys.argv[1]
    coloc_h4_tsv_path = sys.argv[2]
    if len(sys.argv) > 3:
        print("""Two many arguments!
        {}""".format(help_cmd_str))
        sys.exit(1)
except IndexError:
    print("""Argument missing!
    {}""".format(help_cmd_str))
    sys.exit(1)
pathlib.Path(os.path.dirname((coloc_h4_tsv_path))).mkdir(parents=True, exist_ok=True)

#%% input
if not os.path.isfile(coloc_raw_tsv_path):
    print("input file does not exit")
    sys.exit(1)

#%%
coloc_df = pandas.read_csv(coloc_raw_tsv_path, sep="\t", compression='gzip')
coloc_df = coloc_df.loc[coloc_df['PP.H4.abf'] >= h4_cutoff, ]
coloc_df.to_csv(coloc_h4_tsv_path, sep="\t", index=False)
