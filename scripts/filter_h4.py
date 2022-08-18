from eqtl2gwas_pleiotropy.constants import h4_cutoff

import os
import pandas
import pathlib
import sys

#%%
help_cmd_str = "todo"
try:
    annotated_tsv_gz_path = sys.argv[1]
    h4_tsv_gz_path = sys.argv[2]
    if len(sys.argv) > 3:
        print("""Two many arguments!
        {}""".format(help_cmd_str))
        sys.exit(1)
except IndexError:
    print("""Argument missing!
    {}""".format(help_cmd_str))
    sys.exit(1)

pathlib.Path(os.path.dirname((h4_tsv_gz_path))).mkdir(parents=True, exist_ok=True)

#%% input
if not os.path.isfile(annotated_tsv_gz_path):
    print("input file does not exit")
    sys.exit(1)

#%%
coloc_df = pandas.read_csv(annotated_tsv_gz_path, sep="\t")
coloc_df.drop(['id'], inplace=True, axis=1)
coloc_df = coloc_df.loc[coloc_df['PP.H4.abf'] >= h4_cutoff, ]
coloc_df.to_csv(h4_tsv_gz_path, sep="\t", index=False)
