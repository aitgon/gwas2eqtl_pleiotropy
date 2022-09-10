from gwas2eqtl_pleiotropy.EBIeQTLinfo import EBIeQTLinfo
from gwas2eqtl_pleiotropy.constants import h4_cutoff

import os
import pandas
import pathlib
import sys

#%%
help_cmd_str = "todo"
try:
    pp_h4_abf = float(sys.argv[1])
    annotated_tsv_gz_path = sys.argv[2]
    etissue_cat_ods_path = sys.argv[3]
    web_tsv_gz_path = sys.argv[4]
    if len(sys.argv) > 5:
        print("""Two many arguments!
        {}""".format(help_cmd_str))
        sys.exit(1)
except IndexError:
    print("""Argument missing!
    {}""".format(help_cmd_str))
    sys.exit(1)

outdir_path = os.path.dirname(web_tsv_gz_path)
pathlib.Path(outdir_path).mkdir(parents=True, exist_ok=True)

#%% input
if not os.path.isfile(annotated_tsv_gz_path):
    print("input file does not exit")
    sys.exit(1)

#%%
coloc_df = pandas.read_csv(annotated_tsv_gz_path, sep="\t")
coloc_df.drop(['id'], inplace=True, axis=1)
coloc_df = coloc_df.loc[coloc_df['PP.H4.abf'] >= pp_h4_abf, ]

#%%
eqtl_info_df = EBIeQTLinfo().df
etissue_cat_df = pandas.read_excel(etissue_cat_ods_path, engine="odf")
etissue_cat_df.drop(['Unnamed: 8', 'etissue_category.1', 'count'], axis=1, inplace=True)
eqtl_info_df = eqtl_info_df.merge(etissue_cat_df, on=['study', 'qtl_group', 'tissue_ontology_id', 'tissue_ontology_term', 'tissue_label', 'condition_label'])
eqtl_info_df.rename({'identifier': "eqtl_identifier", 'ref': "eqtl_ref"}, axis=1, inplace=True)

#%%
coloc_df = coloc_df.merge(eqtl_info_df[['eqtl_identifier', 'eqtl_ref']].drop_duplicates(), on='eqtl_identifier')

#%%
# import pdb; pdb.set_trace()
# coloc_df.drop(['id'], inplace=True, axis=1)
coloc_df = coloc_df.loc[coloc_df['PP.H4.abf'] >= pp_h4_abf, ]

#%%
coloc_df.to_csv(web_tsv_gz_path, sep="\t", index=True, index_label='id')
