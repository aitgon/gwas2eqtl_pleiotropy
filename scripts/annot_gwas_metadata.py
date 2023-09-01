import os
import pandas
import pathlib
import requests
import sys

from gwas2eqtl_pleiotropy.constants import public_data_dir

#%%
help_cmd_str = "todo"
try:
    gwas_trait_ods = sys.argv[1]
    gwas_category_ods = sys.argv[2]
    gwas_annot_ods_path = sys.argv[3]
    if len(sys.argv) > 4:
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
gwas_trait_df = pandas.read_excel(gwas_trait_ods, engine='odf')
gwas_category_df = pandas.read_excel(gwas_category_ods, engine='odf')

#%%
url = "http://gwas-api.mrcieu.ac.uk/gwasinfo"
json_path = os.path.join(public_data_dir, url.replace('http://', ''))
pathlib.Path(os.path.dirname(json_path)).mkdir(exist_ok=True, parents=True)
if not os.path.isfile(json_path):
    with open(json_path, "wb") as f:  # opening a file handler to create new file
        f.write((requests.get(url)).content)  # writing content to file


#%%
mrcieu_df = pandas.read_json(json_path).T
mrcieu_df = mrcieu_df[['id', 'sample_size', 'ncontrol', 'ncase', 'consortium', 'pmid', 'year', 'author', 'nsnp']]
gwasinfo_tsv_path = os.path.join(outdir_path, "gwasinfo.tsv")
mrcieu_df.to_csv(gwasinfo_tsv_path, sep="\t", header=True, index=False)
mrcieu_df.rename({'id': 'gwas_id'}, axis=1, inplace=True)

#%%
gwas_trait_df.drop(['query_ontology', 'query_term', 'ontology_iri', 'category'], axis=1, inplace=True)
gwas_trait_df.rename({'id': 'gwas_id', 'ontology_id': 'gwas_trait_ontology_id', 'ontology_term': 'gwas_trait_ontology_term'}, axis=1, inplace=True)
gwas_category_df.drop(['query_ontology', 'query_term', 'ontology_iri', 'category'], axis=1, inplace=True)
gwas_category_df.rename({'id': 'gwas_id', 'ontology_id': 'gwas_category_ontology_id', 'ontology_term': 'gwas_category_ontology_term'}, axis=1, inplace=True)
gwas_annot_df = gwas_trait_df.merge(gwas_category_df, on=['gwas_id', 'trait'])
gwas_annot_df = gwas_annot_df.merge(mrcieu_df, on=['gwas_id'])
with pandas.ExcelWriter(gwas_annot_ods_path) as fout:
    gwas_annot_df.to_excel(fout, index=False)
