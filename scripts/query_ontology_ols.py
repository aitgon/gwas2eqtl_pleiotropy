import urllib

import numpy
import requests

from gwas2eqtl_pleiotropy.Logger import Logger
from gwas2eqtl_pleiotropy.URL import URL

import os
import pandas

import pathlib
import sys

from gwas2eqtl_pleiotropy.constants import public_data_dir

gwas_metadata_ods = "/home/gonzalez/Repositories/gwas2eqtl_pleiotropy/config/gwas418.ods"

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
outdir_path = os.path.dirname(gwas_annot_ods_path)
pathlib.Path(outdir_path).mkdir(parents=True, exist_ok=True)

#%%
mygwas_df = pandas.read_excel(gwas_metadata_ods, engine='odf', usecols=['id', 'trait', 'ontology', 'query'])
mygwas_df['ontology_id'] = None
mygwas_df['ontology_label'] = None
mygwas_df['ontology_uri'] = None

for query in mygwas_df['query'].unique():
    ontology = mygwas_df.loc[mygwas_df['query'] == query, 'ontology'].unique()[0]
    ols_url = "https://www.ebi.ac.uk/ols/api/search?q={query}&ontology={ontology}".format(query=urllib.parse.quote(query), ontology=ontology)
    try:
        ols_df = pandas.read_json(ols_url)
    except:
        import pdb; pdb.set_trace()
    if len(ols_df['response']['docs']) > 0:
        if 'obo_id' in ols_df['response']['docs'][0]:
            print(query, ols_df['response']['docs'][0]['obo_id'], ols_df['response']['docs'][0]['label'])
            #print(trait, ols_df['response']['docs'][1]['obo_id'], ols_df['response']['docs'][1]['label'])
            # print()
            mygwas_df.loc[mygwas_df['query'] == query, 'ontology_id'] = ols_df['response']['docs'][0]['obo_id']
            mygwas_df.loc[mygwas_df['query'] == query, 'ontology_label'] = ols_df['response']['docs'][0]['label']
            mygwas_df.loc[mygwas_df['query'] == query, 'ontology_uri'] = ols_df['response']['docs'][0]['iri']
        elif len(ols_df['response']['docs']) > 0:
            if 'obo_id' in ols_df['response']['docs'][1]:
                print(query, ols_df['response']['docs'][1]['obo_id'], ols_df['response']['docs'][1]['label'])
                mygwas_df.loc[mygwas_df['query'] == query, 'ontology_id'] = ols_df['response']['docs'][1]['obo_id']
                mygwas_df.loc[mygwas_df['query'] == query, 'ontology_label'] = ols_df['response']['docs'][1]['label']
                mygwas_df.loc[mygwas_df['query'] == query, 'ontology_uri'] = ols_df['response']['docs'][1]['iri']

mygwas_df = mygwas_df[['id', 'trait', 'ontology_label', 'ontology_id', 'ontology_uri']]
mygwas_df.sort_values(['ontology_label', 'trait', 'id'], inplace=True)
with pandas.ExcelWriter(gwas_annot_ods_path) as fout:
    mygwas_df.to_excel(fout, index=False)
