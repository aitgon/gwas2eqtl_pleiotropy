"""Reads ODS with "id", "query_ontology", "query_term" columns
and queries OLS to add columns "ontology_id", "ontology_term", "ontology_iri"
I use it to normalize GWAS trait names and categories
"""

from gwas2eqtl_pleiotropy.Logger import Logger

import os
import pandas
import pathlib
import sys
import urllib


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

pgs_all_metadata_json = 'https://www.pgscatalog.org/rest/trait_category/all'
pgs_all_metadata_df = pandas.read_json(pgs_all_metadata_json)
pgs_cat_df = pandas.DataFrame()
for rowi, row in pgs_all_metadata_df.iterrows():
    category_label = pgs_all_metadata_df.loc[rowi, 'results']['label']
    category_label_df = pandas.DataFrame(pgs_all_metadata_df.loc[rowi, 'results']['efotraits'])
    category_label_df['category_label'] = category_label
    pgs_cat_df = pandas.concat([pgs_cat_df, category_label_df], axis=0)

#%%
mygwas_df = pandas.read_excel(gwas_metadata_ods, engine='odf', usecols=['id', 'trait', 'query_ontology', 'query_term'])
mygwas_df['ontology_id'] = None
mygwas_df['ontology_term'] = None
mygwas_df['ontology_iri'] = None

for query_term in sorted(mygwas_df['query_term'].unique()):
    query_mask = mygwas_df['query_term'] == query_term
    query_ontology = mygwas_df.loc[query_mask, 'query_ontology'].unique()[0]
    ols_url = "https://www.ebi.ac.uk/ols/api/search?q={query_term}&query_ontology={query_ontology}".format(query_term=urllib.parse.quote(query_term), query_ontology=query_ontology)
    ols_df = pandas.read_json(ols_url)
    obo_id = None
    obo_label = None
    obo_iri = None
    if len(ols_df['response']['docs']) > 0:
        for doc in ols_df['response']['docs']:
            if doc['ontology_name'] == query_ontology:
                if 'obo_id' in doc:
                    obo_id = doc['obo_id']
                elif 'id' in doc:
                    obo_id = doc['id']
                obo_label = doc['label']
                obo_iri = doc['iri']
                break
    if obo_id.startswith("snomed:class"):
        obo_id = obo_id.split('/')[-1]
    mygwas_df.loc[query_mask, 'ontology_id'] = obo_id
    mygwas_df.loc[query_mask, 'ontology_term'] = obo_label.lower()
    mygwas_df.loc[query_mask, 'ontology_iri'] = obo_iri
    if not obo_id is None:
        obo_id2 = obo_id.replace(':', '_')
        Logger.info("{} | {} | {} | {}".format(query_term, obo_label, obo_id, obo_iri))
        if obo_id2 in pgs_cat_df['id'].tolist():
            cat_str = ';'.join(pgs_cat_df.loc[pgs_cat_df['id'] == obo_id2, 'category_label'].tolist())
            mygwas_df.loc[query_mask, 'category'] = cat_str

# mygwas_df = mygwas_df[['id', 'trait', 'query_ontology', 'query_term', 'ontology_id', 'ontology_term', 'ontology_iri']]
mygwas_df.sort_values('id', inplace=True)
with pandas.ExcelWriter(gwas_annot_ods_path) as fout:
    mygwas_df.to_excel(fout, index=False)
