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
mygwas_df = pandas.read_excel(gwas_metadata_ods, engine='odf', usecols=['id', 'trait', 'ontology', 'query'])
mygwas_df['ontology_id'] = None
mygwas_df['ontology_label'] = None
mygwas_df['ontology_uri'] = None

for query in mygwas_df['query'].unique():
    query_mask = mygwas_df['query'] == query
    ontology = mygwas_df.loc[query_mask, 'ontology'].unique()[0]
    ols_url = "https://www.ebi.ac.uk/ols/api/search?q={query}&ontology={ontology}".format(query=urllib.parse.quote(query), ontology=ontology)
    ols_df = pandas.read_json(ols_url)
    obo_id = None
    obo_label = None
    obo_iri = None
    if len(ols_df['response']['docs']) > 0:
        if 'obo_id' in ols_df['response']['docs'][0]:
            obo_id = ols_df['response']['docs'][0]['obo_id']
            obo_label = ols_df['response']['docs'][0]['label']
            obo_iri = ols_df['response']['docs'][0]['iri']
        elif len(ols_df['response']['docs']) > 1:
            if 'obo_id' in ols_df['response']['docs'][1]:
                obo_id = ols_df['response']['docs'][1]['obo_id']
                obo_label = ols_df['response']['docs'][1]['label']
                obo_iri = ols_df['response']['docs'][1]['iri']
    mygwas_df.loc[query_mask, 'ontology_id'] = obo_id
    mygwas_df.loc[query_mask, 'ontology_label'] = obo_label
    mygwas_df.loc[query_mask, 'ontology_uri'] = obo_iri
    if not obo_id is None:
        obo_id2 = obo_id.replace(':', '_')
        print(query, obo_id, obo_label, obo_id2 in pgs_cat_df['id'].tolist())
        if obo_id2 in pgs_cat_df['id'].tolist():
            cat_str = ';'.join(pgs_cat_df.loc[pgs_cat_df['id'] == obo_id2, 'category_label'].tolist())
            mygwas_df.loc[query_mask, 'category'] = cat_str

# import pdb; pdb.set_trace()

mygwas_df = mygwas_df[['id', 'trait', 'ontology_label', 'ontology_id', 'ontology_uri', 'category']]
mygwas_df.sort_values(['ontology_label', 'trait', 'id'], inplace=True)
with pandas.ExcelWriter(gwas_annot_ods_path) as fout:
    mygwas_df.to_excel(fout, index=False)
