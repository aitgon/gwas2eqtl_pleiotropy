import sqlalchemy

from gwas2eqtl_pleiotropy.Logger import Logger
from gwas2eqtl_pleiotropy.db import Base, opengwas2trait_ontology
from sqlalchemy import create_engine

import pandas
import sys


#%%
help_cmd_str = "todo"
try:
    url = sys.argv[1]
    gwas_config = sys.argv[2]
    if len(sys.argv) > 3:
        print("""Two many arguments!
        {}""".format(help_cmd_str))
        sys.exit(1)
except IndexError:
    print("""Argument missing!
    {}""".format(help_cmd_str))
    sys.exit(1)

#%% Create all tables
engine = create_engine(url)
if sqlalchemy.inspect(engine).has_table("opengwas2trait_ontology"):
    opengwas2trait_ontology.__table__.drop(engine)
Base.metadata.create_all(engine)

#%% gwas category
Logger.info("Annotate GWAS categories")

gwas_annot_df = pandas.read_excel(gwas_config, engine="odf")
gwas_annot_df = gwas_annot_df[['id', 'trait', 'ontology_id', 'ontology_term', 'ontology_iri']].drop_duplicates()
gwas_annot_df.rename({'id': 'gwas_id', 'trait': 'gwas_trait', 'ontology_term': 'gwas_ontology_term', 'ontology_id': 'gwas_ontology_id', 'ontology_iri': 'gwas_ontology_iri'}, inplace=True, axis=1)

gwas_annot_df.set_index('gwas_id', inplace=True, verify_integrity=True)
gwas_annot_df.to_sql('opengwas2trait_ontology', con=engine, if_exists='append', index=True, index_label='gwas_id')
