import sqlalchemy

from gwas2eqtl_pleiotropy.Logger import Logger
from gwas2eqtl_pleiotropy.db import Base, gwas_annot
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
if sqlalchemy.inspect(engine).has_table("gwas_annot"):
    gwas_annot.__table__.drop(engine)
Base.metadata.create_all(engine)

#%% gwas category
Logger.info("Annotate GWAS categories")

gwas_annot_df = pandas.read_excel(gwas_config, engine="odf")
gwas_annot_df = gwas_annot_df[['id', 'trait', 'ontology_label', 'ontology_id', 'category_manual']].drop_duplicates()
gwas_annot_df.rename({'id': 'gwas_id', 'trait': 'gwas_trait', 'ontology_label': 'gwas_ontology_label', 'ontology_id': 'gwas_ontology_id', 'category_manual': 'gwas_category'}, inplace=True, axis=1)
# gwas_level1_class_df = gwas_config_df[['icd10_code_level1.1', 'class_pleiotropy']].drop_duplicates()

# gwas_level1_class_df.rename({'icd10_code_level1.1': 'icd10_code_level1', 'class_pleiotropy': 'gwas_class'}, inplace=True, axis=1)
# gwas_annot_df = gwas_id_trait_level1_df.merge(gwas_level1_class_df, on='icd10_code_level1').drop_duplicates(inplace=False)
# gwas_annot_df.drop('icd10_code_level1', inplace=True, axis=1)
gwas_annot_df.set_index('gwas_id', inplace=True, verify_integrity=True)
# import pdb; pdb.set_trace()
gwas_annot_df.to_sql('gwas_annot', con=engine, if_exists='append', index=True, index_label='gwas_id')
