import sqlalchemy

from gwas2eqtl_pleiotropy.Logger import Logger
from gwas2eqtl_pleiotropy.db2 import Base
from sqlalchemy import create_engine

import pandas
import sys


#%%

help_cmd_str = "todo"
try:
    url = sys.argv[1]
    gwas_metadata = sys.argv[2]
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
Base.metadata.create_all(engine)

#%% gwas category
Logger.info("Annotate GWAS categories")

gwas_metadata_df = pandas.read_excel(gwas_metadata, engine="odf")
gwas_id_trait_level1_df = gwas_metadata_df[['id', 'trait', 'icd10_code_level1']].drop_duplicates()
gwas_id_trait_level1_df.rename({'id': 'gwas_id', 'trait': 'gwas_trait'}, inplace=True, axis=1)
gwas_level1_class_df = gwas_metadata_df[['icd10_code_level1.1', 'class_pleiotropy']].drop_duplicates()
gwas_level1_class_df.rename({'icd10_code_level1.1': 'icd10_code_level1', 'class_pleiotropy': 'gwas_class'}, inplace=True, axis=1)
gwas_annot_df = gwas_id_trait_level1_df.merge(gwas_level1_class_df, on='icd10_code_level1').drop_duplicates(inplace=False)
gwas_annot_df.drop('icd10_code_level1', inplace=True, axis=1)
gwas_annot_df.set_index('gwas_id', inplace=True, verify_integrity=True)

engine.execute((Base.metadata.tables['gwas_annot']).delete())
gwas_annot_df.to_sql('gwas_annot', con=engine, if_exists='append', index=True, index_label='gwas_id')
