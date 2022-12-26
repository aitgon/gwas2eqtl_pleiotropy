import os
import pathlib

import pandas
import requests
import sqlalchemy
import sys

from gwas2eqtl_pleiotropy.constants import public_data_dir
from gwas2eqtl_pleiotropy.db import Base, gwascatalog
from sqlalchemy import create_engine


#%%
help_cmd_str = "todo"
try:
    url = sys.argv[1]
    if len(sys.argv) > 2:
        print("""Two many arguments!
        {}""".format(help_cmd_str))
        sys.exit(1)
except IndexError:
    print("""Argument missing!
    {}""".format(help_cmd_str))
    sys.exit(1)

#%% Create all tables
engine = create_engine(url)
meta_data = sqlalchemy.MetaData(bind=engine)
if sqlalchemy.inspect(engine).has_table("gwascatalog"):
    gwascatalog.__table__.drop(engine)
Base.metadata.create_all(engine)

#%% Download
url = "http://ftp.ebi.ac.uk/pub/databases/gwas/releases/2022/12/22/gwas-catalog-associations_ontology-annotated.tsv"
tsv_path = os.path.join(public_data_dir, url.split('http://')[1])
if not os.path.isfile(tsv_path):
    pathlib.Path(os.path.dirname(tsv_path)).mkdir(exist_ok=True, parents=True)
    resp = requests.get(url) # making requests to server
    with open(tsv_path, "wb") as f: # opening a file handler to create new file
        f.write(resp.content) # writing content to file

df = pandas.read_csv(tsv_path, sep="\t", header=0, usecols=['PUBMEDID', 'STUDY', 'DISEASE/TRAIT', 'MAPPED_TRAIT', 'MAPPED_TRAIT_URI', 'STUDY ACCESSION'])
df.drop_duplicates(inplace=True)
df.rename(dict(zip(df.columns, [column.key for column in (gwascatalog.__table__.columns)][1:])), axis=1, inplace=True)
df.to_sql('gwascatalog', con=engine, if_exists='append', index=False)



