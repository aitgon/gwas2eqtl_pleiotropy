import os
import pathlib

import requests
import sqlalchemy

from gwas2eqtl_pleiotropy.constants import public_data_dir
from gwas2eqtl_pleiotropy.db import Base, open_gwas_info
from sqlalchemy import create_engine

import pandas
import sys


#%%
help_cmd_str = "todo"
try:
    sa_url = sys.argv[1]
    if len(sys.argv) > 2:
        print("""Two many arguments!
        {}""".format(help_cmd_str))
        sys.exit(1)
except IndexError:
    print("""Argument missing!
    {}""".format(help_cmd_str))
    sys.exit(1)

#%%
url = "http://gwas-api.mrcieu.ac.uk/gwasinfo"
json_path = os.path.join(public_data_dir, url.replace('http://', ''))
pathlib.Path(os.path.dirname(json_path)).mkdir(exist_ok=True, parents=True)
if not os.path.isfile(json_path):
    with open(json_path, "wb") as f:  # opening a file handler to create new file
        f.write((requests.get(url)).content)  # writing content to file

#%% no select
df = pandas.read_json(json_path).T

#%% Create all tables
engine = create_engine(sa_url)
if sqlalchemy.inspect(engine).has_table("open_gwas_info"):
    open_gwas_info.__table__.drop(engine)
Base.metadata.create_all(engine)

# Delete and insert
engine.execute((Base.metadata.tables['open_gwas_info']).delete())

df.rename({'id': 'gwas_id'}, axis=1, inplace=True)
df.set_index('gwas_id', verify_integrity=True, inplace=True, drop=True)
df.to_sql('open_gwas_info', con=engine, if_exists='append', index=True)
