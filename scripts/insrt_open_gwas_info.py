from gwas2eqtl_pleiotropy.db2 import Base
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
gwasinfo_json_url = "http://gwas-api.mrcieu.ac.uk/gwasinfo"

#%% no select
df = pandas.read_json(gwasinfo_json_url).T

#%% Create all tables
engine = create_engine(sa_url)
Base.metadata.create_all(engine)

# Delete and insert
engine.execute((Base.metadata.tables['open_gwas_info']).delete())

df.rename({'id': 'gwas_id'}, axis=1, inplace=True)
df.set_index('gwas_id', verify_integrity=True, inplace=True, drop=True)
# import pdb; pdb.set_trace()
df.to_sql('open_gwas_info', con=engine, if_exists='append', index=True)
