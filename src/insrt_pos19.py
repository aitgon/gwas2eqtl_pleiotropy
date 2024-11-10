import sqlalchemy

from gwas2eqtl_pleiotropy_db import Base, pos19
from sqlalchemy import create_engine
from liftover import get_lifter

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

#%% Create all tables
engine = create_engine(sa_url)
if sqlalchemy.inspect(engine).has_table(pos19.__tablename__):
    pos19.__table__.drop(engine)
Base.metadata.tables[pos19.__tablename__].create(bind=engine)

#%%
with engine.connect() as con:
    df = pandas.read_sql('coloc', con=con, columns=['chrom', 'pos']).drop_duplicates()
df.sort_values(df.columns.tolist(), inplace=True)
df.reset_index(inplace=True, drop=True)
df.index.rename('id', inplace=True)

converter = get_lifter('hg38', 'hg19')
df['chrom'] = df['chrom'].astype('str')
df['pos19'] = df.apply(lambda x: converter[x['chrom']][x['pos']], axis=1).explode().str[1]

with engine.connect() as connection:
    connection.execute((Base.metadata.tables[pos19.__tablename__]).delete())

# import pdb; pdb.set_trace()
df.to_sql('pos19', con=engine, if_exists='append', index=True)
