from gwas2eqtl_pleiotropy.Logger import Logger
from psycopg2.extras import NumericRange
from gwas2eqtl_pleiotropy.db import Base
from sqlalchemy import create_engine
from liftover import get_lifter

import pandas
import sys

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
Base.metadata.create_all(engine)

#%%
df = pandas.read_sql('coloc', con=engine, columns=['chrom', 'pos']).drop_duplicates()
df.sort_values(df.columns.tolist(), inplace=True)
df.reset_index(inplace=True, drop=True)
df.index.rename('id', inplace=True)

converter = get_lifter('hg38', 'hg19')
df['pos19'] = df.apply(lambda x: converter[x['chrom']][x['pos']], axis=1).explode().str[1]

engine.execute((Base.metadata.tables['pos19']).delete())
# import pdb; pdb.set_trace()
df.to_sql('pos19', con=engine, if_exists='append', index=True)
