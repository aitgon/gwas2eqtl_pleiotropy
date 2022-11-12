import os
import pathlib

import pandas
import sys
from sqlalchemy import delete
import sqlalchemy
from sqlalchemy import MetaData
from sqlalchemy.ext.automap import automap_base
from sqlalchemy import create_engine

from pandas.io import sql
from gwas2eqtl_pleiotropy.db import tophits

# %%
help_cmd_str = "todo"
try:
    engine_url = sys.argv[1]
    tophits_dir_path = sys.argv[2]
    if len(sys.argv) > 3:
        print("""Two many arguments!
        {}""".format(help_cmd_str))
        sys.exit(1)
except IndexError:
    print("""Argument missing!
    {}""".format(help_cmd_str))
    sys.exit(1)


Base = automap_base()
engine = create_engine(engine_url)
Base.prepare(autoload_with=engine)

if sqlalchemy.inspect(engine).has_table('tophits'):
    tophits.drop(engine)
else:
    Base.metadata.create_all(engine, tables=[tophits])

for path in pathlib.Path(tophits_dir_path).rglob('hg38.tsv'):
    fp = path.resolve()
    if os.path.getsize(fp) > 0:
        # print(fp)
        df = pandas.read_csv(fp, sep='\t')
        df.rename({'chr': 'chrom', 'position': 'pos', 'p': 'pval', 'id': 'gwas_id'}, axis=1, inplace=True)
        df.to_sql('tophits', engine_url, if_exists='append', index=False)

