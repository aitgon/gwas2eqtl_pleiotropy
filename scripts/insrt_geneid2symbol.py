from gwas2eqtl_pleiotropy.Logger import Logger
from gwas2eqtl_pleiotropy.db import Base
from sqlalchemy import create_engine

import sqlalchemy
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

#%% download ensg gene ids to to gene symbols
Logger.info("Insert gene symbols")
sql = "select distinct SUBSTRING_INDEX(knownAttrs.geneId, '.', 1) as gene_id, kgXref.geneSymbol as symbol from knownAttrs, kgXref where knownAttrs.kgID=kgXref.kgID;"
user = "genome"
host = "genome-euro-mysql.soe.ucsc.edu"
db = "hg38"
ucsc_engine = sqlalchemy.create_engine("mariadb+mariadbconnector://{user}@{host}/{db}".format(user=user, host=host, db=db))
with ucsc_engine.connect() as ucsc_con:
    ensg2symbol_df = pandas.read_sql(sql=sql, con=ucsc_con, index_col='gene_id')
engine.execute((Base.metadata.tables['ensg2symbol']).delete())
ensg2symbol_df.to_sql('ensg2symbol', con=engine, if_exists='append', index=True)
