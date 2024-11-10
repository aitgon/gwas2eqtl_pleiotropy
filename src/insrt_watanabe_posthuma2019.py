from psycopg2._range import NumericRange

from gwas2eqtl_pleiotropy_logger import Logger
from gwas2eqtl_pleiotropy_db import Base, mysql_ucsc_hg38_ncbirefseq, watanabe_posthuma2019
from sqlalchemy import create_engine
from sqlalchemy import text

import sqlalchemy
import pandas
import sys


#%%

help_cmd_str = "todo"
try:
    sa_url = sys.argv[1]
    input_path = sys.argv[2]
    if len(sys.argv) > 3:
        print("""Two many arguments!
        {}""".format(help_cmd_str))
        sys.exit(1)
except IndexError:
    print("""Argument missing!
    {}""".format(help_cmd_str))
    sys.exit(1)

# %% Create all tables
engine = create_engine(sa_url)
if sqlalchemy.inspect(engine).has_table(mysql_ucsc_hg38_ncbirefseq.__tablename__):
    mysql_ucsc_hg38_ncbirefseq.__table__.drop(engine)
Base.metadata.tables[mysql_ucsc_hg38_ncbirefseq.__tablename__].create(bind=engine)

# parser = argparse.ArgumentParser(description='Insert watanabe posthuma 2019 supplementary table ST12.')
# parser.add_argument('-u', '--url', help='SQLAlchemy URL.')
# parser.add_argument('-i', '--input', help='Input file Supplementary Tables 1–26.')
# args = parser.parse_args()
# sa_url = args.url
# input_path = args.input

# %% Create all tables
engine = create_engine(sa_url)
if sqlalchemy.inspect(engine).has_table(watanabe_posthuma2019.__tablename__):
    watanabe_posthuma2019.__table__.drop(engine)
Base.metadata.create_all(engine)

df = pandas.read_excel(input_path, sheet_name='ST 12', header=1, skiprows=0)
df[['chrom', 'pos19', 'ref', 'alt']] = df['uniqID'].str.split(':', expand=True)
df['chrom'] = df['chrom'].astype(int)
df['pos19'] = df['pos19'].astype(int)
df['rsid'] = df['rsID'].str.split('rs', expand=True)[1]
df.drop(['uniqID', 'rsID'], axis=1, inplace=True)
df.rename({'#traits': 'traits', '#domains': 'domains', 'Ear,_Nose,_Throat': 'Ear_Nose_Throat'}, axis=1, inplace=True)
df['hgvs_id'] = 'NC_' + df['chrom'].astype('str').str.zfill(6) + '.10:g.' + df['pos19'].astype('str') + df[
    'ref'] + '>' + df['alt']

df = df[[c.key for c in watanabe_posthuma2019.__table__.columns]]
df.set_index('hgvs_id', drop=True, inplace=True)
df.to_sql(watanabe_posthuma2019.__tablename__, con=engine, if_exists='append', index=True, index_label='hgvs_id')