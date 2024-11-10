from psycopg2._range import NumericRange

from gwas2eqtl_pleiotropy_logger import Logger
from gwas2eqtl_pleiotropy_db import Base, mysql_ucsc_hg38_ncbirefseq
from sqlalchemy import create_engine
from sqlalchemy import text

import sqlalchemy
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

#%% Create table
engine = create_engine(sa_url)
if sqlalchemy.inspect(engine).has_table(mysql_ucsc_hg38_ncbirefseq.__tablename__):
    mysql_ucsc_hg38_ncbirefseq.__table__.drop(engine)
Base.metadata.tables[mysql_ucsc_hg38_ncbirefseq.__tablename__].create(bind=engine)

print("Query UCSC ncbiRefSeq table in hg38 database.")
# sql = "select n.chrom as chrom, n.txStart as start38, n.txEnd as end38, SUBSTRING_INDEX(n.name, '.', 1) as refseq_acc, n.strand as strand, n.name2 as symbol, SUBSTRING_INDEX(k.kgID, '.', 1) as gene_id from ncbiRefSeq n, kgXref k where (n.name like 'NM_%' or n.name like 'XM_%') and SUBSTRING_INDEX(n.name, '.', 1)=k.mRNA"
sql = "select rs.chrom, rs.txStart as refseq_transcript_start38, rs.txEnd as refseq_transcript_end38, rs.strand as refseq_transcript_strand, substring_index(rs.name, '.', 1) as refseq_transcript_id, kg.geneSymbol as symbol, substring_index(kg.kgID, '.', 1) as transcript_id, substring_index(kn.geneId, '.', 1) as gene_id from ncbiRefSeq rs, kgXref kg, knownAttrs kn where substring_index(rs.name, '.', 1)=kg.mRNA and kg.kgID=kn.kgID"
ncbirefseq_ucsc_mysql_url = "mariadb+mariadbconnector://genome:@genome-euro-mysql.soe.ucsc.edu/hg38"
engine_ucsc = sqlalchemy.create_engine(ncbirefseq_ucsc_mysql_url)
with engine_ucsc.begin() as conn:
    df = pandas.read_sql(sqlalchemy.text(sql), con=conn).drop_duplicates()

print("Insert UCSC ncbiRefSeq table")
df['chrom'] = df['chrom'].str.replace('chr', '')
chrom_lst = [str(chrom) for chrom in range(1, 23)] + ['X', 'Y']
df = df.loc[df['chrom'].isin(chrom_lst)]
df.loc[df['chrom'] == 'X', 'chrom'] = 23
df.loc[df['chrom'] == 'Y', 'chrom'] = 24
df['chrom'] = df['chrom'].astype(int)

df['refseq_transcript_start_end38'] = df.apply(
    lambda x: NumericRange(x['refseq_transcript_start38'] + 1, x['refseq_transcript_end38'], '[)'), axis=1)
# import pdb; pdb.set_trace()
# df.drop(['start', 'end'], axis=1, inplace=True)
df.to_sql(mysql_ucsc_hg38_ncbirefseq.__tablename__, con=engine, if_exists='append', index=True, index_label='id')
