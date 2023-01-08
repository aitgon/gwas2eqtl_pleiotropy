from gwas2eqtl_pleiotropy.db import Base, entrezgene2ensg2symbol
from sqlalchemy import create_engine

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

#%% Create all tables
engine = create_engine(sa_url)
if sqlalchemy.inspect(engine).has_table("entrezgene2ensg2symbol"):
    entrezgene2ensg2symbol.__table__.drop(engine)
Base.metadata.create_all(engine)

input_url = 'https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz'
df = pandas.read_csv(input_url, sep='\t', usecols=['#tax_id', 'GeneID', 'Symbol', 'dbXrefs'])
df = df.loc[df['#tax_id'] == 9606]
df = df.loc[df['dbXrefs'].str.contains('Ensembl:ENSG')]
df['ensg'] = df['dbXrefs'].str.extract(r'Ensembl:(ENSG\d+)')
df.drop(['dbXrefs', '#tax_id'], inplace=True, axis=1)
df.rename({'GeneID': 'entrezgene', 'Symbol': 'gene_symbol', 'ensg': 'gene_id'}, inplace=True, axis=1)

df.to_sql('entrezgene2ensg2symbol', con=engine, if_exists='append', index=False)
