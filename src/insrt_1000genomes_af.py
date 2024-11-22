import glob
import os
from collections import defaultdict

from gwas2eqtl_pleiotropy_db import Base, af_1000genomes
from sqlalchemy import create_engine

import sqlalchemy
import pandas
import sys


#%%

help_cmd_str = "todo"
try:
    sa_url = sys.argv[1]
    data_dir = sys.argv[2]
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
if sqlalchemy.inspect(engine).has_table(af_1000genomes.__tablename__):
    af_1000genomes.__table__.drop(engine)
Base.metadata.tables[af_1000genomes.__tablename__].create(bind=engine)

id = 0  # primary key
print(id)
# import pdb; pdb.set_trace()
for tsv_gz in sorted(glob.glob(os.path.join(data_dir, "ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502", "ALL.chr*.tsv.gz"))):
    print(tsv_gz)

    df = pandas.read_csv(tsv_gz, sep="\t", header=None,
                         names=['chrom', 'pos19', 'rsid', 'ref', 'alt', 'eas_af', 'amr_af', 'afr_af', 'eur_af',
                                'sas_af'],
                         dtype='str')
    # df.columns = ['chrom', 'pos19', 'variantid', 'ref', 'alt', 'maf']

    # Explode list alternative alleles
    df['alt'] = df['alt'].str.split(',')
    df['eas_af'] = df['eas_af'].str.split(',')
    df['amr_af'] = df['amr_af'].str.split(',')
    df['afr_af'] = df['afr_af'].str.split(',')
    df['eur_af'] = df['eur_af'].str.split(',')
    df['sas_af'] = df['sas_af'].str.split(',')

    df = df.set_index(['chrom', 'pos19', 'rsid', 'ref']).apply(pandas.Series.explode).reset_index()
    df = df.query('sas_af!="."')

    df['eas_af'] = df['eas_af'].astype(float)
    df['amr_af'] = df['amr_af'].astype(float)
    df['afr_af'] = df['afr_af'].astype(float)
    df['eur_af'] = df['eur_af'].astype(float)
    df['sas_af'] = df['sas_af'].astype(float)

    # Explode list of variant ids
    df['rsid'] = df['rsid'].str.split(';')
    df = df.explode('rsid')
    df = df.loc[df['rsid'].str.startswith('rs')]  # keep only rsid, not structural variants

    # Put on side duplicated and remove
    dupli_df = df.loc[df.duplicated(keep=False)].sort_values(df.columns.tolist())
    df = df.drop_duplicates(keep=False)

    df['rsid'] = df['rsid'].str.split('rs', expand=True)[1].astype(int)
    df.loc[df['chrom'] == 'X', 'chrom'] = "23"

    df['variant_id'] = 'NC_' + df['chrom'].astype(str).str.zfill(6) + ".10:g." + df['pos19'].astype(str) + df[
        'ref'] + ">" + df['alt']

    # Insert
    df['id'] = [*range(id, id + df.shape[0])]
    df.set_index('id', drop=True, inplace=True)
    df.to_sql(af_1000genomes.__tablename__, con=engine, if_exists='append', index=True)
    id = id + df.shape[0]
