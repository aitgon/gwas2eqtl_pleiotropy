from gwas2eqtl_pleiotropy.Logger import Logger
from psycopg2.extras import NumericRange
from gwas2eqtl_pleiotropy.db import Base
from sqlalchemy import create_engine

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

#%% cytoband
Logger.info("Insert cytobands")
cytoband_url = "http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/cytoBand.txt.gz"
cyto_df = pandas.read_csv(cytoband_url, sep="\t", header=None, usecols=[0, 1, 2, 3])
cyto_df.columns = ['chrom', 'start', 'end', 'cytoband']
cyto_df['chrom'] = cyto_df['chrom'].str.replace('chr', '')
cyto_df = cyto_df.loc[cyto_df['chrom'].isin([str(chrom) for chrom in range(1, 23)])]
cyto_df['chrom'] = cyto_df['chrom'].astype(int)
cyto_df.set_index(cyto_df['chrom'].astype(str) + cyto_df['cytoband'], inplace=True, verify_integrity=True)

engine.execute((Base.metadata.tables['cytoband']).delete())
cyto_df['start_end38'] = cyto_df.apply(lambda x: NumericRange(x['start']+1, x['end'], '[]'), axis=1)
cyto_df.drop(['start', 'end'], axis=1, inplace=True)
cyto_df.to_sql('cytoband', con=engine, if_exists='append', index=True, index_label='id')
