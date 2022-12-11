from gwas2eqtl_pleiotropy.Logger import Logger
from psycopg2.extras import NumericRange
from gwas2eqtl_pleiotropy.db2 import Base
from sqlalchemy import create_engine

import pandas
import sys

#%%

help_cmd_str = "todo"
try:
    url = sys.argv[1]
    etissue_cat_ods_path = sys.argv[2]
    if len(sys.argv) > 3:
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

eqtl_annot_df = pandas.read_excel(etissue_cat_ods_path, index_col=None, header=0)
eqtl_annot_df.drop(['Unnamed: 8', 'etissue_class.1', 'count'], axis=1, inplace=True)

eqtl_url = "https://raw.githubusercontent.com/eQTL-Catalogue/eQTL-Catalogue-resources/master/tabix/tabix_ftp_paths.tsv"
eqtl_df = pandas.read_csv(eqtl_url, sep="\t", header=0)
eqtl_df = eqtl_df.loc[eqtl_df['ftp_path'].str.contains('/ge/|/microarray/',regex=True,na=False), ]
eqtl_df['index'] = (eqtl_df['ftp_path'].str.replace('.all.tsv.gz', '', regex=True)).str.split('/', expand=True)[10].tolist()
eqtl_df.set_index('index', drop=True, verify_integrity=True, inplace=True)

eqtl_annot_df = eqtl_df.merge(eqtl_annot_df, on=['study', 'qtl_group', 'tissue_ontology_id', 'tissue_ontology_term', 'tissue_label', 'condition_label'])
eqtl_annot_df.set_index(eqtl_df.index, verify_integrity=True, inplace=True)
eqtl_annot_df.index.rename('eqtl_id', inplace=True)

Logger.info("Insert eqtl_annot")
# import pdb; pdb.set_trace()
engine.execute((Base.metadata.tables['eqtl_annot']).delete())
eqtl_annot_df.to_sql('eqtl_annot', con=engine, if_exists='append', index=True)
