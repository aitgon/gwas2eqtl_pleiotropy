import sqlalchemy
from sqlalchemy import create_engine
from gwas2eqtl_pleiotropy.db2 import Base

import pandas
import sys
import os

#%%
help_cmd_str = "todo"
try:
    url = sys.argv[1]
    gwas_identifier_ods_path = sys.argv[2]
    eqtl_identifier_tsv_path = sys.argv[3]
    coloc_dir = sys.argv[4]
    if len(sys.argv) > 5:
        print("""Two many arguments!
        {}""".format(help_cmd_str))
        sys.exit(1)
except IndexError:
    print("""Argument missing!
    {}""".format(help_cmd_str))
    sys.exit(1)

# %% Create all tables
engine = create_engine(url)
Base.metadata.create_all(engine)

###############################################################################
#
# select immune cell types
#
###############################################################################

#%%
eqtl_df = pandas.read_csv(eqtl_identifier_tsv_path, sep="\t")
eqtl_df = eqtl_df.loc[eqtl_df['ftp_path'].str.contains('/ge/|/microarray/', regex=True, na=False), ]
eqtl_identifier_lst = sorted((eqtl_df['ftp_path'].str.replace('.all.tsv.gz', '', regex=True)).str.split('/', expand=True)[10].tolist())

#%%
gwas_identifier_df = pandas.read_excel(gwas_identifier_ods_path, header=0)
gwas_identifier_lst = sorted(gwas_identifier_df["id"].tolist())

# Delete data if exists
coloc_tab = Base.metadata.tables['coloc']
stmt = coloc_tab.delete()
engine.execute(stmt)

#%%
# with open(tsv_path, 'w') as fout:
gwas_eqtl_counter = 0
gwas_counter = 0
for gwas_id in gwas_identifier_lst:
    if gwas_counter % 10 == 0:
        print("GWAS counter: " + str(gwas_counter))
    for eqtl_id in eqtl_identifier_lst:
        coloc_tsv_path = coloc_dir.format(**{'gwas_id': gwas_id, 'eqtl_id': eqtl_id})
        if os.path.isfile(coloc_tsv_path):
            df = pandas.read_csv(coloc_tsv_path, sep="\t", header=0)
            df.rename({'egene': 'egene_id'}, axis=1, inplace=True)
            if df.shape[0] > 0:
                df['rsid'] = df['rsid'].str.replace('rs', '').astype(int)
                df['coloc_lead_rsid'] = df['coloc_lead_rsid'].str.replace('rs', '').astype(int)
                df_index = df['chrom'].astype(str) + "_" + df['pos'].astype(str) + "_" + df['egene_id'] + "_" + df['gwas_id'] + "_" + df['eqtl_id'] + "_" + df['coloc_lead_pos'].astype(str)
                df.set_index(df_index, inplace=True, verify_integrity=True)
                df.index.rename('id', inplace=True)
                # Insert
                df.to_sql('coloc', con=engine, if_exists='append', index=True)
            gwas_eqtl_counter = gwas_eqtl_counter + 1
    gwas_counter = gwas_counter + 1
