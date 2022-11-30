import os
import pathlib
import pandas
import sys
import sqlalchemy

from sqlalchemy.ext.automap import automap_base
from sqlalchemy import create_engine
from gwas2eqtl_pleiotropy.db2 import Base
from gwas2eqtl_pleiotropy.db2 import tophits


# %%
help_cmd_str = "todo"
try:
    url = sys.argv[1]
    gwas_metadata = sys.argv[2]
    tophits_path_strf = sys.argv[3]
    if len(sys.argv) > 4:
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

gwas_id_lst = pandas.read_excel(gwas_metadata, engine="odf", usecols=['id'])['id'].to_list()
for gwas_id in gwas_id_lst:
    tophits_path = tophits_path_strf.format(gwas_id=gwas_id)
    if pathlib.Path(tophits_path).stat().st_size == 0:
        continue
    df = pandas.read_csv(tophits_path, sep="\t", header=0)
    df.rename({'chr': 'chrom', 'position': 'pos', 'nea': 'ref', 'ea': 'alt', 'id': 'gwas_id', 'p': 'pval'}, inplace=True, axis=1)
    df['rsid'] = df['rsid'].str.replace('rs', '').astype(int)
    df.index = df['chrom'].astype(str) + "_" + df['pos'].astype(str) + "_" + df['ref'] + "_" + df['alt'] + "_" + df['gwas_id']
    tophits_tbl = tophits.__dict__['__table__']
    tophits_dlt_stmt = tophits_tbl.delete().where(tophits_tbl.c.gwas_id == gwas_id)
    with engine.connect() as con:
        con.execute(tophits_dlt_stmt)
    df.index.rename('id', inplace=True)
    # try:
    df.to_sql('tophits', con=engine, if_exists='append', index=True)
    # except sqlalchemy.exc.IntegrityError:
    #     import pdb; pdb.set_trace()

# for gwas_id in gwas_identifier_lst:
#     if gwas_counter % 10 == 0:
#         print("GWAS counter: " + str(gwas_counter))
#     for eqtl_id in eqtl_identifier_lst:
#         coloc_tsv_path = coloc_dir.format(**{'gwas_id': gwas_id, 'eqtl_id': eqtl_id})
#         if os.path.isfile(coloc_tsv_path):
#             df = pandas.read_csv(coloc_tsv_path, sep="\t", header=0)
#             df.rename({'egene': 'egene_id'}, axis=1, inplace=True)
#             if df.shape[0] > 0:
#                 df['rsid'] = df['rsid'].str.replace('rs', '').astype(int)
#                 df['coloc_lead_rsid'] = df['coloc_lead_rsid'].str.replace('rs', '').astype(int)
#                 df_index = df['chrom'].astype(str) + "_" + df['pos'].astype(str) + "_" + df['egene_id'] + "_" + df['gwas_id'] + "_" + df['eqtl_id'] + "_" + df['coloc_lead_pos'].astype(str)
#                 df.set_index(df_index, inplace=True, verify_integrity=True)
#                 df.index.rename('id', inplace=True)
#                 # Insert
#                 df.to_sql('coloc', con=engine, if_exists='append', index=True)
#             gwas_eqtl_counter = gwas_eqtl_counter + 1
#     gwas_counter = gwas_counter + 1

# Base = automap_base()
# engine = create_engine(engine_url)
# Base.prepare(autoload_with=engine)

# if sqlalchemy.inspect(engine).has_table('tophits'):
#     tophits.drop(engine)
# else:
#     Base.metadata.create_all(engine, tables=[tophits])
#
# for path in pathlib.Path(tophits_dir_strf_path).rglob('hg38.tsv'):
#     fp = path.resolve()
#     if os.path.getsize(fp) > 0:
#         # print(fp)
#         df = pandas.read_csv(fp, sep='\t')
#         df.rename({'chr': 'chrom', 'position': 'pos', 'p': 'pval', 'id': 'gwas_id'}, axis=1, inplace=True)
#         df.to_sql('tophits', engine_url, if_exists='append', index=False)
