from gwas2eqtl_pleiotropy.Logger import Logger
from gwas2eqtl_pleiotropy.UCSC import UCSC
from gwas2eqtl_pleiotropy.URL import URL
from gwas2eqtl_pleiotropy.db2 import Base
from sqlalchemy import create_engine

import pandas
import sys

#%%

help_cmd_str = "todo"
try:
    url = sys.argv[1]
    gwas_metadata = sys.argv[2]
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

#%% download ensg gene ids to to gene symbols
Logger.info("Annotate egene symbols")

ensg2symbol_df = UCSC(database='hg38').gene_id_to_symbol()
# Insert
ensg2symbol_df.to_sql('ensg2symbol', con=engine, if_exists='replace', index=True)

#%% cytoband
Logger.info("Annotate cytobands")
cytoband_url = "http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/cytoBand.txt.gz"
cytoband_path = URL(cytoband_url).download()
cyto_df = pandas.read_csv(cytoband_path, sep="\t", header=None, usecols=[0, 1, 2, 3])
cyto_df.columns = ['chrom', 'start', 'end', 'cytoband']
cyto_df['chrom'] = cyto_df['chrom'].str.replace('chr', '')
cyto_df = cyto_df.loc[cyto_df['chrom'].isin([str(chrom) for chrom in range(1, 23)])]
cyto_df['chrom'] = cyto_df['chrom'].astype(int)

cyto_df.set_index(cyto_df['chrom'].astype(str) + cyto_df['cytoband'], inplace=True, verify_integrity=True)
cyto_df.to_sql('cytoband', con=engine, if_exists='replace', index=True, index_label='id')

# #%%
# sel = select([coloc.c.chrom, coloc.c.pos, coloc.c.rsid]).distinct()
# with engine.connect() as con:
#     fetch_res = con.execute(sel).fetchall()
# rsid2cytoband_df = pandas.DataFrame(fetch_res)
# rsid2cytoband_df['chrom'] = rsid2cytoband_df['chrom'].astype(int)
# rsid2cytoband_df['pos'] = rsid2cytoband_df['pos'].astype(int)
#
# # rsid2cytoband_df['chrom'] = rsid2cytoband_df['chrom'].astype(str)
# # cyto_df['chrom'] = cyto_df['chrom'].astype(str)
# cyto_df = cyto_df.loc[cyto_df['chrom'].isin([str(chrom) for chrom in range(1, 23)])]
# cyto_df['chrom'] = cyto_df['chrom'].astype(int)
# for rowi, row in cyto_df.iterrows():
#     chrom = row['chrom']
#     start = row['start']
#     end = row['end']
#     cytoband = row['cytoband']
#     rsid2cytoband_df.loc[(rsid2cytoband_df['chrom'] == chrom) & (start <= rsid2cytoband_df['pos']) & (rsid2cytoband_df['pos'] <= end), 'cytoband'] = cytoband
# rsid2cytoband_df['cytoband'] = rsid2cytoband_df['chrom'].astype(str) + rsid2cytoband_df['cytoband']
# rsid2cytoband_df = rsid2cytoband_df[['rsid', 'cytoband']]
#
# dlt = rsid2cytoband_tbl.delete()
# ins = rsid2cytoband_tbl.insert()
# with engine.connect() as con:
#     con.execute(dlt)
#     con.execute(ins, rsid2cytoband_df.to_dict('records'))

#%% gwas category
Logger.info("Annotate GWAS categories")

gwas_metadata_df = pandas.read_excel(gwas_metadata, engine="odf")
gwas_id_trait_level1_df = gwas_metadata_df[['id', 'trait', 'icd10_code_level1']].drop_duplicates()
gwas_id_trait_level1_df.rename({'id': 'gwas_id', 'trait': 'gwas_trait'}, inplace=True, axis=1)
gwas_level1_class_df = gwas_metadata_df[['icd10_code_level1.1', 'class_pleiotropy']].drop_duplicates()
gwas_level1_class_df.rename({'icd10_code_level1.1': 'icd10_code_level1', 'class_pleiotropy': 'gwas_class'}, inplace=True, axis=1)
gwas_annot_df = gwas_id_trait_level1_df.merge(gwas_level1_class_df, on='icd10_code_level1').drop_duplicates(inplace=False)
gwas_annot_df.drop('icd10_code_level1', inplace=True, axis=1)
gwas_annot_df.set_index('gwas_id', inplace=True, verify_integrity=True)
gwas_annot_df.to_sql('gwas_annot', con=engine, if_exists='replace', index=True, index_label='gwas_id')
