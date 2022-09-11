from gwas2eqtl_pleiotropy.Logger import Logger
from gwas2eqtl_pleiotropy.UCSC import UCSC
from gwas2eqtl_pleiotropy.URL import URL
from gwas2eqtl_pleiotropy.db import meta, coloc
from gwas2eqtl_pleiotropy.db import ensg2symbol_tbl
from gwas2eqtl_pleiotropy.db import gwas_annot_tbl
from gwas2eqtl_pleiotropy.db import rsid2cytoband_tbl
from sqlalchemy import create_engine, select, Table

import os
import pandas
import pathlib
import sys
import mygene

#%%

help_cmd_str = "todo"
try:
    url_db = sys.argv[1]
    gwas_cat_ods_path = sys.argv[2]
    etissue_cat_ods_path = sys.argv[3]
    annotated_tsv_gz_path = sys.argv[4]
    if len(sys.argv) > 6:
        print("""Two many arguments!
        {}""".format(help_cmd_str))
        sys.exit(1)
except IndexError:
    print("""Argument missing!
    {}""".format(help_cmd_str))
    sys.exit(1)

pathlib.Path(os.path.dirname(annotated_tsv_gz_path)).mkdir(parents=True, exist_ok=True)
engine = create_engine(url_db, echo = False)

#%% Insert into coloc
meta.create_all(engine)

#%% Copy from colocimport
Logger.info("Copy from colocimport")
cols = ["chrom", "pos", "rsid", "ref", "alt", "egene", "gwas_beta", "gwas_pval", "gwas_id", "eqtl_beta", "eqtl_pval", "eqtl_id", "PP.H4.abf", "SNP.PP.H4", "nsnps", "PP.H3.abf", "PP.H2.abf", "PP.H1.abf", "PP.H0.abf", "coloc_lead_pos", "coloc_lead_rsid", "coloc_region"]
colocimport = Table('colocimport', meta, autoload_with=engine)
ins = coloc.insert().from_select(cols, colocimport.select())
with engine.connect() as con:
    con.execute(ins)

#%% cytoband
Logger.info("Annotate cytobands")
cytoband_url = "http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/cytoBand.txt.gz"
cytoband_path = URL(cytoband_url).download()
cyto_df = pandas.read_csv(cytoband_path, sep="\t", header=None, usecols=[0, 1, 2, 3])
cyto_df.columns = ['chrom', 'start', 'end', 'cytoband']
cyto_df['chrom'] = cyto_df['chrom'].str.replace('chr', '')

#%%
sel = select([coloc.c.chrom, coloc.c.pos, coloc.c.rsid]).distinct()
with engine.connect() as con:
    fetch_res = con.execute(sel).fetchall()
rsid2cytoband_df = pandas.DataFrame(fetch_res)

rsid2cytoband_df['chrom'] = rsid2cytoband_df['chrom'].astype(str)
cyto_df['chrom'] = cyto_df['chrom'].astype(str)
for rowi, row in cyto_df.iterrows():
    chrom = row['chrom']
    start = row['start']
    end = row['end']
    cytoband = row['cytoband']
    rsid2cytoband_df.loc[
        (rsid2cytoband_df['chrom'] == chrom) & (start <= rsid2cytoband_df['pos']) & (
                rsid2cytoband_df['pos'] <= end), 'cytoband'] = cytoband
rsid2cytoband_df['cytoband'] = rsid2cytoband_df['chrom'].astype(str) + rsid2cytoband_df['cytoband']
rsid2cytoband_df = rsid2cytoband_df[['rsid', 'cytoband']]

dlt = rsid2cytoband_tbl.delete()
ins = rsid2cytoband_tbl.insert()
with engine.connect() as con:
    con.execute(dlt)
    con.execute(ins, rsid2cytoband_df.to_dict('records'))

#%% gwas category
Logger.info("Annotate GWAS categories")
gwas_annot_df = pandas.read_excel(gwas_cat_ods_path, engine="odf")
gwas_annot_df = gwas_annot_df[['id', 'manual_category', 'trait']]
gwas_annot_df.rename({'id': "gwas_id", 'manual_category': 'gwas_category', 'trait': 'gwas_trait'}, axis=1, inplace=True)
gwas_annot_df.drop_duplicates(inplace=True)

dlt = gwas_annot_tbl.delete()
ins = gwas_annot_tbl.insert()
with engine.connect() as con:
    con.execute(dlt)
    con.execute(ins, gwas_annot_df.to_dict('records'))

#%% download gene symbols
Logger.info("Annotate egene symbols")
# ensg2symbol_df = UCSC().gene_id_to_symbol()

# egene_lst = None
sel = select(coloc.c.egene).distinct()
with engine.connect() as con:
    # egene_lst = [egene[0] for egene in con.execute(sel).fetchall()]
    fetch_res = con.execute(sel).fetchall()

ensg2symbol_df = pandas.DataFrame(fetch_res).drop_duplicates()
egene_lst = ensg2symbol_df['egene'].tolist()
mg = mygene.MyGeneInfo()
query_df = mg.getgenes(egene_lst, fields='symbol', as_dataframe=True)
query_df['egene'] = query_df.index
query_df = query_df[['egene', 'symbol']].drop_duplicates()
ensg2symbol_df = ensg2symbol_df.merge(query_df, left_on='egene', right_on='egene', how='left')
ensg2symbol_df.loc[ensg2symbol_df['symbol'].isna(), 'symbol'] = ensg2symbol_df.loc[ensg2symbol_df['symbol'].isna(), 'egene']
ensg2symbol_df.rename({'egene': 'gene_id'}, axis=1, inplace=True)
# import pdb; pdb.set_trace()
# ensg2symbol_df = ensg2symbol_df.merge(pandas.Series(egene_lst, name='gene_id'), left_on='gene_id', right_on='gene_id')

dlt = ensg2symbol_tbl.delete()
ins = ensg2symbol_tbl.insert()
with engine.connect() as con:
    con.execute(dlt)
    con.execute(ins, ensg2symbol_df.to_dict('records'))

# #%% delete table import to save space
Logger.info("Drop colocimport to save space")
dlt = colocimport.delete()
with engine.connect() as con:
    con.execute(dlt)

