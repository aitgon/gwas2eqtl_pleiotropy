import math

from sqlalchemy import create_engine, MetaData, select, Table
from sqlalchemy.ext.automap import automap_base

from gwas2eqtl_pleiotropy.EBIeQTLinfo import EBIeQTLinfo
from gwas2eqtl_pleiotropy.Logger import Logger

import os
import pandas
import pathlib
import sys

from gwas2eqtl_pleiotropy.UCSC import UCSC
#%%
from gwas2eqtl_pleiotropy.URL import URL
from gwas2eqtl_pleiotropy.db import meta, coloc

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

# pathlib.Path(os.path.dirname(annotated_tsv_gz_path)).mkdir(parents=True, exist_ok=True)

#%% Insert into coloc
engine = create_engine(url_db, echo = False)
meta.create_all(engine)
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
for rowi, row in cyto_df.iterrows():
    chrom = row['chrom']
    if chrom in [str(c) for c in range(1, 23)]:
        start = row['start']
        end = row['end']
        cytoband = "{}".format(chrom) + row['cytoband']
        upd = coloc.update()
        upd = upd.where(coloc.c.chrom == chrom)
        upd = upd.where(coloc.c.pos >= start)
        upd = upd.where(coloc.c.pos <= end)
        upd = upd.values(cytoband=cytoband)
        with engine.connect() as con:
            con.execute(upd)

#%% download gene symbols
ensg2symbol_df = UCSC().gene_id_to_symbol()
ensg2symbol_df = ensg2symbol_df.set_index('gene_id')

egene_lst = None
sel = select(coloc.c.egene).distinct()
with engine.connect() as con:
    egene_lst = [egene[0] for egene in con.execute(sel).fetchall()]

with engine.connect() as con:
    for egene in egene_lst:
        symbol = ensg2symbol_df.loc[egene]['symbol']
        upd = coloc.update().where(coloc.c.egene == egene).values(egene_symbol=symbol)
        con.execute(upd)

#%% gwas category
gwas_annot_df = pandas.read_excel(gwas_cat_ods_path, engine="odf")
gwas_annot_df = gwas_annot_df[['id', 'manual_category', 'trait']]
gwas_annot_df.rename({'id': "gwas_id", 'manual_category': 'gwas_category', 'trait': 'gwas_trait'}, axis=1, inplace=True)

with engine.connect() as con:
    for gwas_id in gwas_annot_df['gwas_id'].unique():
        # print(gwas_id)
        gwas_category, gwas_trait = \
        gwas_annot_df.loc[gwas_annot_df['gwas_id'] == gwas_id, ['gwas_category', 'gwas_trait']].values[0]
        upd = coloc.update().where(coloc.c.gwas_id == gwas_id).values(gwas_category=gwas_category, gwas_trait=gwas_trait)
        con.execute(upd)


