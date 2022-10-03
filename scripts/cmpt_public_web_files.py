from sqlalchemy import create_engine, select, asc
from gwas2eqtl_pleiotropy.Logger import Logger
from gwas2eqtl_pleiotropy.db import meta, coloc, gwas_annot_tbl, ensg2symbol_tbl, rsid2cytoband_tbl

import os
import pandas
import pathlib
import sys


#%%
help_cmd_str = "todo"
try:
    db_sqlite = sys.argv[1]
    coloc_tsv_gz_path = sys.argv[2]
    if len(sys.argv) > 3:
        print("""Two many arguments!
        {}""".format(help_cmd_str))
        sys.exit(1)
except IndexError:
    print("""Argument missing!
    {}""".format(help_cmd_str))
    sys.exit(1)

outdir_path = os.path.dirname(coloc_tsv_gz_path)
pathlib.Path(outdir_path).mkdir(parents=True, exist_ok=True)

#%% Insert into coloc
engine = create_engine("sqlite:///{}".format(db_sqlite), echo = False)
meta.create_all(engine)

#%%
with engine.connect() as con:
    for chrom in range(1, 23):  # loop over chroms to avoid memory problems
        # print(chrom)
        sel_cols = [coloc, gwas_annot_tbl.c.gwas_trait, gwas_annot_tbl.c.gwas_class, ensg2symbol_tbl.c.symbol, rsid2cytoband_tbl.c.cytoband]
        sel = select(sel_cols).distinct()
        sel = sel.where(coloc.c.gwas_id == gwas_annot_tbl.c.gwas_id)
        sel = sel.where(coloc.c.egene == ensg2symbol_tbl.c.gene_id)
        sel = sel.where(coloc.c.rsid == rsid2cytoband_tbl.c.rsid)
        sel = sel.where(coloc.c.chrom == chrom)
        sel = sel.order_by(asc(coloc.c.pos))

        Logger.info("SQL query: chrom {}".format(chrom))
        result_lst = con.execute(sel).fetchall()
        if len(result_lst) == 0:
            continue
        coloc_df = pandas.DataFrame(result_lst)
        coloc_df.drop(['id'], inplace=True, axis=1)
        coloc_df.rename({'symbol': 'egene_symbol'}, axis=1, inplace=True)
        coloc_df = coloc_df[['chrom', 'pos', 'cytoband', 'rsid', 'ref', 'alt', 'gwas_id', 'gwas_trait', 'gwas_beta', 'gwas_pval', 'egene_symbol', 'egene', 'egene_symbol', 'eqtl_id', 'eqtl_beta', 'eqtl_pval', 'PP.H4.abf', 'SNP.PP.H4', 'PP.H3.abf', 'PP.H2.abf', 'PP.H1.abf', 'PP.H0.abf', 'nsnps', 'coloc_lead_pos', 'coloc_lead_rsid', 'coloc_region', 'gwas_class']]
        if chrom == 1:  # write new
            coloc_df.to_csv(coloc_tsv_gz_path, sep='\t', index=False, header=True)
        elif chrom > 1:  # append
            coloc_df.to_csv(coloc_tsv_gz_path, sep='\t', index=False, header=False, mode='a')
