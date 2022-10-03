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
        coloc_cols = []
        coloc_cols.append('chrom')
        coloc_cols.append('pos')
        coloc_cols.append('cytoband')
        coloc_cols.append('rsid')
        coloc_cols.append('ref')
        coloc_cols.append('alt')
        coloc_cols.append('gwas_trait')
        coloc_cols.append('gwas_class')
        coloc_cols.append('gwas_beta')
        coloc_cols.append('egene_symbol')
        coloc_cols.append('eqtl_beta')
        coloc_cols.append('eqtl_id')
        coloc_cols.append('egene')
        coloc_cols.append('gwas_id')
        coloc_cols.append('gwas_pval')
        coloc_cols.append('eqtl_pval')
        coloc_cols.append('PP.H4.abf')
        coloc_cols.append('SNP.PP.H4')
        coloc_cols.append('PP.H3.abf')
        coloc_cols.append('PP.H2.abf')
        coloc_cols.append('PP.H1.abf')
        coloc_cols.append('PP.H0.abf')
        coloc_cols.append('coloc_region')
        coloc_cols.append('coloc_lead_pos')
        coloc_cols.append('coloc_lead_rsid')
        coloc_cols.append('nsnps')
        coloc_df = coloc_df[coloc_cols]
        if chrom == 1:  # write new
            coloc_df.to_csv(coloc_tsv_gz_path, sep='\t', index=False, header=True)
        elif chrom > 1:  # append
            coloc_df.to_csv(coloc_tsv_gz_path, sep='\t', index=False, header=False, mode='a')
