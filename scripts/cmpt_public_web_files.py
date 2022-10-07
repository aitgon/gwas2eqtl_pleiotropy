from sqlalchemy import create_engine, select, asc
from gwas2eqtl_pleiotropy.Logger import Logger
from gwas2eqtl_pleiotropy.db import meta, coloc, gwas_annot_tbl, ensg2symbol_tbl, rsid2cytoband_tbl

import os
import pandas
import pathlib
import sys

# %%
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

# %% Insert into coloc
engine = create_engine("sqlite:///{}".format(db_sqlite), echo=False)
meta.create_all(engine)

# %%
with engine.connect() as con:
    for chrom in range(1, 23):  # loop over chroms to avoid memory problems
        # print(chrom)
        sel_cols = [coloc, gwas_annot_tbl.c.gwas_trait, gwas_annot_tbl.c.gwas_class, ensg2symbol_tbl.c.symbol,
                    rsid2cytoband_tbl.c.cytoband]
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
        col_name1_lst = []
        col_name1_lst.append('chrom')
        col_name1_lst.append('pos')
        col_name1_lst.append('cytoband')
        col_name1_lst.append('rsid')
        col_name1_lst.append('ref')
        col_name1_lst.append('alt')
        col_name1_lst.append('gwas_trait')
        col_name1_lst.append('gwas_class')
        col_name1_lst.append('gwas_beta')
        col_name1_lst.append('egene_symbol')
        col_name1_lst.append('eqtl_beta')
        col_name1_lst.append('eqtl_id')
        col_name1_lst.append('egene')
        col_name1_lst.append('gwas_id')
        col_name1_lst.append('gwas_pval')
        col_name1_lst.append('eqtl_pval')
        col_name1_lst.append('PP.H4.abf')
        col_name1_lst.append('SNP.PP.H4')
        col_name1_lst.append('PP.H3.abf')
        col_name1_lst.append('PP.H2.abf')
        col_name1_lst.append('PP.H1.abf')
        col_name1_lst.append('PP.H0.abf')
        col_name1_lst.append('coloc_region')
        col_name1_lst.append('coloc_lead_pos')
        col_name1_lst.append('coloc_lead_rsid')
        col_name1_lst.append('nsnps')
        coloc_df = coloc_df[col_name1_lst]
        # remove dots
        col_name2_lst = [col.replace('.', '_') for col in col_name1_lst]
        coloc_df.rename(dict(zip(col_name1_lst, col_name2_lst)), axis=1, inplace=True)
        if chrom == 1:  # write new
            coloc_df.to_csv(coloc_tsv_gz_path, sep='\t', index=False, header=True)
        elif chrom > 1:  # append
            coloc_df.to_csv(coloc_tsv_gz_path, sep='\t', index=False, header=False, mode='a')
