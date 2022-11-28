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
    coloc5_tsv_gz_path = sys.argv[3]
    coloc8_tsv_gz_path = sys.argv[4]
    coloc95_tsv_gz_path = sys.argv[5]
    if len(sys.argv) > 6:
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
pathlib.Path(coloc5_tsv_gz_path).unlink(missing_ok=True)
pathlib.Path(coloc8_tsv_gz_path).unlink(missing_ok=True)
pathlib.Path(coloc95_tsv_gz_path).unlink(missing_ok=True)

with engine.connect() as con:
    for chrom in range(1, 23):  # loop over chroms to avoid memory problems
        for start_mb in range(250):  # loop over positions in 1mb position windows to avoid memory problems
            start = start_mb*1e6
            end = start+1e6-1
            sql_sel_cols = [coloc, gwas_annot_tbl.c.gwas_trait, gwas_annot_tbl.c.gwas_class, ensg2symbol_tbl.c.symbol,
                            rsid2cytoband_tbl.c.cytoband]
            sql_sel = select(sql_sel_cols).distinct()
            sql_sel = sql_sel.where(coloc.c.gwas_id == gwas_annot_tbl.c.gwas_id)
            sql_sel = sql_sel.where(coloc.c.egene == ensg2symbol_tbl.c.gene_id)
            sql_sel = sql_sel.where(coloc.c.rsid == rsid2cytoband_tbl.c.rsid)
            sql_sel = sql_sel.where(coloc.c.chrom == chrom)
            sql_sel = sql_sel.where(coloc.c.pos >= start)
            sql_sel = sql_sel.where(coloc.c.pos <= end)
            sql_sel = sql_sel.order_by(asc(coloc.c.pos))

            Logger.info("SQL query: chrom {} start {} end {}".format(chrom, start, end))
            result_lst = con.execute(sql_sel).fetchall()
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
            # if chrom == 1 and start == 0:  # write new
            #     print("if")
            #     coloc_df.to_csv(coloc_tsv_gz_path, sep='\t', index=False, header=True)
            #     (coloc_df.loc[coloc_df['PP_H4_abf'] >= 0.5]).to_csv(coloc5_tsv_gz_path, sep='\t', index=False, header=True)
            #     (coloc_df.loc[coloc_df['PP_H4_abf'] >= 0.8]).to_csv(coloc8_tsv_gz_path, sep='\t', index=False, header=True)
            #     (coloc_df.loc[coloc_df['PP_H4_abf'] >= 0.95]).to_csv(coloc95_tsv_gz_path, sep='\t', index=False, header=True)
            # else:  # append
            coloc_df.to_csv(coloc_tsv_gz_path, sep='\t', index=False, header=False, mode='a')
            (coloc_df.loc[coloc_df['PP_H4_abf'] >= 0.5]).to_csv(coloc5_tsv_gz_path, sep='\t', index=False, header=False, mode='a')
            (coloc_df.loc[coloc_df['PP_H4_abf'] >= 0.8]).to_csv(coloc8_tsv_gz_path, sep='\t', index=False, header=False, mode='a')
            (coloc_df.loc[coloc_df['PP_H4_abf'] >= 0.95]).to_csv(coloc95_tsv_gz_path, sep='\t', index=False, header=False, mode='a')
