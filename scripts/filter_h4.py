import pdb

from sqlalchemy import create_engine, select

from gwas2eqtl_pleiotropy.EBIeQTLinfo import EBIeQTLinfo
from gwas2eqtl_pleiotropy.Logger import Logger
from gwas2eqtl_pleiotropy.constants import h4_cutoff, snp_h4_cutoff

import os
import pandas
import pathlib
import sys

from gwas2eqtl_pleiotropy.db import meta, coloc, gwas_annot_tbl, ensg2symbol_tbl, rsid2cytoband_tbl

#%%
help_cmd_str = "todo"
try:
    url_db = sys.argv[1]
    etissue_cat_ods_path = sys.argv[2]
    h4_tsv_gz_path = sys.argv[3]
    if len(sys.argv) > 4:
        print("""Two many arguments!
        {}""".format(help_cmd_str))
        sys.exit(1)
except IndexError:
    print("""Argument missing!
    {}""".format(help_cmd_str))
    sys.exit(1)

outdir_path = os.path.dirname(h4_tsv_gz_path)
pathlib.Path(outdir_path).mkdir(parents=True, exist_ok=True)

#%% Insert into coloc
engine = create_engine(url_db, echo = False)
meta.create_all(engine)

#%%
Logger.info("Select colocalization joins")

sel_cols = [coloc, gwas_annot_tbl.c.gwas_class, ensg2symbol_tbl.c.symbol, rsid2cytoband_tbl.c.cytoband]
sel = select(sel_cols).distinct()
sel = sel.where(coloc.c.gwas_id == gwas_annot_tbl.c.gwas_id)
sel = sel.where(coloc.c.egene == ensg2symbol_tbl.c.gene_id)
sel = sel.where(coloc.c.rsid == rsid2cytoband_tbl.c.rsid)
sel = sel.where(coloc.c['PP.H4.abf'] >= h4_cutoff)
with engine.connect() as con:
    result_lst = con.execute(sel).fetchall()
Logger.info("End select colocalization joins")

h4_df = pandas.DataFrame(result_lst)
h4_df.drop(['id'], inplace=True, axis=1)
h4_df.rename({'symbol': 'egene_symbol'}, axis=1, inplace=1)

#%% Download eQTL annotations
eqtl_info_df = EBIeQTLinfo().df
etissue_cat_df = pandas.read_excel(etissue_cat_ods_path, engine="odf")
eqtl_info_df = eqtl_info_df.merge(etissue_cat_df, on=['study', 'qtl_group', 'tissue_ontology_id', 'tissue_ontology_term', 'tissue_label', 'condition_label'])
eqtl_info_df.rename({'identifier': "eqtl_id"}, axis=1, inplace=True)
h4_df = h4_df.merge(eqtl_info_df[['eqtl_id', 'etissue_category']].drop_duplicates(), on='eqtl_id')

#%% how many coloc loci?
loci_h4_df = h4_df.sort_values(by=['PP.H4.abf', 'SNP.PP.H4'], ascending=False)
loci_h4_df = loci_h4_df.drop_duplicates(subset=['PP.H4.abf', 'coloc_region', 'nsnps'], keep='first')
loci_coloc_h4_path = os.path.join(outdir_path, "loci_h4.tsv")
loci_h4_df.to_csv(loci_coloc_h4_path, sep="\t", index=False)

#%% how many causal variants, SNP.PP.H4>=0.5?
variants_h4_df = h4_df.sort_values(by=['SNP.PP.H4'], ascending=False)
variants_h4_df = variants_h4_df.drop_duplicates(subset=['chrom', 'pos', 'cytoband', 'rsid'], keep='first')
variants_h4_df = variants_h4_df.loc[variants_h4_df['SNP.PP.H4'] >= snp_h4_cutoff, ]
variants_h4_path = os.path.join(outdir_path, "variants_h4.tsv")
variants_h4_df.to_csv(variants_h4_path, sep="\t", index=False)

#%%
h4_df.sort_values(['chrom', 'pos'], ascending=True, inplace=True)
h4_df.to_csv(h4_tsv_gz_path, sep="\t", index=False)
