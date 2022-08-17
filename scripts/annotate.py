from eqtl2gwas_pleiotropy.EBIeQTLinfo import EBIeQTLinfo
from eqtl2gwas_pleiotropy.Logger import Logger

import os
import pandas
import pathlib
import sys


#%%
from eqtl2gwas_pleiotropy.URL import URL

help_cmd_str = "todo"
try:
    pp_h4_abf = float(sys.argv[1])
    coloc_tsv_gz_path = sys.argv[2]
    gwas_cat_ods_path = sys.argv[3]
    etissue_cat_ods_path = sys.argv[4]
    annotated_tsv_gz_path = sys.argv[5]
    if len(sys.argv) > 6:
        print("""Two many arguments!
        {}""".format(help_cmd_str))
        sys.exit(1)
except IndexError:
    print("""Argument missing!
    {}""".format(help_cmd_str))
    sys.exit(1)

pathlib.Path(os.path.dirname(annotated_tsv_gz_path)).mkdir(parents=True, exist_ok=True)

#%% Download eQTL annotations
eqtl_info_df = EBIeQTLinfo().df
etissue_cat_df = pandas.read_excel(etissue_cat_ods_path, engine="odf")
eqtl_info_df = eqtl_info_df.merge(etissue_cat_df, on=['study', 'qtl_group', 'tissue_ontology_id', 'tissue_ontology_term', 'tissue_label', 'condition_label'])
eqtl_info_df.rename({'identifier': "eqtl_identifier"}, axis=1, inplace=True)

#%% gwas category
gwas_annot_df = pandas.read_excel(gwas_cat_ods_path, engine="odf")
gwas_annot_df.rename({'id': "gwas_identifier", 'subcategory': 'gwas_category_mrcieu', 'manual_category': 'gwas_category', 'trait': 'gwas_trait'}, axis=1, inplace=True)

#%%
Logger.info("Reading {}".format(coloc_tsv_gz_path))
coloc_df = pandas.read_csv(coloc_tsv_gz_path, sep="\t")
coloc_df = coloc_df.loc[coloc_df['PP.H4.abf'] >= pp_h4_abf]
coloc_df.rename({'pp_h4': "pp_h4", }, axis=1, inplace=True)

#%%
coloc_df = coloc_df.merge(gwas_annot_df[["gwas_identifier", "gwas_trait", 'gwas_category_mrcieu', 'gwas_category']], on='gwas_identifier')
coloc_df = coloc_df.merge(eqtl_info_df[['eqtl_identifier', 'etissue_category']].drop_duplicates(), on='eqtl_identifier')

#%% cytoband
Logger.info("Annotate cytobands")
cytoband_url = "http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/cytoBand.txt.gz"
cytoband_path = URL(cytoband_url).download()
cyto_df = pandas.read_csv(cytoband_path, sep="\t", header=None, usecols=[0, 1, 2, 3])
cyto_df.columns = ['chrom', 'start', 'end', 'cytoband']
cyto_df['chrom'] = cyto_df['chrom'].str.replace('chr', '')
coloc_df['chrom'] = coloc_df['chrom'].copy().astype('str')
coloc_df['cytoband'] = None
for rowi, row in cyto_df.iterrows():
    chrom = row['chrom']
    start = row['start']
    end = row['end']
    cytoband = row['cytoband']
    coloc_df.loc[
        (coloc_df['chrom'] == chrom) & (start <= coloc_df['pos']) & (
                coloc_df['pos'] <= end), 'cytoband'] = cytoband

#%%
columns = ['chrom', 'pos', 'cytoband', 'rsid', 'ref', 'alt', 'egene', 'egene_symbol',
           'eqtl_beta', 'eqtl_pvalue', 'eqtl_identifier', 'etissue_category', 'gwas_beta',
           'gwas_pvalue', 'gwas_identifier', 'gwas_trait', 'gwas_category',
           'SNP.PP.H4', 'PP.H4.abf', 'coloc_window', 'nsnps', 'PP.H3.abf', 'PP.H2.abf', 'PP.H1.abf', 'PP.H0.abf']
coloc_df = coloc_df[columns]

#%% TSV
coloc_df.sort_values(by=coloc_df.columns.tolist(), inplace=True)
Logger.info("Writing {}".format(annotated_tsv_gz_path))
coloc_df.to_csv(annotated_tsv_gz_path, sep="\t", index=True, index_label='id', na_rep='na')
