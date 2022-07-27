from eqtl2gwas_pleiotropy.EBIeQTLinfo import EBIeQTLinfo
from eqtl2gwas_pleiotropy.Logger import Logger
from eqtl2gwas_pleiotropy.OpenGWASinfo import OpenGWASinfo
from eqtl2gwas_pleiotropy.PathManager import PathManager
from eqtl2gwas_pleiotropy.constants import coloc_h4_tsv_path

import os
import pandas
import pathlib
import sys

#%%
help_cmd_str = "todo"
try:
    coloc_h4_tsv_path = sys.argv[1]
    gwas_cat_ods_path = sys.argv[2]
    etissue_cat_ods_path = sys.argv[3]
    h4_annotated_tsv_path = sys.argv[4]
    h4_annotated_ods_path = sys.argv[5]
    h4_annotated_bed_path = sys.argv[6]
    if len(sys.argv) > 7:
        print("""Two many arguments!
        {}""".format(help_cmd_str))
        sys.exit(1)
except IndexError:
    print("""Argument missing!
    {}""".format(help_cmd_str))
    sys.exit(1)
pathlib.Path(os.path.dirname(h4_annotated_ods_path)).mkdir(parents=True, exist_ok=True)
pathlib.Path(os.path.dirname(h4_annotated_bed_path)).mkdir(parents=True, exist_ok=True)

# h4_annotated_tsv_path = os.path.join(outdir_path, 'h4_annotated.tsv')
#%% Download eQTL annotations
eqtl_info_df = EBIeQTLinfo().df
etissue_cat_df = pandas.read_excel(etissue_cat_ods_path, engine="odf")
eqtl_info_df = eqtl_info_df.merge(etissue_cat_df, on=['study', 'qtl_group', 'tissue_ontology_id', 'tissue_ontology_term', 'tissue_label', 'condition_label'])
eqtl_info_df.rename({'identifier': "eqtl_identifier"}, axis=1, inplace=True)

#%% gwas category
gwas_annot_df = pandas.read_excel(gwas_cat_ods_path, engine="odf")
gwas_annot_df.rename({'id': "gwas_identifier", 'subcategory': 'gwas_category_mrcieu', 'manual_category': 'gwas_category_eqtl2gwas', 'trait': 'gwas_trait'}, axis=1, inplace=True)

#%%
coloc_h4_df = pandas.read_csv(coloc_h4_tsv_path, sep="\t")

# import pdb; pdb.set_trace()
#%%
coloc_h4_df = coloc_h4_df.merge(gwas_annot_df[["gwas_identifier", "gwas_trait", 'gwas_category_mrcieu', 'gwas_category_eqtl2gwas']], on='gwas_identifier')
coloc_h4_df = coloc_h4_df.merge(eqtl_info_df[['eqtl_identifier', 'etissue_category']].drop_duplicates(), on='eqtl_identifier')

#%%
columns = ['chrom', 'pos', 'rsid', 'ref', 'alt', 'egene', 'egene_symbol',
           'eqtl_beta', 'eqtl_pvalue', 'eqtl_identifier', 'etissue_category', 'gwas_beta',
           'gwas_pvalue', 'gwas_identifier', 'gwas_trait', 'gwas_category_eqtl2gwas',
           'pp_h4', 'PP.H4.abf', 'coloc_window', 'nsnps', 'PP.H3.abf', 'PP.H2.abf', 'PP.H1.abf', 'PP.H0.abf', 'gwas_category_mrcieu']
coloc_h4_df = coloc_h4_df[columns]

#%% TSV
coloc_h4_df.sort_values(by=coloc_h4_df.columns.tolist(), inplace=True)
Logger.info("Writing {}".format(h4_annotated_tsv_path))
coloc_h4_df.to_csv(h4_annotated_tsv_path, sep="\t", index=False, na_rep='na')
Logger.info("Writing {}".format(h4_annotated_ods_path))
with pandas.ExcelWriter(h4_annotated_ods_path) as fout:
    coloc_h4_df.to_excel(fout, index=False)

#%% BED
coloc_h4_bed_df = coloc_h4_df.rename({'chrom': '#chrom', 'pos': 'end'}, axis=1)
coloc_h4_bed_df['#chrom'] = 'chr' + coloc_h4_bed_df['#chrom'].astype('str')
coloc_h4_bed_df.insert(1, 'start', coloc_h4_bed_df['end'] - 1)
coloc_h4_bed_df.sort_values(by=coloc_h4_bed_df.columns.tolist(), inplace=True)
Logger.info("Writing {}".format(h4_annotated_bed_path))
coloc_h4_bed_df.to_csv(h4_annotated_bed_path, sep="\t", index=False, na_rep='na')
