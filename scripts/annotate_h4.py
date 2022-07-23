from eqtl2gwas_pleiotropy.EBIeQTLinfo import EBIeQTLinfo
from eqtl2gwas_pleiotropy.Logger import Logger
from eqtl2gwas_pleiotropy.OpenGWASinfo import OpenGWASinfo
from eqtl2gwas_pleiotropy.PathManager import PathManager
from eqtl2gwas_pleiotropy.constants import coloc_h4_tsv_path

import os
import pandas
import pathlib
import sys


#%% Parameters
# if not '__file__' in locals():
#     __file__ = "annotate_h4.py"
# if not os.path.isfile(coloc_h4_tsv_path):
#     print("input file does not exit")
#     sys.exit(1)
#
# outdir_path = os.path.join(PathManager.get_outdir_path(), os.path.basename(__file__))
# pathlib.Path(outdir_path).mkdir(parents=True, exist_ok=True)

#%%
help_cmd_str = "todo"
try:
    coloc_h4_tsv_path = sys.argv[1]
    manual_annot_ods_path = sys.argv[2]
    h4_annotated_tsv_path = sys.argv[3]
    h4_annotated_ods_path = sys.argv[4]
    h4_annotated_bed_path = sys.argv[5]
    if len(sys.argv) > 6:
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
eqtl_info_df.rename({'identifier': "eqtl_identifier", 'tissue_label': 'etissue_label'}, axis=1, inplace=True)

# #%% Download OpenGWAS annotations
# open_gwas_df = OpenGWASinfo().df
# open_gwas_df.rename({'subcategory': 'gwas_subcategory'}, axis=1, inplace=True)
# open_gwas_df = open_gwas_df[['gwas_identifier', 'gwas_subcategory']].drop_duplicates()
# # remove subcat with nan
# open_gwas_df = open_gwas_df.loc[~open_gwas_df["gwas_subcategory"].isna(), ]
gwas_annot_df = pandas.read_excel(manual_annot_ods_path, engine="odf")
gwas_annot_df.rename({'id': "gwas_identifier", 'subcategory': 'gwas_category_mrcieu', 'manual_category': 'gwas_category_eqtl2gwas', 'trait': 'gwas_trait'}, axis=1, inplace=True)

#%%
coloc_h4_df = pandas.read_csv(coloc_h4_tsv_path, sep="\t")

# import pdb; pdb.set_trace()
#%%
coloc_h4_df = coloc_h4_df.merge(gwas_annot_df[["gwas_identifier", "gwas_trait", 'gwas_category_mrcieu', 'gwas_category_eqtl2gwas']], on='gwas_identifier')
coloc_h4_df = coloc_h4_df.merge(eqtl_info_df[['eqtl_identifier', 'etissue_label']].drop_duplicates(), on='eqtl_identifier')

#%%
columns = ['chrom', 'pos', 'rsid', 'ref', 'alt', 'egene', 'egene_symbol',
           'eqtl_beta', 'eqtl_pvalue', 'eqtl_identifier', 'gwas_beta',
           'gwas_pvalue', 'gwas_identifier', 'gwas_trait', 'gwas_category_eqtl2gwas',
           'pp_h4', 'PP.H4.abf', 'coloc_window', 'nsnps', 'PP.H3.abf', 'PP.H2.abf', 'PP.H1.abf', 'PP.H0.abf', 'gwas_category_mrcieu', 'etissue_label']
coloc_h4_df = coloc_h4_df[columns]

#%% TSV
coloc_h4_df.sort_values(by=coloc_h4_df.columns.tolist(), inplace=True)
coloc_h4_df.to_csv(h4_annotated_tsv_path, sep="\t", index=False)
Logger.info("Writing {}".format(h4_annotated_tsv_path))
with pandas.ExcelWriter(h4_annotated_ods_path) as fout:
    coloc_h4_df.to_excel(fout, index=False)

#%% BED
coloc_h4_bed_df = coloc_h4_df.rename({'chrom': '#chrom', 'pos': 'end'}, axis=1)
coloc_h4_bed_df['#chrom'] = 'chr' + coloc_h4_bed_df['#chrom'].astype('str')
# coloc_h4_bed_df['start'] = coloc_h4_bed_df['end'] - 1
coloc_h4_bed_df.insert(1, 'start', coloc_h4_bed_df['end'] - 1)
# h4_annotated_bed_path = os.path.join(outdir_path, 'h4_annotated.bed')
coloc_h4_bed_df.sort_values(by=coloc_h4_bed_df.columns.tolist(), inplace=True)
coloc_h4_bed_df.to_csv(h4_annotated_bed_path, sep="\t", index=False)
