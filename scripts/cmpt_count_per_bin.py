from eqtl2gwas_pleiotropy.EBIeQTLinfo import EBIeQTLinfo
from eqtl2gwas_pleiotropy.OpenGWASinfo import OpenGWASinfo
from eqtl2gwas_pleiotropy.PathManager import PathManager
from eqtl2gwas_pleiotropy.constants import coloc_raw_tsv_path

import os
import pandas
import pathlib
import sys


#%% Parameters
if not '__file__' in locals():
    __file__ = "cmpt_count_per_bin.py"
if not os.path.isfile(coloc_raw_tsv_path):
    print("input file does not exit")
    sys.exit(1)

outdir_path = os.path.join(PathManager.get_project_path(), "out", os.path.basename(__file__))
pathlib.Path(outdir_path).mkdir(parents=True, exist_ok=True)

coloc_h4_tsv_path = os.path.join(PathManager.get_project_path(), "out", 'filter_h4.py/coloc_h4.tsv')

#%% Input1
h4_annot_tsv_path = os.path.join(PathManager.get_outdir_path(), "annotate_h4.py", "h4_annotated.tsv")
if not os.path.isfile(h4_annot_tsv_path):
    print("input file does not exit")
    sys.exit(1)

#%%
coloc_df = pandas.read_csv(h4_annot_tsv_path, sep="\t")

#%% bin
coloc_df['pos100k'] = round(coloc_df['pos']/100000).astype('int')

#%% rsid count in bins
bin_count_rsid_df = coloc_df[['chrom', 'pos100k', 'rsid']].drop_duplicates().groupby(['chrom', 'pos100k']).agg({'rsid': ['size', (lambda x: (",".join(sorted(set(x)))))]}).reset_index()
bin_count_rsid_df.columns = ['chrom', 'pos100k', 'rsid_count', 'rsid_lst']
bin_count_rsid_df = bin_count_rsid_df.sort_values(by=['rsid_count', 'rsid_lst', 'chrom', 'pos100k'], ascending=[False, True, True, True])
tsv_path = os.path.join(outdir_path, "bin_rsid_count.tsv")
bin_count_rsid_df.to_csv(tsv_path, sep="\t", index=False)

#%% etissue pleiotropy
pleio_etissue_df = coloc_df[['chrom', 'pos100k', 'etissue_subcategory']].drop_duplicates().groupby(['chrom', 'pos100k']).agg({'etissue_subcategory': ['size', (lambda x: (",".join(sorted(set(x)))))]}).reset_index()
pleio_etissue_df.columns = ['chrom', 'pos100k', 'etissue_subcategory_count', 'etissue_subcategory_lst']
pleio_etissue_df = pleio_etissue_df.sort_values(by=['etissue_subcategory_count', 'etissue_subcategory_lst', 'chrom', 'pos100k'], ascending=[False, True, True, True])
tsv_path = os.path.join(outdir_path, "bin_etissue_count.tsv")
pleio_etissue_df.to_csv(tsv_path, sep="\t", index=False)

#%% gwas pleiotropy
pleio_gwas_df = coloc_df[['chrom', 'pos100k', 'gwas_subcategory']].drop_duplicates().groupby(['chrom', 'pos100k']).agg({'gwas_subcategory': ['size', (lambda x: (",".join(sorted(set(x)))))]}).reset_index()
pleio_gwas_df.columns = ['chrom', 'pos100k', 'gwas_subcategory_count', 'gwas_subcategory_lst']
pleio_gwas_df = pleio_gwas_df.sort_values(by=['gwas_subcategory_count', 'gwas_subcategory_lst', 'chrom', 'pos100k'], ascending=[False, True, True, True])
tsv_path = os.path.join(outdir_path, "bin_gwas_count.tsv")
pleio_gwas_df.to_csv(tsv_path, sep="\t", index=False)

#%% egene pleiotropy
pleio_egene_df = coloc_df[['chrom', 'pos100k', 'egene', 'egene_symbol']].drop_duplicates().groupby(['chrom', 'pos100k']).agg({'egene': ['size', (lambda x: (",".join(sorted(set(x)))))], 'egene_symbol': (lambda x: (",".join(sorted([str(i) for i in x]))))}).reset_index()
pleio_egene_df.columns = ['chrom', 'pos100k', 'egene_count', 'egene_lst', 'egene_symbol_lst']
pleio_egene_df = pleio_egene_df.sort_values(by=['egene_count', 'egene_lst', 'egene_symbol_lst', 'chrom', 'pos100k'], ascending=[False, True, True, True, True])
tsv_path = os.path.join(outdir_path, "bin_egene_count.tsv")
pleio_egene_df.to_csv(tsv_path, sep="\t", index=False)
