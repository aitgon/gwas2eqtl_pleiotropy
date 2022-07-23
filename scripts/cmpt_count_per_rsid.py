import shlex
import subprocess

from eqtl2gwas_pleiotropy.EBIeQTLinfo import EBIeQTLinfo
from eqtl2gwas_pleiotropy.OpenGWASinfo import OpenGWASinfo
from eqtl2gwas_pleiotropy.PathManager import PathManager
from eqtl2gwas_pleiotropy.URL import URL
from eqtl2gwas_pleiotropy.constants import coloc_raw_tsv_path

import os
import pandas
import pathlib
import sys


# #%% Parameters
# if not '__file__' in locals():
#     __file__ = "cmpt_count_per_rsid.py"
# if not os.path.isfile(coloc_raw_tsv_path):
#     print("input file does not exit")
#     sys.exit(1)
#
# outdir_path = os.path.join(PathManager.get_project_path(), "out", os.path.basename(__file__))
# pathlib.Path(outdir_path).mkdir(parents=True, exist_ok=True)

# coloc_h4_tsv_path = os.path.join(PathManager.get_outdir_path(), 'filter_h4.py", "coloc_h4.tsv')

#%%
help_cmd_str = "todo"
try:
    # coloc_h4_tsv_path = sys.argv[1]
    annotated_tsv_path = sys.argv[1]
    outdir_path = sys.argv[2]
    # h4_annotated_ods_path = sys.argv[3]
    # h4_annotated_bed_path = sys.argv[4]
    if len(sys.argv) > 3:
        print("""Two many arguments!
        {}""".format(help_cmd_str))
        sys.exit(1)
except IndexError:
    print("""Argument missing!
    {}""".format(help_cmd_str))
    sys.exit(1)
pathlib.Path(outdir_path).mkdir(parents=True, exist_ok=True)

# #%% Input1
# h4_annot_tsv_path = os.path.join(PathManager.get_outdir_path(), "annotate_h4.py", "h4_annotated.tsv")
# if not os.path.isfile(h4_annot_tsv_path):
#     print("input file does not exit")
#     sys.exit(1)

#%%
coloc_df = pandas.read_csv(annotated_tsv_path, sep="\t")

#%% etissue pleiotropy
pleio_etissue_df = coloc_df[['chrom', 'pos', 'rsid', 'etissue_label']].drop_duplicates().groupby(['chrom', 'pos', 'rsid']).agg({'etissue_label': ['size', (lambda x: (",".join(sorted(x))))]}).reset_index()
pleio_etissue_df.columns = ['chrom', 'pos', 'rsid', 'etissue_label_count', 'etissue_subcategory_lst']
pleio_etissue_df = pleio_etissue_df.sort_values(by=['etissue_label_count', 'etissue_subcategory_lst', 'rsid'], ascending=[False, True, True])
tsv_path = os.path.join(outdir_path, "count_per_rsid_etissue.tsv")
pleio_etissue_df.to_csv(tsv_path, sep="\t", index=False)

#%% egene pleiotropy
pleio_egene_df = coloc_df[['chrom', 'pos', 'rsid', 'egene', 'egene_symbol']].drop_duplicates().groupby(['chrom', 'pos', 'rsid']).agg({'egene': ['size', (lambda x: (",".join(sorted(x))))], 'egene_symbol': (lambda x: (",".join(sorted([str(i) for i in x]))))}).reset_index()
pleio_egene_df.columns = ['chrom', 'pos', 'rsid', 'egene_count', 'egene_lst', 'egene_symbol_lst']
pleio_egene_df = pleio_egene_df.sort_values(by=['egene_count', 'egene_lst', 'egene_symbol_lst', 'rsid'], ascending=[False, True, True, True])
tsv_path = os.path.join(outdir_path, "count_per_rsid_egene.tsv")
pleio_egene_df.to_csv(tsv_path, sep="\t", index=False)

#%% gwas pleiotropy
pleio_gwas_df = coloc_df[['chrom', 'pos', 'rsid', 'gwas_category_eqtl2gwas']].drop_duplicates().groupby(['chrom', 'pos', 'rsid']).agg({'gwas_category_eqtl2gwas': ['size', (lambda x: (",".join(sorted(x))))]}).reset_index()
# import pdb; pdb.set_trace()
# pleio_gwas_df.columns = ['chrom', 'pos', 'rsid', 'gwas_subcategory_count', 'gwas_subcategory_lst']
pleio_gwas_df.columns = ['chrom', 'pos', 'rsid', 'gwas_category_count', 'gwas_category_lst']
pleio_gwas_df = pleio_gwas_df.sort_values(by=['gwas_category_count', 'gwas_category_lst', 'rsid'], ascending=[False, True, True])
tsv_path = os.path.join(outdir_path, "count_per_rsid_gwas.tsv")
pleio_gwas_df.to_csv(tsv_path, sep="\t", index=False)

#%########################################### bed files, flanking=0
# bed files of variants splitted by gwas categories
flank = 0
variant_bed_df = pleio_gwas_df.copy()
variant_bed_df['chrom'] = 'chr' + variant_bed_df['chrom'].astype(str)
variant_bed_df['start'] = variant_bed_df['pos'] - 1 - flank
variant_bed_df['end'] = variant_bed_df['pos'] + flank
variant_bed_df = variant_bed_df[['chrom', 'start', 'end', 'rsid', 'gwas_category_count', 'gwas_category_lst']]

for count_pleio in range(1, 6):
    variant_pleio_i_bed_path = os.path.join(outdir_path, "variant_pleio_{}_flank_{}.bed".format(count_pleio, flank))
    variant_pleio_i_bed_df = variant_bed_df.loc[variant_bed_df['gwas_category_count'] == count_pleio, ]
    variant_pleio_i_bed_df = variant_pleio_i_bed_df.sort_values(by=['chrom', 'start', 'end'])
    variant_pleio_i_bed_df.to_csv(variant_pleio_i_bed_path, sep="\t", index=False, header=False)

    # liftover to hg19
    hg38ToHg19_chain_path = 'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz'

    # %% Chain
    variant_hg19_pleio_i_bed_path = os.path.join(outdir_path, "variant_hg19_pleio_{}_flank_{}.bed".format(count_pleio, flank))
    urlobj = URL(url=hg38ToHg19_chain_path)
    hg19tohg38_path = urlobj.download()

    crossmap_cmd = "CrossMap.py bed {chain} {inbed} {outbed}".format(
        chain=hg19tohg38_path, inbed=variant_pleio_i_bed_path,
        outbed=variant_hg19_pleio_i_bed_path)
    # %%
    crossmap_cmd_lst = shlex.split(crossmap_cmd)
    subprocess.run(crossmap_cmd_lst)

#%########################################### bed files, flanking=100
flank = 50
# bed files of variants splitted by gwas categories
variant_bed_df = pleio_gwas_df.copy()
variant_bed_df['chrom'] = 'chr' + variant_bed_df['chrom'].astype(str)
variant_bed_df['start'] = variant_bed_df['pos'] - 1 - flank
variant_bed_df['end'] = variant_bed_df['pos'] + flank
variant_bed_df = variant_bed_df[['chrom', 'start', 'end', 'rsid', 'gwas_category_count', 'gwas_category_lst']]

for count_pleio in range(1, 6):
    variant_pleio_i_bed_path = os.path.join(outdir_path, "variant_pleio_{}_flank_{}.bed".format(count_pleio, flank))
    variant_pleio_i_bed_df = variant_bed_df.loc[variant_bed_df['gwas_category_count'] == count_pleio, ]
    variant_pleio_i_bed_df = variant_pleio_i_bed_df.sort_values(by=['chrom', 'start', 'end'])
    variant_pleio_i_bed_df.to_csv(variant_pleio_i_bed_path, sep="\t", index=False, header=False)
