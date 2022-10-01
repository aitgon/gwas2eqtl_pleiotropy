import os
import pandas
import pathlib
import shlex
import subprocess

from gwas2eqtl_pleiotropy.Logger import Logger
from gwas2eqtl_pleiotropy.PathManager import PathManager


#%% outdir path
if not '__file__' in locals():
    __file__ = "cmpt_variant_to_region.py"
outdir_path = os.path.join(PathManager.get_project_path(), "out", os.path.basename(__file__))
pathlib.Path(outdir_path).mkdir(parents=True, exist_ok=True)

#%% input h4_annotated
h4_annotated_bed_path = os.path.join(PathManager.get_outdir_path(), "annotate_db.py", "h4_annotated.bed")

#%% input region window
region_window_bed_path = os.path.join(PathManager.get_outdir_path(), "cmpt_pleiotropic_regions.py", "region_window_100000.bed")

#%%
variant_to_region_bed_path = os.path.join(outdir_path, "variant_to_regions.bed")

# %% intersect
cmd_fmt = "bedtools intersect -a {} -b {} -wb"
cmd_str = cmd_fmt.format(h4_annotated_bed_path, region_window_bed_path)
Logger.info(cmd_str)
with open(variant_to_region_bed_path, 'w') as fout:
    result = subprocess.run(shlex.split(cmd_str), stdout=fout)

#%% TSV
variant_to_region_bed_df = pandas.read_csv(variant_to_region_bed_path, sep="\t", header=None)
variant_to_region_tsv_path = os.path.join(outdir_path, "variant_to_regions.tsv")
variant_to_region_df = variant_to_region_bed_df.drop([18, 19, 20, 21, 22, 23, 24, 26], axis=1)
variant_to_region_df.columns = ['chrom', 'variant_start', 'variant_end', 'rsid', 'ref', 'alt', 'egene_symbol', 'egene', 'eqtl_beta', 'eqtl_pvalue', 'eqtl_identifier', 'gwas_beta', 'gwas_pvalue', 'gwas_identifier', 'gwas_trait_name', 'pp_h4', 'PP.H4.abf', 'gwas_subcategory', 'etissue_subcategory', 'region_start', 'region_end', 'gwas_category_count', 'gwas_category_lst']
variant_to_region_df.drop_duplicates(inplace=True)
variant_to_region_df.to_csv(variant_to_region_tsv_path, sep="\t", header=True, index=False)
