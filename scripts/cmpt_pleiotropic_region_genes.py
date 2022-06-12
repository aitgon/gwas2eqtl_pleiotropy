import os
import pandas
import pathlib
import shlex
import subprocess

from eqtl2gwas_pleiotropy.Logger import Logger
from eqtl2gwas_pleiotropy.PathManager import PathManager


#%% outdir path
if not '__file__' in locals():
    __file__ = "cmpt_pleiotropic_region_genes.py"
outdir_path = os.path.join(PathManager.get_project_path(), "out", os.path.basename(__file__))
pathlib.Path(outdir_path).mkdir(parents=True, exist_ok=True)

#%% input variant_to_regions
variant_to_region_tsv_path = os.path.join(PathManager.get_outdir_path(), "cmpt_variant_to_region.py", "variant_to_regions.tsv")
variant_to_region_df = pandas.read_csv(variant_to_region_tsv_path, sep="\t")

# import pdb; pdb.set_trace()

while variant_to_region_df.shape[0] > 0:
    region_pleio_count = variant_to_region_df['gwas_category_count'].max()
    print(region_pleio_count)
    egene_pleio_lst = sorted(variant_to_region_df.loc[variant_to_region_df['gwas_category_count'] == region_pleio_count, 'egene'].unique().tolist())
    # try:
    egene_symbol_pleio_lst = sorted(variant_to_region_df.loc[(variant_to_region_df['gwas_category_count'] == region_pleio_count) & (~variant_to_region_df['egene_symbol'].isna()), 'egene_symbol'].unique().tolist())
    # except:
    #     import pdb; pdb.set_trace()
    variant_to_region_df.drop(variant_to_region_df.loc[variant_to_region_df['gwas_category_count'] == region_pleio_count].index, inplace=True)

    pleio_egene_txt_path = os.path.join(outdir_path, "egene_pleio_{}.txt".format(region_pleio_count))
    with open(pleio_egene_txt_path, "w") as output:
        for item in egene_pleio_lst:
            output.write(item + "\n")
    pleio_egene_symbol_txt_path = os.path.join(outdir_path, "egene_symbol_pleio_{}.txt".format(region_pleio_count))
    with open(pleio_egene_symbol_txt_path, "w") as output:
        for item in egene_symbol_pleio_lst:
            output.write(item + "\n")
