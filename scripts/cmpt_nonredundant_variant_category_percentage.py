import shlex
import subprocess
import pandas
from eqtl2gwas_pleiotropy.Logger import Logger

#%%
pleiotropy = 1

#%%
region_window_bed_path = "/home/gonzalez/Repositories/eqtl2gwas_pleiotropy/out/cmpt_pleiotropic_regions.py/region_window_100000.bed"

#%%
variant_pleio_flank_0_bed_path = "/home/gonzalez/Repositories/eqtl2gwas_pleiotropy/out/cmpt_count_per_rsid.py/variant_pleio_{}_flank_0.bed".format(pleiotropy)

#%%
variant_region_bed = "variant_region.bed"

#%% intersect
cmd = "bedtools intersect -a {} -b {} -wb"
cmd = cmd.format(variant_pleio_flank_0_bed_path, region_window_bed_path)
Logger.info(cmd)
with open(variant_region_bed, 'w') as fout:
    result = subprocess.run(shlex.split(cmd), stdout=fout)

#%%
df_variant_region = pandas.read_csv(variant_region_bed, header=None, sep="\t")[range(9)]

#%%
df_variant_region.columns = ['chrom', 'start', 'end', 'rsid', 'gwas_subcategory_count', 'gwas_subcategory_lst', 'region_chrom', 'region_start', 'region_end']

#%%
df_variant_region = df_variant_region.drop_duplicates(['region_chrom', 'region_start', 'region_end'])
loci_variant_count = df_variant_region.shape[0]

#%%
df_long = df_variant_region.copy()
df_long["gwas_subcategory_lst"] = df_long["gwas_subcategory_lst"].str.split(",")
df_long = df_long.explode("gwas_subcategory_lst")

for category in df_long['gwas_subcategory_lst'].unique():
    cat_variant_count = df_long.loc[(df_long['gwas_subcategory_lst'] == category)].shape[0]
    print(category, loci_variant_count, cat_variant_count, round(cat_variant_count/loci_variant_count*100))
