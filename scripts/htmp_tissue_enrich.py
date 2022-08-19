"""This script will compute a tissue enrichment in a given GWAS"""
from eqtl2gwas_pleiotropy.OpenGWASinfo import OpenGWASinfo
from eqtl2gwas_pleiotropy.PathManager import PathManager
from matplotlib import pyplot as plt

import os
import pandas
import pathlib
import seaborn


#%% Outdir
if not '__file__' in locals():
    __file__ = "htmp_tissue_enrich.py"
outdir_path = os.path.join(PathManager.get_project_path(), "out", os.path.basename(__file__))
pathlib.Path(outdir_path).mkdir(parents=True, exist_ok=True)

#%%
seaborn.set_theme(color_codes=True)
tissue_enrich_tsv_path = "/home/gonzalez/Repositories/eqtl2gwas_pleiotropy/out/gwas413/genome/5e-08/1000000/cmpt_tissue_enrich_in_gwas.py/tissue_enrich.tsv"
# gwas_info_df = pandas.read_csv("out/download/gwas-api.mrcieu.ac.uk/gwasinfo/opengwas.tsv", sep="\t")
gwas_category_path = "/home/gonzalez/Repositories/eqtl2gwas_pleiotropy/out/gwas413/genome/5e-08/1000000/filter_h4.py/h4.tsv"
gwas_category_df = pandas.read_csv(gwas_category_path, sep="\t", header=0, usecols=["gwas_identifier", "gwas_category"]).drop_duplicates()

#%%
df_raw = pandas.read_csv(tissue_enrich_tsv_path, sep="\t")
df = df_raw[['gwas_identifier', 'gwas_trait', 'eqtl_identifier', 'log2fc_prop_gwas_in_eqtl']]

#%%
# df['gwas'] = df['gwas_identifier'].str.cat(df['gwas_trait'], sep="|")
# df.drop(['gwas_trait'], axis=1, inplace=True)

#%% gwas category
# gwas_identifier_subcategory_df = gwas_info_df[['gwas_identifier', 'gwas_category']]
# gwas_identifier_subcategory_df = gwas_identifier_subcategory_df.loc[gwas_identifier_subcategory_df['gwas_identifier'].isin(df.gwas_identifier.tolist())]
# import pdb; pdb.set_trace()
df = df.merge(gwas_category_df[['gwas_identifier', 'gwas_category']], on="gwas_identifier")

#%%
for i, gwas_category in enumerate(df['gwas_category'].unique()):
    if isinstance(gwas_category, str):
        subcategory_fn = gwas_category.replace(" ", "_").replace("/", "")
        df_gwas = df.loc[df['gwas_category'] == gwas_category, ["eqtl_identifier", "gwas_identifier", "log2fc_prop_gwas_in_eqtl"]]
    else:
        subcategory_fn = "NA"
        df_gwas = df.loc[df['gwas_category'].isna(), ["eqtl_identifier", "gwas_identifier", "log2fc_prop_gwas_in_eqtl"]]
    wide_df = df_gwas.pivot_table(index="eqtl_identifier", columns="gwas_identifier", values="log2fc_prop_gwas_in_eqtl",
                             fill_value=0)
    try:
        wide_df.columns = df[['gwas_identifier', 'gwas_trait']].drop_duplicates().merge(
            pandas.DataFrame({'gwas_identifier': wide_df.columns.tolist()}), on='gwas_identifier')[
            'gwas_trait'].tolist()
        g = seaborn.clustermap(wide_df, center=0, cmap="vlag", mask=(wide_df == 0))
        g.fig.suptitle(gwas_category)
        heatmap_gwas_eqtl_png_path = os.path.join(outdir_path, subcategory_fn)
        plt.savefig(heatmap_gwas_eqtl_png_path)
    except:
        continue

