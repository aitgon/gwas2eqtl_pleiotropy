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
tissue_enrich_tsv_path = "out/cmpt_tissue_enrich_in_gwas.py/tissue_enrich.tsv"
# gwas_info_df = pandas.read_csv("out/download/gwas-api.mrcieu.ac.uk/gwasinfo/opengwas.tsv", sep="\t")
gwas_info_df = OpenGWASinfo().df

#%%
df_raw = pandas.read_csv(tissue_enrich_tsv_path, sep="\t")
df = df_raw[['gwas_identifier', 'gwas_trait_name', 'eqtl_identifier', 'delta_p_gwas_eqtl']]

#%%
# df['gwas'] = df['gwas_identifier'].str.cat(df['gwas_trait_name'], sep="|")
# df.drop(['gwas_trait_name'], axis=1, inplace=True)

#%% gwas subcategories
# gwas_identifier_subcategory_df = gwas_info_df[['gwas_identifier', 'subcategory']]
# gwas_identifier_subcategory_df = gwas_identifier_subcategory_df.loc[gwas_identifier_subcategory_df['gwas_identifier'].isin(df.gwas_identifier.tolist())]
df = df.merge(gwas_info_df[['gwas_identifier', 'subcategory', 'pmid']], on="gwas_identifier")

#%%
for i, subcategory in enumerate(df['subcategory'].unique()):
    if isinstance(subcategory, str):
        subcategory_fn = subcategory.replace(" ", "_").replace("/", "")
        df_gwas = df.loc[df['subcategory'] == subcategory, ["eqtl_identifier", "gwas_identifier", "delta_p_gwas_eqtl"]]
    else:
        subcategory_fn = "NA"
        df_gwas = df.loc[df['subcategory'].isna(), ["eqtl_identifier", "gwas_identifier", "delta_p_gwas_eqtl"]]
    wide_df = df_gwas.pivot_table(index="eqtl_identifier", columns="gwas_identifier", values="delta_p_gwas_eqtl",
                             fill_value=0)
    wide_df.columns = df[['gwas_identifier', 'gwas_trait_name']].drop_duplicates().merge(
        pandas.DataFrame({'gwas_identifier': wide_df.columns.tolist()}), on='gwas_identifier')[
        'gwas_trait_name'].tolist()
    g = seaborn.clustermap(wide_df, center=0, cmap="vlag", mask=(wide_df == 0))
    g.fig.suptitle(subcategory)
    heatmap_gwas_eqtl_png_path = os.path.join(outdir_path, subcategory_fn)
    plt.savefig(heatmap_gwas_eqtl_png_path)

