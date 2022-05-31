"""This script will compute a tissue enrichment in a given GWAS"""
import os
import pathlib
import pandas

from scipy.stats import binomtest
from statsmodels.stats.multitest import multipletests
from eqtl2gwas_pleiotropy.PathManager import PathManager
from eqtl2gwas_pleiotropy.constants import coloc_h4_tsv_path

#%% Parameters
if not '__file__' in locals():
    __file__ = "cmpt_tissue_enrich.py"
outdir_path = os.path.join(PathManager.get_project_path(), "out", os.path.basename(__file__))
pathlib.Path(outdir_path).mkdir(parents=True, exist_ok=True)

#%%
df_raw = pandas.read_csv(coloc_h4_tsv_path, sep="\t")

#%%
columns = ['chrom', 'pos', 'rsid', 'eqtl_identifier', 'gwas_identifier', 'gwas_trait_name']
df = df_raw[columns].drop_duplicates()

#%% nb loci for each gwas
ngwas_df = df[['chrom', 'pos', 'rsid', 'gwas_identifier', 'gwas_trait_name']].drop_duplicates().groupby(['gwas_identifier', 'gwas_trait_name']).size().reset_index()
ngwas_df.rename({0: 'count_gwas'}, axis=1, inplace=1)
ngwas_df = ngwas_df.sort_values(by=['count_gwas'], ascending=False)
ngwas_df = ngwas_df.loc[ngwas_df['count_gwas'] >= 50, ]

#%%
df = df.merge(ngwas_df[['gwas_identifier', 'gwas_trait_name']], on=['gwas_identifier', 'gwas_trait_name'])

#%% nb loci for each eqtl
neqtl_df = df[['chrom', 'pos', 'rsid', 'eqtl_identifier']].drop_duplicates().groupby(['eqtl_identifier']).size().reset_index()
neqtl_df.rename({0: 'count_eqtl'}, axis=1, inplace=1)
neqtl_df = neqtl_df.sort_values(by=['count_eqtl'], ascending=False)
neqtl_df = neqtl_df.loc[neqtl_df['count_eqtl'] >= 50,]

#%%
df = df.merge(neqtl_df['eqtl_identifier'], on=['eqtl_identifier'])

#%%
n = df[['chrom', 'pos', 'rsid']].drop_duplicates().shape[0]
neqtl_df['p_eqtl'] = neqtl_df['count_eqtl'] / n

#%% nb loci with both gwas1 and tissue
ngwas_eqtl_df = df.groupby(['gwas_identifier', 'gwas_trait_name', 'eqtl_identifier']).size().reset_index()
ngwas_eqtl_df.rename({0: 'count_gwas_eqtl'}, axis=1, inplace=1)
ngwas_eqtl_df = ngwas_eqtl_df.sort_values(by=['count_gwas_eqtl'], ascending=False)
ngwas_eqtl_df = ngwas_eqtl_df.loc[ngwas_eqtl_df['count_gwas_eqtl'] >= 5, ]

#%%
df2 = ngwas_eqtl_df.merge(ngwas_df, on=['gwas_identifier', 'gwas_trait_name'])
df2 = df2.merge(neqtl_df[['eqtl_identifier', 'p_eqtl']], on=['eqtl_identifier'])
df2['p_gwas_eqtl'] = df2['count_gwas_eqtl'] / df2['count_gwas']

#%%
df2['delta_p_gwas_eqtl'] = df2['p_gwas_eqtl']-df2['p_eqtl']

#%%
df2['pvalue'] = df2.apply(lambda x: binomtest(k=x['count_gwas_eqtl'], n=x['count_gwas'], p=x['p_eqtl']).pvalue, axis=1)
df2['p_bonferroni'] = multipletests(df2['pvalue'], method='bonferroni')[1]
df2 = df2.loc[df2['p_bonferroni'] < 0.05, ]
df2 = df2.sort_values(by=['eqtl_identifier', 'gwas_identifier', 'gwas_trait_name'])
tsv_path = os.path.join(outdir_path, "tissue_enrich.tsv")
df2.to_csv(tsv_path, sep='\t', index=False)