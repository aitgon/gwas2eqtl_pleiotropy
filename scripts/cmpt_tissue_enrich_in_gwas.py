"""This script will compute a tissue enrichment in a given GWAS with the following algorithm.
Count all loci of each gwas and of each eqtl
Then count the proportion of loci of a given eqtl in a given gwas
Binomial test to evaluate significance
FDR 0.05"""
import os

import numpy
import pandas
import pathlib
import sys

from scipy.stats import binomtest
from statsmodels.stats.multitest import fdrcorrection
from eqtl2gwas_pleiotropy.PathManager import PathManager


#%%
help_cmd_str = "todo"
try:
    h4_tsv_path = sys.argv[1]
    tsv_path = sys.argv[2]
    if len(sys.argv) > 3:
        print("""Two many arguments!
        {}""".format(help_cmd_str))
        sys.exit(1)
except IndexError:
    print("""Argument missing!
    {}""".format(help_cmd_str))
    sys.exit(1)


# #%% Parameters
# if not '__file__' in locals():
#     __file__ = "cmpt_tissue_enrich_in_gwas.py"
outdir_path = os.path.dirname(tsv_path)
pathlib.Path(outdir_path).mkdir(parents=True, exist_ok=True)

#%%
# h4_tsv_path = "/home/gonzalez/Repositories/eqtl2gwas_pleiotropy/out/gwas413/genome/5e-08/1000000/filter_h4.py/h4.tsv"
df_raw = pandas.read_csv(h4_tsv_path, sep="\t")

#%%
columns = ['chrom', 'pos', 'rsid', 'eqtl_identifier', 'etissue_category', 'gwas_identifier', 'gwas_trait', 'gwas_category']
df = df_raw[columns].drop_duplicates()

#%% nb loci for each gwas
ngwas_df = df[['chrom', 'pos', 'rsid', 'gwas_identifier', 'gwas_trait']].drop_duplicates().groupby(['gwas_identifier', 'gwas_trait']).size().reset_index()
ngwas_df.rename({0: 'count_gwas'}, axis=1, inplace=1)
ngwas_df = ngwas_df.sort_values(by=['count_gwas'], ascending=False)
ngwas_df = ngwas_df.loc[ngwas_df['count_gwas'] >= 50, ]  # keep only if more than 50 loci

#%%
df = df.merge(ngwas_df[['gwas_identifier', 'gwas_trait']], on=['gwas_identifier', 'gwas_trait'])

#%% nb loci for each eqtl
neqtl_df = df[['chrom', 'pos', 'rsid', 'eqtl_identifier']].drop_duplicates().groupby(['eqtl_identifier']).size().reset_index()
neqtl_df.rename({0: 'count_eqtl'}, axis=1, inplace=1)
neqtl_df = neqtl_df.sort_values(by=['count_eqtl'], ascending=False)
neqtl_df = neqtl_df.loc[neqtl_df['count_eqtl'] >= 50, ]  # keep only if more than 50 loci

#%%
df = df.merge(neqtl_df['eqtl_identifier'], on=['eqtl_identifier'])

#%%
n = df[['chrom', 'pos', 'rsid']].drop_duplicates().shape[0]
neqtl_df['prop_eqtl'] = neqtl_df['count_eqtl'] / n

#%% nb loci with both gwas and tissue
ngwas_eqtl_df = df.groupby(['gwas_identifier', 'gwas_trait', 'eqtl_identifier']).size().reset_index()
ngwas_eqtl_df.rename({0: 'count_eqtl_in_gwas'}, axis=1, inplace=1)
ngwas_eqtl_df = ngwas_eqtl_df.sort_values(by=['count_eqtl_in_gwas'], ascending=False)
ngwas_eqtl_df = ngwas_eqtl_df.loc[ngwas_eqtl_df['count_eqtl_in_gwas'] >= 5, ]

#%%
df2 = ngwas_eqtl_df.merge(ngwas_df, on=['gwas_identifier', 'gwas_trait'])
df2 = df2.merge(neqtl_df[['eqtl_identifier', 'prop_eqtl']], on=['eqtl_identifier'])
df2['prop_eqtl_in_gwas'] = df2['count_eqtl_in_gwas'] / df2['count_gwas']

#%%
df2['log2fc_prop_gwas_in_eqtl'] = numpy.log2(df2['prop_eqtl_in_gwas']/df2['prop_eqtl'])

#%%
df2['p_binom'] = df2.apply(lambda x: binomtest(k=x['count_eqtl_in_gwas'], n=x['count_gwas'], p=x['prop_eqtl']).pvalue, axis=1)
df2['p_fdr5'] = fdrcorrection(df2['p_binom'])[1]
df2 = df2.loc[df2['p_fdr5'] < 0.05, ]
df2 = df2.sort_values(by=['p_fdr5', 'eqtl_identifier'])
# tsv_path = os.path.join(outdir_path, "tissue_enrich.tsv")
df2.to_csv(tsv_path, sep='\t', index=False)
