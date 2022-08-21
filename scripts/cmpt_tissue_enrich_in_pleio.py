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

from scipy.stats import fisher_exact
from statsmodels.stats.multitest import fdrcorrection


#%%
help_cmd_str = "todo"
try:
    h4_tsv_path = sys.argv[1]
    count_per_rsid_gwas_tsv_path = sys.argv[2]
    upper_var_gwas_cat_count = int(sys.argv[3])
    tsv_path = sys.argv[4]
    if len(sys.argv) > 5:
        print("""Two many arguments!
        {}""".format(help_cmd_str))
        sys.exit(1)
except IndexError:
    print("""Argument missing!
    {}""".format(help_cmd_str))
    sys.exit(1)


# #%% Parameters
outdir_path = os.path.dirname(tsv_path)
pathlib.Path(outdir_path).mkdir(parents=True, exist_ok=True)

#%%
h4_df = pandas.read_csv(h4_tsv_path, sep="\t")

#%%
count_per_rsid_gwas_df = pandas.read_csv(count_per_rsid_gwas_tsv_path, sep="\t")

#%%
m_df = h4_df.merge(count_per_rsid_gwas_df, on=['chrom', 'pos', 'cytoband', 'rsid'])

#%%
m2_df = m_df[['rsid', 'gwas_category_count', 'etissue_category']].drop_duplicates().sort_values('gwas_category_count', ascending=False)
m2_df.loc[m2_df['gwas_category_count'] >= upper_var_gwas_cat_count, 'gwas_category_count'] = upper_var_gwas_cat_count

#%%
g_target_df = m2_df.groupby(['gwas_category_count', 'etissue_category']).count().reset_index()

#%%
g_backg_df = m2_df.groupby(['gwas_category_count']).count().reset_index()[['gwas_category_count', 'rsid']]

#%%
g_target_df = g_target_df.merge(g_backg_df, on='gwas_category_count')
g_target_df.rename({'rsid_x': 'fisher_a_b', 'rsid_y': 'fisher_c_d'}, axis=1, inplace=1)
g_target_df['fisher_a_b'] = g_target_df['fisher_c_d'] - g_target_df['fisher_a_b']

#%%
pleio_1_target_df = g_target_df.loc[g_target_df['gwas_category_count'] == 1]
pleio_1_target_df = pleio_1_target_df[['etissue_category', 'fisher_a_b', 'fisher_c_d']].rename({'fisher_a_b': 'fisher_b', 'fisher_c_d': 'fisher_d'}, axis=1)

#%%
pleio_n_target_df = g_target_df.loc[g_target_df['gwas_category_count'] > 1]
pleio_n_target_df = pleio_n_target_df[['gwas_category_count', 'etissue_category', 'fisher_a_b', 'fisher_c_d']].rename({'fisher_a_b': 'fisher_a', 'fisher_c_d': 'fisher_c'}, axis=1)

#%%
fisher_df = pleio_n_target_df.merge(pleio_1_target_df, on='etissue_category')
oddsr_pval = fisher_df.apply(lambda x: fisher_exact(numpy.array([[x['fisher_a'], x['fisher_b']], [x['fisher_c'], x['fisher_d']]])), axis=1)
oddsr_pval_df = pandas.DataFrame(oddsr_pval.tolist(), columns=['oddsr', 'pvalue'])

#%%
fisher_df = pandas.concat([fisher_df, oddsr_pval_df], axis=1)
fisher_df.sort_values(by='pvalue', ascending=True, inplace=True)

#%%
fisher_df.to_csv(tsv_path, sep='\t', index=False)
