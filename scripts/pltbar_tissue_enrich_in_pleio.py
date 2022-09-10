"""This script will compute tissue enrichment in pleio n vs pleio 1 using fisher exact
FDR 0.05"""

import os
import numpy
import pandas
import pathlib
import seaborn
import sys

from scipy.stats import fisher_exact
from statsmodels.stats.multitest import fdrcorrection
from matplotlib import pyplot as plt
from gwas2eqtl_pleiotropy.constants import dpi, label_fontsize, tick_fontsize

#%%
plt.rcParams["figure.figsize"] = (16, 6)

help_cmd_str = "todo"
try:
    gwas_count_etissue_tsv_path = sys.argv[1]
    upper_var_gwas_cat_count = int(sys.argv[2])
    tsv_path = sys.argv[3]
    if len(sys.argv) > 4:
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
# h4_df = pandas.read_csv(h4_tsv_path, sep="\t")

#%%
df = pandas.read_csv(gwas_count_etissue_tsv_path, sep="\t")

#%%
df = df[['rsid', 'gwas_category_count', 'etissue_subcategory_lst']].drop_duplicates()
df['etissue_subcategory_lst'] = df['etissue_subcategory_lst'].str.split(',')
df = df.explode('etissue_subcategory_lst')
df.rename({'etissue_subcategory_lst': 'etissue_category'}, axis=1, inplace=True)

#%%
# df = m_df[['rsid', 'gwas_category_count', 'etissue_category']].drop_duplicates().sort_values('gwas_category_count', ascending=False)
df.loc[df['gwas_category_count'] >= upper_var_gwas_cat_count, 'gwas_category_count'] = upper_var_gwas_cat_count

#%%
g_target_df = df.groupby(['gwas_category_count', 'etissue_category']).count().reset_index()

#%%
# g_backg_df = df.groupby(['gwas_category_count']).count().reset_index()[['gwas_category_count', 'rsid']]
g_backg_df = df[['rsid', 'gwas_category_count']].drop_duplicates().groupby(['gwas_category_count']).count()

#%%
g_target_df = g_target_df.merge(g_backg_df, on='gwas_category_count')

#%% ########################################################################
# g_target_df = g_target_df.loc[g_target_df['etissue_category']=="ImmuneCell"]

#%%
g_target_df.rename({'rsid_x': 'fisher_a_b', 'rsid_y': 'fisher_c_d'}, axis=1, inplace=1)
g_target_df['fisher_c_d'] = g_target_df['fisher_c_d'] - g_target_df['fisher_a_b']

#%%
pleio_1_target_df = g_target_df.loc[g_target_df['gwas_category_count'] == 1]
pleio_1_target_df = pleio_1_target_df[['etissue_category', 'fisher_a_b', 'fisher_c_d']].rename({'fisher_a_b': 'fisher_b', 'fisher_c_d': 'fisher_d'}, axis=1)

#%%
pleio_n_target_df = g_target_df.loc[g_target_df['gwas_category_count'] > 1]
pleio_n_target_df = pleio_n_target_df[['gwas_category_count', 'etissue_category', 'fisher_a_b', 'fisher_c_d']].rename({'fisher_a_b': 'fisher_a', 'fisher_c_d': 'fisher_c'}, axis=1)

#%%
fisher_df = pleio_n_target_df.merge(pleio_1_target_df, on='etissue_category')
oddsr_pval = fisher_df.apply(lambda x: fisher_exact(numpy.array([[x['fisher_a'], x['fisher_b']], [x['fisher_c'], x['fisher_d']]])), axis=1)
oddsr_pval_df = pandas.DataFrame(oddsr_pval.tolist(), columns=['oddsratio', 'p_value'])

#%%
fisher_df = pandas.concat([fisher_df, oddsr_pval_df], axis=1)
fisher_df.sort_values(by=['etissue_category', 'gwas_category_count'], ascending=True, inplace=True)

#%%
fisher_df['fdr5'] = fdrcorrection(fisher_df['p_value'])[1]


#%%
fisher_df.to_csv(tsv_path, sep='\t', index=False)

#%% plot
g = seaborn.catplot(height=8.27, aspect=6/8.27, data=fisher_df, kind="bar", x="oddsratio", y="etissue_category", hue="gwas_category_count", ci="sd", palette="rocket_r", orient='h')

label_fontsize=14
tick_fontsize=14
plt.title("eTissue enrichment in variants", fontsize=label_fontsize)
plt.xlabel("Odds ratio - category count k vs 1", fontsize=label_fontsize)
plt.xticks(fontsize=tick_fontsize, rotation=0)
plt.ylabel("eTissue category", fontsize=label_fontsize)
plt.yticks(fontsize=tick_fontsize)
plt.grid(axis='x')

plt.tight_layout()
png_path = os.path.join(outdir_path, "plt.png")
plt.savefig(png_path, dpi=dpi)
plt.close()
