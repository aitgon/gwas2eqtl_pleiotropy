from statannot import add_stat_annotation

from gwas2eqtl_pleiotropy.UCSC import UCSC
from gwas2eqtl_pleiotropy.constants import label_fontsize, tick_fontsize

import matplotlib.pyplot as plt
import os
import pandas
import pathlib
import seaborn
import sys


# Plot parameters
plt.rcParams["figure.figsize"] = (8, 6)
from gwas2eqtl_pleiotropy.constants import seaborn_theme_dic
seaborn.set_theme(**seaborn_theme_dic)

#%%
help_cmd_str = "todo"
try:
    h4_annot_tsv_path = sys.argv[1]
    count_per_rsid_gwas_tsv_path = sys.argv[2]
    interactome_tsv_path = sys.argv[3]
    upper_var_gwas_cat_count = int(sys.argv[4])
    vlnplt_png_path = sys.argv[5]
    if len(sys.argv) > 6:
        print("""Two many arguments!
        {}""".format(help_cmd_str))
        sys.exit(1)
except IndexError:
    print("""Argument missing!
    {}""".format(help_cmd_str))
    sys.exit(1)

#%% Input1
if not os.path.isfile(h4_annot_tsv_path):
    print("input file does not exit")
    sys.exit(1)

#%% Input2
if not os.path.isfile(count_per_rsid_gwas_tsv_path):
    print("input file does not exit")
    sys.exit(1)

outdir_path = os.path.dirname(vlnplt_png_path)
pathlib.Path(outdir_path).mkdir(parents=True, exist_ok=True)

#%%
h4_df = pandas.read_csv(h4_annot_tsv_path, sep="\t")

#%%
count_per_rsid_gwas_df = pandas.read_csv(count_per_rsid_gwas_tsv_path, sep="\t")
gwas_category_count_max_int = count_per_rsid_gwas_df['gwas_category_count'].max()

#%%
m_df = h4_df.merge(count_per_rsid_gwas_df, on=['chrom', 'pos', 'rsid'])
m_df = m_df[['chrom', 'pos', 'rsid', 'egene', 'egene_symbol', 'gwas_category_count']].drop_duplicates()


#%%
ensg2enst2ensp_df = UCSC(database='hg18').get_ensg2enst2ensp()
m_df = m_df.merge(ensg2enst2ensp_df[['gene', 'protein']], left_on='egene', right_on='gene').drop_duplicates()

#%%
score_cutoff=700
inter_df = pandas.read_csv(interactome_tsv_path, sep=" ", compression='gzip')
inter_df = inter_df.loc[inter_df['combined_score'] >= score_cutoff, ]
inter_df['protein1']=inter_df['protein1'].str.replace('9606.', '', regex=False)
inter_df['protein2']=inter_df['protein2'].str.replace('9606.', '', regex=False)
inter_df = pandas.concat([inter_df, inter_df[['protein2', 'protein1', 'combined_score']]], axis=0).drop_duplicates()

#%%
m_df = m_df.merge(inter_df, left_on="protein", right_on='protein1')
m_df = m_df[['chrom', 'pos', 'rsid', 'egene', 'protein2', 'gwas_category_count']].drop_duplicates()

# %%
sel_cols = ['egene', 'protein2']

# %%
m2_df = m_df[['chrom', 'pos', 'rsid'] + sel_cols + ['gwas_category_count']].drop_duplicates()
m2_df.sort_values(['gwas_category_count', 'chrom', 'pos'], inplace=True, ascending=[False, True, True])
tsv_path = os.path.join(outdir_path, 'variants2egenes.tsv')
m2_df.to_csv(tsv_path, header=True, index=False, sep='\t')

#%% set upper_var_gwas_cat_count
m_df = m_df[sel_cols + ['gwas_category_count']]
m_df.loc[m_df['gwas_category_count'] >= upper_var_gwas_cat_count, "gwas_category_count"] = upper_var_gwas_cat_count

#%%
m_df = m_df.drop_duplicates()
m_df = m_df.groupby(['egene', 'gwas_category_count']).count()
m_df = m_df.reset_index()
m_df.columns = ['egene', 'gwas_category_count', 'interactor_count']

#%%
describe_tsv_path = os.path.join(outdir_path, "describe.tsv")
m_df.groupby('gwas_category_count')['interactor_count'].apply(lambda x: x.describe()).to_csv(describe_tsv_path, sep="\t")

#%%
order = [*range(1, upper_var_gwas_cat_count+1)]
xticklabels = order.copy()
xticklabels[-1] = '≥{}'.format(order[-1])
title = "Interactors per egene"
xlabel = "GWAS category count"
ylabel = "Phys. interactor count"
y = "interactor_count"

#%%
box_pairs = [(1, i) for i in range(2, upper_var_gwas_cat_count+1) ]
ax = seaborn.violinplot(x="gwas_category_count", y=y, data=m_df, order=order, palette="rocket_r")
test_results = add_stat_annotation(ax, data=m_df, x="gwas_category_count", y=y, order=order,
                                   box_pairs=box_pairs,
                                   test='Mann-Whitney', text_format='star',
                                   loc='inside', verbose=2)

plt.title(title, fontsize=label_fontsize)
plt.xlabel(xlabel, fontsize=label_fontsize)
plt.xticks(fontsize=tick_fontsize, rotation=0)
plt.ylabel(ylabel, fontsize=label_fontsize)
plt.yticks(fontsize=tick_fontsize)
ax.set_xticklabels(xticklabels)

plt.tight_layout()
plt.savefig(vlnplt_png_path)
plt.close()
