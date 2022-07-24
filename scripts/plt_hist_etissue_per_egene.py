from eqtl2gwas_pleiotropy.PathManager import PathManager

import os
import pandas
import pathlib
import sys
import matplotlib.pyplot as plt
from eqtl2gwas_pleiotropy.constants import label_fontsize, dpi

plt.rcParams["figure.figsize"] = (8, 6)
ylabel = "Prob. Density"

#%%
help_cmd_str = "todo"
try:
    h4_annot_tsv_path = sys.argv[1]
    count_per_rsid_gwas_tsv_path = sys.argv[2]
    upper_var_gwas_cat_count = int(sys.argv[3])
    hist_density_True_distr_egene_per_rsid_p2_png_path = sys.argv[4]
    if len(sys.argv) > 5:
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

outdir_path = os.path.dirname(hist_density_True_distr_egene_per_rsid_p2_png_path)
pathlib.Path(outdir_path).mkdir(parents=True, exist_ok=True)

#%%
h4_df = pandas.read_csv(h4_annot_tsv_path, sep="\t")

#%%
count_per_rsid_gwas_df = pandas.read_csv(count_per_rsid_gwas_tsv_path, sep="\t")
gwas_category_count_max_int = count_per_rsid_gwas_df['gwas_category_count'].max()

#%%
m_df = h4_df.merge(count_per_rsid_gwas_df, on=['chrom', 'pos', 'rsid'])

# %%
density = True
legend_label_ref = "GWAS cat. nb. 1"
sel_cols = ['egene', 'etissue_category']  # tissue per egene
title = "eTissues per eGene"
xlabel = "# eTissues"
ylabel = "Probability Density"
ylim = [0, 0.6]

#%% set upper_var_gwas_cat_count
m_df = m_df[sel_cols + ['gwas_category_count']]
m_df.loc[m_df['gwas_category_count'] >= upper_var_gwas_cat_count, "gwas_category_count"] = upper_var_gwas_cat_count

#%%
pleio1_df = m_df.loc[m_df['gwas_category_count'] == 1, sel_cols].drop_duplicates()
count_1_lst = pleio1_df.groupby(sel_cols[0]).count()[sel_cols[1]].tolist()

for p_count in range(2, upper_var_gwas_cat_count+1):

    pleio_df = m_df.loc[m_df['gwas_category_count'] == p_count, sel_cols].drop_duplicates()
    legend_label_pleio = "GWAS cat. nb. {}".format(p_count)

    count_i_lst = pleio_df.groupby(sel_cols[0]).count()[sel_cols[1]].tolist()
    bins = range(max(count_1_lst)+1)
    plt.hist(count_i_lst, alpha=0.5, density=density, label=legend_label_pleio, bins=range(22))  # density=False would make counts
    plt.hist(count_1_lst, alpha=0.5, density=density, label=legend_label_ref, bins=range(22))  # density=False would make counts

    plt.grid(axis='y')
    plt.legend(loc='upper right', fontsize=label_fontsize)
    plt.title(title, fontsize=label_fontsize)
    plt.xlabel(xlabel, fontsize=label_fontsize)
    plt.ylabel(ylabel, fontsize=label_fontsize)
    plt.ylim(ylim)

    png_path = os.path.join(outdir_path, 'hist_density_{}_distr_{}_per_{}_p{}.png'.format(density, sel_cols[1], sel_cols[0], p_count))
    plt.savefig(png_path, dpi=dpi)
    plt.clf()
    plt.close()

