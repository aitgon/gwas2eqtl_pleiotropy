from eqtl2gwas_pleiotropy.PathManager import PathManager

import os
import pandas
import pathlib
import sys
import matplotlib.pyplot as plt
from eqtl2gwas_pleiotropy.constants import label_fontsize

plt.rcParams["figure.figsize"] = (8, 6)
ylabel = "Prob. Density"

#%% Input1


h4_annot_tsv_path = os.path.join(PathManager.get_outdir_path(), "annotate_h4.py", "h4_annotated.tsv")
if not os.path.isfile(h4_annot_tsv_path):
    print("input file does not exit")
    sys.exit(1)

#%% Input2
count_per_rsid_gwas_tsv_path = os.path.join(PathManager.get_outdir_path(), "cmpt_count_per_rsid.py", "count_per_rsid_gwas.tsv")
if not os.path.isfile(count_per_rsid_gwas_tsv_path):
    print("input file does not exit")
    sys.exit(1)

#%% Output
if not '__file__' in locals():
    outdir_path = os.path.join(PathManager.get_outdir_path(), "plt_hist_egenes_per_gwas.py")
else:
    outdir_path = os.path.join(PathManager.get_outdir_path(), os.path.basename(__file__))
pathlib.Path(outdir_path).mkdir(parents=True, exist_ok=True)

out_tsv_path = os.path.join(outdir_path, "todo.tsv")

#%%
h4_df = pandas.read_csv(h4_annot_tsv_path, sep="\t")

#%%
count_per_rsid_gwas_df = pandas.read_csv(count_per_rsid_gwas_tsv_path, sep="\t")

#%%
m_df = h4_df.merge(count_per_rsid_gwas_df, on=['chrom', 'pos', 'rsid'])

###############################################################################
# Hypothesis 1: How many gwas per egene
###############################################################################

# %%
sel_cols = ['egene', 'gwas_trait_name']  # gwas per egene
# sel_cols = ['gwas_trait_name', 'egene']
ylabel = "Prob. Density"
density = False

#%% segregate gwas
segregated_df = m_df[[sel_cols[0], 'gwas_subcategory_count']].sort_values(by='gwas_subcategory_count', ascending=False)
segregated_df = segregated_df.drop_duplicates(sel_cols[0], keep='first')

#%%
label_fontsize = 30
tick_fontsize = 24
xlabel = "# GWAS traits"
ylabel = "Probability density"

distr_back_egene_to_gwas_lst = m_df[sel_cols].drop_duplicates().groupby([sel_cols[0]]).count()[sel_cols[1]].to_list()  # background
for density in [False, True]:
    for p_count in range(1, 6):
        title = "GWAS per eGene - Pleiotropy {}".format(p_count)
        this_variable_p_df = segregated_df.loc[segregated_df['gwas_subcategory_count'] == p_count, ]
        selected_df = m_df.merge(this_variable_p_df, on=sel_cols[0])[sel_cols].drop_duplicates()
        count_gwas_per_egene_lst = selected_df.groupby([sel_cols[0]]).count()[sel_cols[1]].to_list()
        bins = range(max(distr_back_egene_to_gwas_lst)+1)
        plt.title(title, fontsize=label_fontsize)
        plt.xlabel(xlabel, fontsize=label_fontsize)
        if density:
            plt.ylim([0, 0.65])
            plt.ylabel(ylabel, fontsize=label_fontsize)
        else:
            plt.ylabel("# {}".format(sel_cols[0]), fontsize=label_fontsize)
        plt.hist(count_gwas_per_egene_lst, alpha=0.5, density=density, label='eGenes: Pleio. {}'.format(p_count), bins=range(22))  # density=False would make counts
        plt.hist(distr_back_egene_to_gwas_lst, alpha=0.5, density=density, label='All eGenes', bins=range(22))  # density=False would make counts
        plt.grid(axis='y')
        plt.legend(loc='upper right', fontsize=label_fontsize)
        png_path = os.path.join(outdir_path, 'hist_density_{}_distr_{}_per_{}_p{}.png'.format(density, sel_cols[1], sel_cols[0], p_count))
        if density:
            plt.savefig(png_path, dpi=600)
        else:
            plt.savefig(png_path)
        plt.clf()
        plt.close()
