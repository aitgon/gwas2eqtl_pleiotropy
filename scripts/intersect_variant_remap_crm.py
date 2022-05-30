from eqtl2gwas_pleiotropy.PathManager import PathManager

import os
import pandas
import pathlib
import sys
import matplotlib.pyplot as plt

from eqtl2gwas_pleiotropy.ReMapCRM import ReMapCRM
from eqtl2gwas_pleiotropy.URL import URL
from eqtl2gwas_pleiotropy.constants import label_fontsize

#%%
plt.rcParams["figure.figsize"] = (8, 6)
ylabel = "Prob. Density"

#%% Input
variant_pleio_flank_0_bed_path = os.path.join(PathManager.get_outdir_path(), "cmpt_count_per_rsid.py", "variant_pleio_{}_flank_0.bed")
# if not os.path.isfile(variant_pleio_flank_0_bed_path):
#     print("input file does not exit")
#     sys.exit(1)

#%% Output
if not '__file__' in locals():
    outdir_path = os.path.join(PathManager.get_outdir_path(), "intersect_variant_remap_crm.py")
else:
    outdir_path = os.path.join(PathManager.get_outdir_path(), os.path.basename(__file__))
pathlib.Path(outdir_path).mkdir(parents=True, exist_ok=True)

out_tsv_path = os.path.join(outdir_path, "todo.tsv")

#%%
# count_per_rsid_gwas_df = pandas.read_csv(variant_pleio_flank_0_bed_path, sep="\t")

#%%
# import requests
# # Define the remote file to retrieve
# remote_url = 'https://www.google.com/robots.txt'
# # Define the local filename to save data
# local_file = 'local_copy.txt'
# # Make http request for remote file data
# data = requests.get(remote_url)
# # Save file data to local copy
# with open(local_file, 'wb')as file:
# file.write(data.content)

remap_crm_path = URL("http://remap.univ-amu.fr/storage/remap2022/hg38/MACS2/remap2022_crm_macs2_hg38_v1_0.bed.gz", data_public_dir="/home/gonzalez/Software/public").download()

#%%
# ReMapCRM().download()
remap_crm_path = "/home/gonzalez/Software2/public/remap.univ-amu.fr/storage/remap2022/hg38/MACS2/remap2022_crm_macs2_hg38_v1_0.bed.gz"

#%%
# m_df = h4_df.merge(count_per_rsid_gwas_df, on=['chrom', 'pos', 'rsid'])

###############################################################################
# Hypothesis 1: How many gwas per egene
###############################################################################

# # %%
# sel_cols = ['rsid', 'egene']  # egene per variant
# # sel_cols = ['gwas_trait_name', 'egene']
# ylabel = "Prob. Density"
# density = False
#
# #%% segregate gwas
# segregated_df = m_df[[sel_cols[0], 'gwas_subcategory_count']].sort_values(by='gwas_subcategory_count', ascending=False)
# segregated_df = segregated_df.drop_duplicates(sel_cols[0], keep='first')
#
# #%%
# distr_back_egene_to_gwas_lst = m_df[sel_cols].drop_duplicates().groupby([sel_cols[0]]).count()[sel_cols[1]].to_list()  # background
# for density in [False, True]:
#     for p_count in range(1, 6):
#         title = "Distrib. eGene per variant - {} Categories".format(p_count)
#         this_variable_p_df = segregated_df.loc[segregated_df['gwas_subcategory_count'] == p_count, ]
#         selected_df = m_df.merge(this_variable_p_df, on=sel_cols[0])[sel_cols].drop_duplicates()
#         count_gwas_per_egene_lst = selected_df.groupby([sel_cols[0]]).count()[sel_cols[1]].to_list()
#         bins = range(max(distr_back_egene_to_gwas_lst)+1)
#         plt.title(title, fontsize=label_fontsize)
#         plt.xlabel("# {}".format(sel_cols[1]), fontsize=label_fontsize)
#         if density:
#             plt.ylim([0, 1])
#             plt.ylabel("Probability {}".format(sel_cols[0]), fontsize=label_fontsize)
#         else:
#             plt.ylabel("# {}".format(sel_cols[0]), fontsize=label_fontsize)
#         plt.hist(count_gwas_per_egene_lst, alpha=0.5, density=density, label='{} pleio {}'.format(sel_cols[0], p_count), bins=range(22))  # density=False would make counts
#         plt.hist(distr_back_egene_to_gwas_lst, alpha=0.5, density=density, label='{} all'.format(sel_cols[0]), bins=range(22))  # density=False would make counts
#         plt.grid(axis='y')
#         plt.legend(loc='upper right', fontsize=label_fontsize)
#         png_path = os.path.join(outdir_path, 'hist_density_{}_distr_{}_per_{}_p{}.png'.format(density, sel_cols[1], sel_cols[0], p_count))
#         plt.savefig(png_path)
#         plt.clf()
#         plt.close()
