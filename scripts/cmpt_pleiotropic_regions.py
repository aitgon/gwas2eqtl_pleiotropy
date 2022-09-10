from gwas2eqtl_pleiotropy.constants import region_bin, label_fontsize, tick_fontsize
from gwas2eqtl_pleiotropy.constants import seaborn_theme_dic
from matplotlib import pyplot as plt

import math
import os
import pandas
import pathlib
import seaborn
import sys


#%%
seaborn.set_theme(**seaborn_theme_dic)

#%%
help_cmd_str = "todo"
try:
    count_per_rsid_gwas_tsv_path = sys.argv[1]
    pleio_tsv_path = sys.argv[2]
    png_path = sys.argv[3]
    if len(sys.argv) > 4:
        print("""Two many arguments!
        {}""".format(help_cmd_str))
        sys.exit(1)
except IndexError:
    print("""Argument missing!
    {}""".format(help_cmd_str))
    sys.exit(1)
outdir_path = os.path.dirname(pleio_tsv_path)
pathlib.Path(outdir_path).mkdir(parents=True, exist_ok=True)

#%%
# count_per_rsid_gwas_tsv_path = os.path.join(PathManager.get_outdir_path(), "cmpt_count_per_rsid.py", "count_per_rsid_gwas.tsv")
df = pandas.read_csv(count_per_rsid_gwas_tsv_path, sep="\t")

#%% Do it forward
df.sort_values(by=['chrom', 'pos'], inplace=True, ascending=[True, True])
more_than_1 = 0
chrom_pleio_latest = 0
pos_pleio_latest = -99999
df['region_pleio_fwd'] = False

#%%
for i, row in df.iterrows():
    chrom = row['chrom']
    pos = row['pos']
    # import pdb; pdb.set_trace()
    if row['gwas_category_count'] > 1:
        chrom_pleio_latest = chrom
        pos_pleio_latest = pos
    if chrom == chrom_pleio_latest and (pos - pos_pleio_latest) <= region_bin:
        df.loc[i, 'region_pleio_fwd'] = True
df.to_csv("df.tsv", sep="\t", index=False)

#%% Do it reversed
df.sort_values(by=['chrom', 'pos'], inplace=True, ascending=[True, False])
more_than_1 = 0
chrom_pleio_latest = 0
pos_pleio_latest = 999999999
df['region_pleio_rev'] = False

#%%
for i, row in df.iterrows():
    chrom = row['chrom']
    pos = row['pos']
    df.loc[i, 'region_pleio'] = False
    if row['gwas_category_count'] > 1:  # update pos_pleio_latest
        chrom_pleio_latest = chrom
        pos_pleio_latest = pos
    if chrom == chrom_pleio_latest and (pos_pleio_latest - pos) <= region_bin:
        if df.loc[i, 'region_pleio_fwd']:
            df.loc[i, 'region_pleio'] = True

#%%
df.drop(['region_pleio_fwd', 'region_pleio_rev'], axis=1, inplace=True)
df.sort_values(by=['chrom', 'pos'], inplace=True, ascending=[True, True])

#%% bed file
region_lst = []
pos_prev = 0
region_pleio_prev = False
start = math.nan
end = math.nan
gwas_category_count = 0
gwas_category_lst = math.nan
category_lst = []

for i, row in df.iterrows():
    # beginning of region, set start, start category list, store category
    if row['region_pleio'] and not region_pleio_prev:
        cytoband = row['cytoband']
        start = row['pos']
        gwas_category_count = row['gwas_category_count']
        gwas_category_lst = row['gwas_category_lst']
        category_lst = row['gwas_category_lst'].split(",")
    # end of region, set end, store category
    elif not row['region_pleio'] and region_pleio_prev:
        end = pos_prev
    # middle of region, store categories
    if row['region_pleio'] and row['gwas_category_count'] > gwas_category_count:
        gwas_category_count = row['gwas_category_count']
        gwas_category_lst = row['gwas_category_lst']
        category_lst = category_lst + row['gwas_category_lst'].split(",")
    # reset start and end, store region
    if not math.isnan(start) and not math.isnan(end):
        # region_lst.append([row['chrom'], start, end, gwas_subcategory_count, gwas_subcategory_lst, category_lst])
        # import pdb; pdb.set_trace()
        category_lst = sorted([*set(category_lst)])
        category_str = ','.join(category_lst)
        region_lst.append([row['chrom'], cytoband, start, end, len(category_lst), category_str])
        start = math.nan
        end = math.nan
        gwas_subcategory_count = 0
    pos_prev = row['pos']
    region_pleio_prev = row['region_pleio']

#%% tsv
regions_pleio_df = pandas.DataFrame(region_lst, columns=['chrom', 'cytoband', 'start', 'end', 'gwas_category_count', 'gwas_category_lst'])
regions_pleio_df.to_csv(pleio_tsv_path, sep="\t", index=False, header=True)

#%% bed
regions_pleio_df['start'] = regions_pleio_df['start'] - 1
regions_pleio_df['chrom'] = 'chr' + regions_pleio_df['chrom'].astype('str')
regions_pleio_df.drop(['cytoband'], axis=1, inplace=True)
pleio_bed_path = os.path.join(outdir_path, "region_window_{}.bed".format(region_bin))
regions_pleio_df.to_csv(pleio_bed_path, sep="\t", index=False, header=False)

#%########################################### bed files
for count_pleio in range(1, 6):
    region_pleio_i_bed_path = os.path.join(outdir_path, "region_window_{}_pleio_{}.bed".format(region_bin, count_pleio))
    region_pleio_i_df = regions_pleio_df.loc[regions_pleio_df['gwas_category_count'] == count_pleio,]
    region_pleio_i_df.to_csv(region_pleio_i_bed_path, sep="\t", index=False, header=False)


#%##############
# histogram
plt.rcParams["figure.figsize"] = (8, 6)
region_lenght_ser = (regions_pleio_df['end'] - regions_pleio_df['start'])

ylabel = "Proportion"
title = "Lengths of pleiotropic regions"
ylim=[1e-10, 1]
edgecolor='k'
linewidth = 2

hist_kwargs = {'density': False, 'edgecolor': edgecolor, 'linewidth': linewidth}

#%%
data_ser = region_lenght_ser/100000
seaborn.histplot(data_ser, stat='percent', discrete=True)

plt.grid(visible=True, axis='y')
plt.title(title, fontsize=label_fontsize)
plt.xlabel("Region length [1e5 bp]", fontsize=label_fontsize)
plt.xticks(fontsize=tick_fontsize)
plt.ylabel(ylabel, fontsize=label_fontsize)
# plt.yscale('log')
plt.yticks(fontsize=tick_fontsize)

plt.tight_layout()
plt.savefig(png_path)
plt.close()
