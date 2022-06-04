from eqtl2gwas_pleiotropy.PathManager import PathManager
from eqtl2gwas_pleiotropy.constants import region_bin
from matplotlib import pyplot as plt

import math
import numpy
import os
import pandas
import pathlib


#%% Outdir
if not '__file__' in locals():
    __file__ = "cmpt_pleiotropic_regions.py"
outdir_path = os.path.join(PathManager.get_outdir_path(), os.path.basename(__file__))
pathlib.Path(outdir_path).mkdir(parents=True, exist_ok=True)

#%%
count_per_rsid_gwas_tsv_path = os.path.join(PathManager.get_outdir_path(), "cmpt_count_per_rsid.py", "count_per_rsid_gwas.tsv")
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
    if row['gwas_subcategory_count'] > 1:
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
    if row['gwas_subcategory_count'] > 1:  # update pos_pleio_latest
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
gwas_subcategory_count = 0
gwas_subcategory_lst = math.nan
category_lst = []

for i, row in df.iterrows():
    # beginning of region, set start, start category list, store category
    if row['region_pleio'] and not region_pleio_prev:
        start = row['pos']
        gwas_subcategory_count = row['gwas_subcategory_count']
        gwas_subcategory_lst = row['gwas_subcategory_lst']
        category_lst = row['gwas_subcategory_lst'].split(",")
    # end of region, set end, store category
    elif not row['region_pleio'] and region_pleio_prev:
        end = pos_prev
    # middle of region, store categories
    if row['region_pleio'] and row['gwas_subcategory_count'] > gwas_subcategory_count:
        gwas_subcategory_count = row['gwas_subcategory_count']
        gwas_subcategory_lst = row['gwas_subcategory_lst']
        category_lst = category_lst + row['gwas_subcategory_lst'].split(",")
    # reset start and end, store region
    if not math.isnan(start) and not math.isnan(end):
        # region_lst.append([row['chrom'], start, end, gwas_subcategory_count, gwas_subcategory_lst, category_lst])
        # import pdb; pdb.set_trace()
        category_lst = sorted([*set(category_lst)])
        category_str = ','.join(category_lst)
        region_lst.append([row['chrom'], start, end, len(category_lst), category_str])
        start = math.nan
        end = math.nan
        gwas_subcategory_count = 0
    pos_prev = row['pos']
    region_pleio_prev = row['region_pleio']

# tsv
regions_pleio_df = pandas.DataFrame(region_lst, columns=['chrom', 'start', 'end', 'gwas_category_count', 'gwas_category_lst'])
pleio_tsv_path = os.path.join(outdir_path, "region_window_{}.tsv".format(region_bin))
regions_pleio_df.to_csv(pleio_tsv_path, sep="\t", index=False, header=True)

# bed
regions_pleio_df['start'] = regions_pleio_df['start'] - 1
regions_pleio_df['chrom'] = 'chr' + regions_pleio_df['chrom'].astype('str')
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

ylabel = "# Regions"
title = "Length Distribution of Pleiotropic Regions"
ylim=[1e-10, 1]
edgecolor='k'
label_fontsize = 20
tick_fontsize = 10
linewidth = 2

hist_kwargs = {'density': False, 'edgecolor': edgecolor, 'linewidth': linewidth}

#%%
bins = numpy.array(range(11))*100000
ax = region_lenght_ser.hist(**hist_kwargs, bins=bins)

ax.set_xlabel("Region Length [bp]", fontsize=label_fontsize)
ax.set_ylabel(ylabel, fontsize=label_fontsize)
ax.set_yscale('log')
plt.title(title, fontsize=label_fontsize)
plt.xticks(fontsize=tick_fontsize)
plt.yticks(fontsize=tick_fontsize)

fig = ax.get_figure()
png_path = os.path.join(outdir_path, "regions_{}_length_hist.png".format(region_bin))
fig.savefig(png_path)
plt.clf()
plt.close()
