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
    wlength = int(sys.argv[1])
    count_per_rsid_gwas_tsv_path = sys.argv[2]
    pleio_tsv_path = sys.argv[3]
    png_path = sys.argv[4]
    if len(sys.argv) > 5:
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
df = pandas.read_csv(count_per_rsid_gwas_tsv_path, sep="\t")

# #%% Do it forward
# df.sort_values(by=['chrom', 'pos38'], inplace=True, ascending=[True, True])
# more_than_1 = 0
# chrom_pleio_latest = 0
# pos_pleio_latest = -99999
# df['region_pleio_fwd'] = False
#
# #%%
# for i, row in df.iterrows():
#     chrom = row['chrom']
#     pos38 = row['pos38']
#     # import pdb; pdb.set_trace()
#     if row['gwas_class_count'] > 1:
#         chrom_pleio_latest = chrom
#         pos_pleio_latest = pos38
#     if chrom == chrom_pleio_latest and (pos38 - pos_pleio_latest) <= region_bin:
#         df.loc[i, 'region_pleio_fwd'] = True
#
# #%% Do it reversed
# df.sort_values(by=['chrom', 'pos38'], inplace=True, ascending=[True, False])
# more_than_1 = 0
# chrom_pleio_latest = 0
# pos_pleio_latest = 999999999
# df['region_pleio_rev'] = False
#
# #%%
# for i, row in df.iterrows():
#     chrom = row['chrom']
#     pos38 = row['pos38']
#     df.loc[i, 'region_pleio'] = False
#     if row['gwas_class_count'] > 1:  # update pos_pleio_latest
#         chrom_pleio_latest = chrom
#         pos_pleio_latest = pos38
#     if chrom == chrom_pleio_latest and (pos_pleio_latest - pos38) <= region_bin:
#         if df.loc[i, 'region_pleio_fwd']:
#             df.loc[i, 'region_pleio'] = True

#%%
# df.drop(['region_pleio_fwd', 'region_pleio_rev'], axis=1, inplace=True)
# df.sort_values(by=['chrom', 'pos38'], inplace=True, ascending=[True, True])

pleio_df = pandas.DataFrame()
df.sort_values(['chrom', 'pos38'], inplace=True)
df['gwas_class_lst'] = df['gwas_class_lst'].str.split(',')
df['gwas_class_count'] = 1

########################
# window = 1e3
# window = 3e3
# window = 1e4
# window = 3e4
# window = 1e5
for chrom in sorted(df['chrom'].unique()):
    print(chrom)
    df_pleio_chrom = df.loc[df['chrom'] == chrom].copy()
    df_pleio_chrom['dist'] = [0] + [df_pleio_chrom.iloc[i]['pos38'] - df_pleio_chrom.iloc[i - 1]['pos38'] for i in range(1, df_pleio_chrom.shape[0])]
    region_pleio = 0
    df_pleio_chrom['region_pleio'] = region_pleio
    for i, row in df_pleio_chrom.iloc[1:].iterrows():
        if row['dist'] <= wlength:
            df_pleio_chrom.loc[i, 'region_pleio'] = region_pleio
        else:
            region_pleio = region_pleio + 1
            df_pleio_chrom.loc[i, 'region_pleio'] = region_pleio

    df_pleio_chrom.drop('dist', inplace=True, axis=1)
    df_pleio_chrom_grouped = df_pleio_chrom.groupby(['chrom', 'region_pleio']).agg({'cytoband': lambda x: sorted(x.unique()),
                                                 'pos38': lambda x: str(min(x)) + "-" + str(max(x)),
                                                 'rsid': lambda x: sorted(x.unique()),
                                                 'gwas_class_lst': lambda x: sorted(set(x.sum()))
                                                           })
    df_pleio_chrom_grouped['gwas_class_count'] = df_pleio_chrom_grouped['gwas_class_lst'].apply(lambda x: len(x))
    pleio_df = pandas.concat([pleio_df, df_pleio_chrom_grouped], axis=0, verify_integrity=False)

pleio_df.reset_index(inplace=True)
pleio_df = pleio_df.sort_values(['gwas_class_count', 'chrom', 'pos38'], ascending=[False, True, True])
# import pdb; pdb.set_trace()
pleio_df[['start', 'end']] = pleio_df['pos38'].str.split('-', expand=True)
pleio_df['start'] = pleio_df['start'].astype(int)
pleio_df['end'] = pleio_df['end'].astype(int)
import pdb; pdb.set_trace()
pleio_df = pleio_df[['chrom', 'start', 'end', 'cytoband', 'rsid', 'gwas_class_count', 'gwas_class_lst']]
pleio_df.to_csv(pleio_tsv_path, sep="\t", index=False, header=True)


# ########################
# #%% bed file
# region_lst = []
# pos_prev = 0
# region_pleio_prev = False
# start = math.nan
# end = math.nan
# gwas_class_count = 0
# gwas_class_lst = math.nan
# category_lst = []
#
# for i, row in df.iterrows():
#     # beginning of region, set start, start category list, store category
#     if row['region_pleio'] and not region_pleio_prev:
#         cytoband = row['cytoband']
#         start = row['pos38']
#         gwas_class_count = row['gwas_class_count']
#         gwas_class_lst = row['gwas_class_lst']
#         category_lst = row['gwas_class_lst'].split(",")
#     # end of region, set end, store category
#     elif not row['region_pleio'] and region_pleio_prev:
#         end = pos_prev
#     # middle of region, store categories
#     if row['region_pleio'] and row['gwas_class_count'] > gwas_class_count:
#         gwas_class_count = row['gwas_class_count']
#         gwas_class_lst = row['gwas_class_lst']
#         category_lst = category_lst + row['gwas_class_lst'].split(",")
#     # reset start and end, store region
#     if not math.isnan(start) and not math.isnan(end):
#         category_lst = sorted([*set(category_lst)])
#         category_str = ','.join(category_lst)
#         region_lst.append([row['chrom'], cytoband, start, end, len(category_lst), category_str])
#         start = math.nan
#         end = math.nan
#         gwas_subcategory_count = 0
#     pos_prev = row['pos38']
#     region_pleio_prev = row['region_pleio']
#
# #%% tsv
# df = pandas.DataFrame(region_lst, columns=['chrom', 'cytoband', 'start', 'end', 'gwas_class_count', 'gwas_class_lst'])
# df.sort_values(['gwas_class_count', 'chrom', 'start'], ascending=[False, True, True], inplace=True)
# df.to_csv(pleio_tsv_path, sep="\t", index=False, header=True)
#
# #%##############
# #%% GWAS region pleiotropy for MS
# regions_pleio_ms_df = df.copy()
# regions_pleio_ms_df = regions_pleio_ms_df.sort_values(by=['gwas_class_count', 'chrom', 'start', 'gwas_class_lst'], ascending=[False, True, True, True])
# regions_pleio_ms_df = regions_pleio_ms_df.drop_duplicates('cytoband', keep='first')
# regions_pleio_ms_df = regions_pleio_ms_df.loc[regions_pleio_ms_df['gwas_class_count'] >= 6]
# # format output
# regions_pleio_ms_df.drop(['gwas_class_count'], inplace=True, axis=1)
# regions_pleio_ms_df['gwas_class_lst'] = regions_pleio_ms_df['gwas_class_lst'].str.replace(',', ', ')
# regions_pleio_ms_df['start'] = regions_pleio_ms_df['start'].apply(lambda x : '{0:,}'.format(x))
# regions_pleio_ms_df['end'] = regions_pleio_ms_df['end'].apply(lambda x : '{0:,}'.format(x))
# tsv_path = os.path.join(outdir_path, "region_window_ms_100000.tsv")
# regions_pleio_ms_df.to_csv(tsv_path, sep="\t", index=False)
#
# #%% bed
# df['start'] = df['start'] - 1
# df['chrom'] = 'chr' + df['chrom'].astype('str')
# df.drop(['cytoband'], axis=1, inplace=True)
# pleio_bed_path = os.path.join(outdir_path, "region_window_{}.bed".format(region_bin))
# df.to_csv(pleio_bed_path, sep="\t", index=False, header=False)
#
# #%########################################### bed files
# for count_pleio in range(1, 6):
#     region_pleio_i_bed_path = os.path.join(outdir_path, "region_window_{}_pleio_{}.bed".format(region_bin, count_pleio))
#     region_pleio_i_df = df.loc[df['gwas_class_count'] == count_pleio,]
#     region_pleio_i_df.to_csv(region_pleio_i_bed_path, sep="\t", index=False, header=False)

#############################################
# Distribution of region by length
# % region <= 100kb, % >100kb and <=200kb, ...
region_length_ser = (pleio_df['end'] - pleio_df['start'])
data_ser = region_length_ser / wlength

ax = plt.hist(data_ser, bins=range(50), density=True)
plt.xlim(0, 10)
plt.tight_layout()
plt.savefig(png_path)
plt.close()

#%##############
# histogram
plt.rcParams["figure.figsize"] = (8, 6)
region_length_ser = (pleio_df['end'] - pleio_df['start'])

ylabel = "Cumulative proportion"
title = "Length of pleiotropic regions"
ylim=[0, 1]
edgecolor='k'
linewidth = 2

hist_kwargs = {'density': False, 'edgecolor': edgecolor, 'linewidth': linewidth}

#%%
data_ser = region_length_ser / wlength

# import pdb; pdb.set_trace()
# seaborn.histplot(data_ser, stat='percent', discrete=True)
ax = plt.hist(data_ser, bins=range(50), density=True, cumulative=True)
plt.xlim(0, 10)

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
