import pdb
from itertools import chain

from gwas2eqtl_pleiotropy.constants import region_bin, label_fontsize, tick_fontsize
from gwas2eqtl_pleiotropy.constants import seaborn_theme_dic
from matplotlib import pyplot as plt


import os
import pandas
import pathlib
import seaborn
import sys
import ast


#%%
seaborn.set_theme(**seaborn_theme_dic)

#%%
help_cmd_str = "todo"
try:
    wlength = int(sys.argv[1])
    manuscript_pleio_cutoff = int(sys.argv[2])
    count_per_rsid_gwas_egene_etissue_ods = sys.argv[3]
    pleio_tsv_path = sys.argv[4]
    png_path = sys.argv[5]
    if len(sys.argv) > 6:
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
df = pandas.read_excel(count_per_rsid_gwas_egene_etissue_ods, engine='odf')
df['gwas_category_lst'] = df['gwas_category_lst'].apply(lambda x: ast.literal_eval(x))
df['egene_lst']= df['egene_lst'].apply(lambda x: ast.literal_eval(x))
df['eqtl_gene_symbol_lst'] = df['eqtl_gene_symbol_lst'].apply(lambda x: ast.literal_eval(x))
df['etissue_category_term_lst'] = df['etissue_category_term_lst'].apply(lambda x: ast.literal_eval(x))

region_df = pandas.DataFrame()
df.sort_values(['chrom', 'pos38'], inplace=True)
# df['gwas_category_lst'] = df['gwas_category_lst'].str.split(';')
df['gwas_category_count'] = 1

########################
# window = 1e3
# window = 3e3
# window = 1e4
# window = 3e4
# window = 1e5
for chrom in sorted(df['chrom'].unique()):  # problem split by chromosome
    df_chrom = df.loc[df['chrom'] == chrom].copy()

    # compute distance between consecutive eQTLs
    df_chrom['dist'] = [0] + [df_chrom.iloc[i]['pos38'] - df_chrom.iloc[i - 1]['pos38'] for i in range(1, df_chrom.shape[0])]
    region_pleio = 0  # init region pleio label
    df_chrom['region_pleio'] = region_pleio
    for i, row in df_chrom.iloc[1:].iterrows():  # Set or increase region label according to distance and wlength
        if row['dist'] <= wlength:
            df_chrom.loc[i, 'region_pleio'] = region_pleio
        else:
            region_pleio = region_pleio + 1
            df_chrom.loc[i, 'region_pleio'] = region_pleio

    df_chrom.drop('dist', inplace=True, axis=1)
    df_chrom_groupby = df_chrom.groupby(['chrom', 'region_pleio'])
    df_chrom_groupby.agg(start=('pos38', min), end=('pos38', max))
    df_chrom_region = df_chrom_groupby.agg(start=('pos38', min), end=('pos38', max),
                         gwas_category_lst=('gwas_category_lst', lambda x: set(chain.from_iterable(x))),
                         eqtl_gene_symbol_lst=('eqtl_gene_symbol_lst', lambda x: set(chain.from_iterable(x))),
                         etissue_category_term_lst=('etissue_category_term_lst', lambda x: set(chain.from_iterable(x))),
                         )
    df_pleio_chrom_grouped = df_chrom_groupby.agg({'cytoband': lambda x: sorted(x.unique()),
                                                 'pos38': lambda x: str(min(x)) + "-" + str(max(x)),
                                                 'rsid': lambda x: sorted(x.unique()),
                                                 'gwas_category_lst': lambda x: sorted(set(x.sum()))
                                                   })
    df_pleio_chrom_grouped['gwas_category_count'] = df_pleio_chrom_grouped['gwas_category_lst'].apply(lambda x: len(x))
    region_df = pandas.concat([region_df, df_chrom_region], axis=0, verify_integrity=False)
    import pdb; pdb.set_trace()

region_df.reset_index(inplace=True)
region_df = region_df.sort_values(['gwas_category_count', 'chrom', 'pos38'], ascending=[False, True, True])

# df = pandas.read_csv(count_per_rsid_gwas_egene_etissue_ods, sep="\t")

region_df[['start', 'end']] = region_df['pos38'].str.split('-', expand=True)
region_df['start'] = region_df['start'].astype(int)
region_df['end'] = region_df['end'].astype(int)


region_df['pleio_rsid'] = None
region_df['pleio_pos38'] = None
# pleio_df['name'] = None
for i, row in region_df.iterrows():
    start = row['start']
    end = row['end']
    variant_pleio = (df.loc[(df['pos38'] >= start) & (df['pos38'] <= end)]).head(1)
    # region_name = \
    # (variant_lead['cytoband'] + "_" + variant_lead['pos38'].astype(str) + "_" + variant_lead['rsid']).values[0]
    region_df.loc[i, 'pleio_rsid'] = variant_pleio['rsid'].values[0]
    region_df.loc[i, 'pleio_pos38'] = variant_pleio['pos38'].values[0]

region_df = region_df[['chrom', 'start', 'end', 'cytoband', 'pleio_rsid', 'pleio_pos38', 'gwas_category_count', 'gwas_category_lst']]
region_df.to_csv(pleio_tsv_path, sep="\t", index=False, header=True)


# ########################
# #%% bed file
# region_lst = []
# pos_prev = 0
# region_pleio_prev = False
# start = math.nan
# end = math.nan
# gwas_category_count = 0
# gwas_category_lst = math.nan
# category_lst = []
#
# for i, row in df.iterrows():
#     # beginning of region, set start, start category list, store category
#     if row['region_pleio'] and not region_pleio_prev:
#         cytoband = row['cytoband']
#         start = row['pos38']
#         gwas_category_count = row['gwas_category_count']
#         gwas_category_lst = row['gwas_category_lst']
#         category_lst = row['gwas_category_lst'].split(",")
#     # end of region, set end, store category
#     elif not row['region_pleio'] and region_pleio_prev:
#         end = pos_prev
#     # middle of region, store categories
#     if row['region_pleio'] and row['gwas_category_count'] > gwas_category_count:
#         gwas_category_count = row['gwas_category_count']
#         gwas_category_lst = row['gwas_category_lst']
#         category_lst = category_lst + row['gwas_category_lst'].split(",")
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
# df = pandas.DataFrame(region_lst, columns=['chrom', 'cytoband', 'start', 'end', 'gwas_category_count', 'gwas_category_lst'])
# df.sort_values(['gwas_category_count', 'chrom', 'start'], ascending=[False, True, True], inplace=True)
# df.to_csv(pleio_tsv_path, sep="\t", index=False, header=True)
#
#%##############
#%% GWAS region pleiotropy for MS
# import pdb; pdb.set_trace()
regions_pleio_ms_df = region_df.copy()

regions_pleio_ms_df = regions_pleio_ms_df.sort_values(by=['gwas_category_count', 'chrom', 'start'], ascending=[False, True, True])
# regions_pleio_ms_df = regions_pleio_ms_df.drop_duplicates('cytoband', keep='first')
regions_pleio_ms_df = regions_pleio_ms_df.loc[regions_pleio_ms_df['gwas_category_count'] >= manuscript_pleio_cutoff]

# format output
# import pdb; pdb.set_trace()
regions_pleio_ms_df.drop(['gwas_category_count'], inplace=True, axis=1)
regions_pleio_ms_df['gwas_category_lst'] = regions_pleio_ms_df['gwas_category_lst'].apply(lambda x: '; '.join(x))
regions_pleio_ms_df['start'] = regions_pleio_ms_df['start'].apply(lambda x: '{0:,}'.format(x))
regions_pleio_ms_df['end'] = regions_pleio_ms_df['end'].apply(lambda x: '{0:,}'.format(x))
regions_pleio_ms_df['cytoband'] = regions_pleio_ms_df['cytoband'].apply(lambda x: '; '.join(x) )
regions_pleio_ms_df = regions_pleio_ms_df[['chrom', 'pleio_rsid', 'pleio_pos38', 'cytoband', 'start', 'end', 'gwas_category_lst']]
regions_pleio_ms_df.rename({'gwas_category_lst': 'GWAS categoryes'}, axis=1, inplace=True)

tsv_path = os.path.join(outdir_path, "region_window_ms.tsv")
regions_pleio_ms_df.to_csv(tsv_path, sep="\t", index=False)

#%% bed
region_df['start'] = region_df['start'] - 1
region_df['chrom'] = 'chr' + region_df['chrom'].astype('str')
region_df.drop(['cytoband'], axis=1, inplace=True)
pleio_bed_path = os.path.join(outdir_path, "region_window_{}.bed".format(region_bin))
region_df.to_csv(pleio_bed_path, sep="\t", index=False, header=False)
#
# #%########################################### bed files
# for count_pleio in range(1, 6):
#     region_pleio_i_bed_path = os.path.join(outdir_path, "region_window_{}_pleio_{}.bed".format(region_bin, count_pleio))
#     region_pleio_i_df = pleio_df.loc[pleio_df['gwas_category_count'] == count_pleio,]
#     region_pleio_i_df.to_csv(region_pleio_i_bed_path, sep="\t", index=False, header=False)

#############################################
# Distribution of region by length
# % region <= 100kb, % >100kb and <=200kb, ...
region_length_ser = (region_df['end'] - region_df['start'])
data_ser = region_length_ser / wlength

ax = plt.hist(data_ser, bins=range(50), density=True)
plt.xlim(0, 10)
plt.tight_layout()
plt.savefig(png_path)
plt.close()

#%##############
# histogram
plt.rcParams["figure.figsize"] = (8, 6)
region_length_ser = (region_df['end'] - region_df['start'])

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
