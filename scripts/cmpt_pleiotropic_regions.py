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
    pleio_ods_path = sys.argv[4]
    png_path = sys.argv[5]
    if len(sys.argv) > 6:
        print("""Two many arguments!
        {}""".format(help_cmd_str))
        sys.exit(1)
except IndexError:
    print("""Argument missing!
    {}""".format(help_cmd_str))
    sys.exit(1)
outdir_path = os.path.dirname(pleio_ods_path)
pathlib.Path(outdir_path).mkdir(parents=True, exist_ok=True)

#%%
df = pandas.read_excel(count_per_rsid_gwas_egene_etissue_ods, engine='odf')
df['gwas_category_lst'] = df['gwas_category_lst'].str.split(';')
df['egene_lst']= df['egene_lst'].str.split(';')
df['eqtl_gene_symbol_lst'] = df['eqtl_gene_symbol_lst'].str.split(';')
df['etissue_category_term_lst'] = df['etissue_category_term_lst'].str.split(';')

region_df = pandas.DataFrame()
df.sort_values(['chrom', 'pos38'], inplace=True)
df['gwas_category_count'] = 1

################################################################################
for chrom in sorted(df['chrom'].unique()):  # problem split by chromosome
    # for chrom in [22]:  # problem split by chromosome
    chrom_df = df.loc[df['chrom'] == chrom].copy()

    # compute distance between consecutive eQTLs
    chrom_df['dist'] = [0] + [chrom_df.iloc[i]['pos38'] - chrom_df.iloc[i - 1]['pos38'] for i in range(1, chrom_df.shape[0])]
    region_pleio = 0  # init region pleio label
    chrom_df['region_pleio'] = region_pleio
    for i, row in chrom_df.iloc[1:].iterrows():  # Set or increase region label according to distance and wlength
        if row['dist'] <= wlength:
            chrom_df.loc[i, 'region_pleio'] = region_pleio
        else:
            region_pleio = region_pleio + 1
            chrom_df.loc[i, 'region_pleio'] = region_pleio
    chrom_df.drop('dist', inplace=True, axis=1)
    # groupby regions
    df_chrom_groupby = chrom_df.groupby(['chrom', 'region_pleio'])
    # region min and max
    chrom_region_df = df_chrom_groupby.agg(start=('pos38', min), end=('pos38', max),
                                           cytoband=('cytoband', lambda x: sorted(set(x))),
                                           rsid=('rsid', lambda x: sorted(x)),
                                           gwas_category=('gwas_category_lst', lambda x: sorted(set(chain.from_iterable(x)))),
                                           eqtl_gene_symbol=('eqtl_gene_symbol_lst', lambda x: set(chain.from_iterable([str(i) for i in x]))),
                                           etissue_category_term=('etissue_category_term_lst', lambda x: sorted(set(chain.from_iterable(x)))),
                                           )
    # add gene marker
    pubmed_df = chrom_df[
        ['chrom', 'eqtl_gene_marker_symbol', 'eqtl_gene_marker_id', 'pubmed_count', 'region_pleio']].sort_values(
        'pubmed_count', ascending=False).drop_duplicates(['chrom', 'region_pleio'], keep='first')
    pubmed_df.set_index(['chrom', 'region_pleio'], verify_integrity=True, inplace=True)
    # import pdb; pdb.set_trace()
    # pubmed_df.drop('pubmed_count', axis=1, inplace=True)
    # Merge region-level and variant-level information
    pubmed_df.rename({'eqtl_gene_marker_id': 'region_marker_id', 'eqtl_gene_marker_symbol': 'region_marker_symbol', 'pubmed_count': 'region_pubmed_count'},
                     axis=1, inplace=True)
    chrom_region_df = chrom_region_df.merge(pubmed_df, left_index=True, right_on=['chrom', 'region_pleio'])

    region_df = pandas.concat([region_df, chrom_region_df], axis=0, verify_integrity=False)

# count
region_df['gwas_category_count'] = region_df['gwas_category'].apply(len)
region_df['eqtl_gene_symbol_count'] = region_df['eqtl_gene_symbol'].apply(len)
region_df['etissue_category_term_count'] = region_df['etissue_category_term'].apply(len)

# list joined with ;
region_df['cytoband'] = region_df['cytoband'].apply(lambda x: ";".join(sorted([i for i in x if not i is None])))
region_df['gwas_category'] = region_df['gwas_category'].apply(lambda x: ";".join(sorted([i for i in x if not i is None])))
region_df['eqtl_gene_symbol'] = region_df['eqtl_gene_symbol'].apply(lambda x: ";".join(sorted([i for i in x if not i is None])))
region_df['etissue_category_term'] = region_df['etissue_category_term'].apply(lambda x: ";".join(sorted([i for i in x if not i is None])))
region_df['rsid'] = region_df['rsid'].apply(lambda x: ";".join(sorted([str(i) for i in x if not i is None])))

# write to ods
region_df.reset_index(inplace=True)
region_df.drop('region_pleio', axis=1, inplace=True)
columns = ['chrom', 'start', 'end', 'cytoband', 'region_marker_symbol',
           'gwas_category_count', 'gwas_category', 'eqtl_gene_symbol_count',
           'eqtl_gene_symbol', 'etissue_category_term_count', 'etissue_category_term', 'region_marker_id', 'rsid', 'region_pubmed_count']
region_df = region_df[columns]
region_df = region_df.sort_values(['gwas_category_count', 'chrom', 'start'], ascending=[False, True, True])

with pandas.ExcelWriter(pleio_ods_path, engine="odf") as writer:
    region_df.to_excel(writer, index=False)

#%% bed
bed_df = region_df.copy()
bed_df['start'] = bed_df['start'] - 1
bed_df['chrom'] = 'chr' + bed_df['chrom'].astype('str')
pleio_bed_path = os.path.join(outdir_path, "region_window_{}.bed".format(region_bin))
bed_df.to_csv(pleio_bed_path, sep="\t", index=False, header=False)

#%##############
#%% GWAS region pleiotropy for MS
# import pdb; pdb.set_trace()
regions_pleio_ms_df = region_df.copy()
regions_pleio_ms_df = regions_pleio_ms_df.sort_values(by=['gwas_category_count', 'chrom', 'start'], ascending=[False, True, True])
regions_pleio_ms_df = regions_pleio_ms_df.loc[regions_pleio_ms_df['gwas_category_count'] >= manuscript_pleio_cutoff]

# format output
regions_pleio_ms_df.drop(['gwas_category_count'], inplace=True, axis=1)
regions_pleio_ms_df['start'] = regions_pleio_ms_df['start'].apply(lambda x: '{0:,}'.format(x))
regions_pleio_ms_df['end'] = regions_pleio_ms_df['end'].apply(lambda x: '{0:,}'.format(x))
regions_pleio_ms_df['gwas_category'] = regions_pleio_ms_df['gwas_category'].str.replace(';', '; ')
regions_pleio_ms_df = regions_pleio_ms_df[['chrom', 'start', 'end', 'cytoband', 'region_marker_symbol', 'gwas_category']]
tsv_path = os.path.join(outdir_path, "region_window_ms.tsv")
regions_pleio_ms_df.to_csv(tsv_path, sep="\t", index=False)

#%% bed
region_df['start'] = region_df['start'] - 1
region_df['chrom'] = 'chr' + region_df['chrom'].astype('str')
region_df.drop(['cytoband'], axis=1, inplace=True)
pleio_bed_path = os.path.join(outdir_path, "region_window_{}.bed".format(region_bin))
region_df.to_csv(pleio_bed_path, sep="\t", index=False, header=False)

#%##############################################################################
#
# histogram
#
################################################################################

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

# ax = plt.hist(data_ser, bins=range(50), density=True, cumulative=True)
ax = plt.hist(data_ser, bins=[i/100 for i in range(0,1001)], density=True, cumulative=True)
plt.xlim(0, 1)

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
