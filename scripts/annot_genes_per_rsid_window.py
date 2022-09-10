#%% Load refseq gene coordinates
import shlex
import subprocess

import pandas
import seaborn
from matplotlib import pyplot as plt
from statannot import add_stat_annotation

from gwas2eqtl_pleiotropy.Logger import Logger
from gwas2eqtl_pleiotropy.UCSC import UCSC

#%%
from gwas2eqtl_pleiotropy.constants import label_fontsize, tick_fontsize

refseq_df = UCSC().get_ref_gene_table(force=True)

#%%
refseq_df['name'] = refseq_df.index

#%%
refseq_df.index = range(refseq_df.shape[0])

#%%
refseq_bed_df = refseq_df.copy()

#%%
refseq_bed_df['txStart'] = refseq_bed_df['txStart'] - 1

#%%
refseq_bed_df = refseq_bed_df.sort_values(refseq_bed_df.columns.tolist())
refseq_bed_df.to_csv('refseq.bed', sep="\t", header=None, index=False)

#%%
count_per_rsid_gwas_tsv = "out/gwas413/genome/5e-08/1000000/cmpt_count_per_rsid.py/count_per_rsid_gwas.tsv"
rsid_df = pandas.read_csv(count_per_rsid_gwas_tsv, sep="\t", header=0)

#%%
rsid_1mb_bed_df = rsid_df.copy()

#%%
w = 10000
rsid_1mb_bed_df['chromStart'] = rsid_1mb_bed_df['pos'].astype(int) - int(w/2) - 1
rsid_1mb_bed_df.loc[rsid_1mb_bed_df['chromStart'] < 0, 'chromStart'] = 0
rsid_1mb_bed_df['chromEnd'] = rsid_1mb_bed_df['pos'] + int(w/2)
rsid_1mb_bed_df['chrom'] = 'chr' + rsid_1mb_bed_df['chrom'].astype('str')

#%%
rsid_1mb_bed_df = rsid_1mb_bed_df[['chrom', 'chromStart', 'chromEnd', 'rsid', 'gwas_category_count', 'gwas_category_lst', 'pos']]
rsid_1mb_bed_df = rsid_1mb_bed_df.sort_values(rsid_1mb_bed_df.columns.tolist())
rsid_1mb_bed_df.to_csv('rsid_1mb.bed', sep="\t", header=None, index=False)

#%%
intersect_bed_path = "intersected.bed"
cmd_stf = "bedtools intersect -sorted -a {rsid_1mb_bed} -b {refseq_bed} -wb"
cmd = cmd_stf.format(**{'refseq_bed': 'refseq.bed', 'rsid_1mb_bed': 'rsid_1mb.bed', 'intersect_bed_path': intersect_bed_path})
Logger.info(cmd)
with open(intersect_bed_path, 'w') as fout:
    result = subprocess.run(shlex.split(cmd), stdout=fout)

#%%
upper_var_gwas_cat_count = 5

#%%
df = pandas.read_csv(intersect_bed_path, sep="\t", header=None)
df = df[[3, 4, 12]].drop_duplicates()
df.columns = ['rsid', 'gwas_category_count', 'gene']
df.loc[df['gwas_category_count'] > upper_var_gwas_cat_count, 'gwas_category_count'] = upper_var_gwas_cat_count

df = df.groupby(['rsid', 'gwas_category_count']).gwas_cat_count().reset_index()

#%%
order = [*range(1, upper_var_gwas_cat_count+1)]
seaborn.set_theme(style="whitegrid")
xticklabels = order.copy()
xticklabels[-1] = '≥{}'.format(order[-1])
title = "ReMap TFs per var."
xlabel = "GWAS category count"
ylabel = "TF count"
y = "gene"

#%%
box_pairs = [(1, i) for i in range(2, upper_var_gwas_cat_count+1) ]
ax = seaborn.violinplot(x="gwas_category_count", y=y, data=df, order=order, palette="rocket_r")
test_results = add_stat_annotation(ax, data=df, x="gwas_category_count", y=y, order=order,
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
plt.savefig("gene_density.png")
plt.close()
