"""At each variant pleiotropy level, this code computes and plots the Fisher
odds ratio and p-value of CRM vs. Non-CRM variants"""

import numpy
import os
import pandas
import pathlib
import seaborn
import shlex
import subprocess
import sys

from gwas2eqtl_pleiotropy.Logger import Logger
from gwas2eqtl_pleiotropy.constants import public_data_dir, annotator_config_dic, label_fontsize, tick_fontsize, dpi, seaborn_theme_dic
from matplotlib import pyplot as plt
from scipy.stats import fisher_exact
from statannotations.Annotator import Annotator

seaborn.set_theme(**seaborn_theme_dic)

#%%
help_cmd_str = "todo"
try:
    count_per_rsid_gwas_tsv_path = sys.argv[1]
    variant_pleio_1_flank_10_hg38_bed = sys.argv[2]
    max_gwas_class_count = int(sys.argv[3])
    remap_crm_png = sys.argv[4]
    if len(sys.argv) > 5:
        print("""Two many arguments!
        {}""".format(help_cmd_str))
        sys.exit(1)
except IndexError:
    print("""Argument missing!
    {}""".format(help_cmd_str))
    sys.exit(1)

outdir_path = os.path.join(os.path.dirname(remap_crm_png))
pathlib.Path(outdir_path).mkdir(parents=True, exist_ok=True)

#%% input dir cmpt_count_per_rsid
indir_path = os.path.dirname(count_per_rsid_gwas_tsv_path)

remap_crm_path = os.path.join(public_data_dir, "remap.univ-amu.fr/storage/remap2022/hg38/MACS2/remap2022_crm_macs2_hg38_v1_0.bed.gz")

out_df_columns = ['pleio_count', 'pleio_n_crm_count', 'pleio_1_crm_count', 'pleio_n_nocrm_count', 'pleio_1_nocrm_count', 'oddsr', 'p']
out_df = pandas.DataFrame(columns = out_df_columns)

#%% bedtools intersect
flank = 10
for count_pleio in range(1, max_gwas_class_count+1):
    bed_path = os.path.join(indir_path, "variant_pleio_{}_flank_{}_hg38.bed".format(count_pleio, flank))
    intersect_bed_path = os.path.join(outdir_path, "remap_crm_pleio_{}.bed".format(count_pleio))
    cmd_stf = "bedtools intersect -sorted -a {bed_path} -b {remap_crm_path} -loj -wb"
    cmd = cmd_stf.format(**{'bed_path': bed_path, 'remap_crm_path': remap_crm_path, 'output_bed': intersect_bed_path})
    Logger.info(cmd)
    with open(intersect_bed_path, 'w') as fout:
        result = subprocess.run(shlex.split(cmd), stdout=fout)

    crm_pleio_df = pandas.read_csv(intersect_bed_path, sep="\t", header=None, usecols=[0, 1, 2, 3, 4, 5, 6], names=['chrom', 'start', 'end', 'rsid', 'gwas_class_count', 'gwas_class_lst', 'crm'])

    pleio_n_nocrm_count = (crm_pleio_df['crm'] == '.').sum()
    pleio_n_crm_count = (crm_pleio_df['crm'] != '.').sum()

    if count_pleio == 1:  # fisher test

        pleio_1_nocrm_count = pleio_n_nocrm_count
        pleio_1_crm_count = pleio_n_crm_count

    a = pleio_n_crm_count
    b = pleio_1_crm_count
    c = pleio_n_nocrm_count
    d = pleio_1_nocrm_count
    table = numpy.array([[a, b], [c, d]])
    oddsr, p = fisher_exact(table, alternative='greater')
    out_row_lst = [count_pleio, pleio_n_crm_count, pleio_1_crm_count,
                   pleio_n_nocrm_count, pleio_1_nocrm_count, oddsr, p]
    out_df = pandas.concat([out_df, pandas.DataFrame(
        dict(zip(out_df_columns, out_row_lst)), index=[count_pleio])], axis=0)

####################################################################

#%% set signif symbols
out_df['signif'] = "ns"
out_df.loc[out_df['p'] <= 5.00e-02, 'signif'] = '*'
out_df.loc[out_df['p'] <= 1.00e-02, 'signif'] = '**'
out_df.loc[out_df['p'] <= 1.00e-03, 'signif'] = '***'
out_df.loc[out_df['p'] <= 1.00e-04, 'signif'] = '****'
out_df.rename({'pleio_count': 'gwas_class_count'}, axis=1, inplace=1)

#%%
# out_df = in_df.loc[in_df['consequence'] == consequence, ['gwas_class_count', 'oddsr', 'p', 'signif']]

#%%
out_df['gwas_class_count'] = [str(i) for i in out_df['gwas_class_count']]
order = out_df['gwas_class_count'].tolist()
xticklabels = order.copy()
xticklabels[-1] = '≥{}'.format(order[-1])
title = "CRM annotation"
xlabel = "GWAS category count"
ylabel = "Odds ratio"
y = "oddsr"
x = "gwas_class_count"

#%%
ax = seaborn.barplot(x=x, y=y, data=out_df, order=order, palette="rocket_r")

#%%
pairs = [('1', x) for x in out_df['gwas_class_count'] if x != "1"]
formatted_pvalues = out_df['signif'].tolist()[1:]

annotator = Annotator(ax, pairs, data=out_df, x=x, y=y, order=order, size=label_fontsize)
annotator.set_custom_annotations(formatted_pvalues)
annotator.configure(**annotator_config_dic)
annotator.annotate()

ax.set_xticklabels(xticklabels)
plt.grid(axis="y")
plt.title(title, fontsize=label_fontsize)
plt.xlabel(xlabel, fontsize=label_fontsize)
plt.xticks(fontsize=tick_fontsize, rotation=0)
plt.ylabel(ylabel, fontsize=label_fontsize)
plt.yticks(fontsize=tick_fontsize)

plt.tight_layout()
plt.savefig(remap_crm_png, dpi=dpi)
plt.close()
