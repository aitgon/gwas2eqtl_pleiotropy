"""This code takes a list of bed files with the regions of eQTLs with different pleiotropies and a region bed file, for instance CRMs from remap
It creates a barplot with the odds ration of the enrichment"""

import argparse

from gwas2eqtl_pleiotropy_logger import Logger
from constants import annotator_config_dic, label_fontsize, tick_fontsize, dpi, seaborn_theme_dic, \
    palette
from matplotlib import pyplot as plt
from scipy.stats import fisher_exact
from statannotations.Annotator import Annotator

import numpy
import os
import pandas
import pathlib
import seaborn
import shlex
import subprocess

seaborn.set_theme(**seaborn_theme_dic)

help_cmd_str = "todo."

def parse_arguments():
    parser = argparse.ArgumentParser(description="Process file paths for eQTL and CRM analysis.")
    parser.add_argument(
        "--title",
        help="Title of the barplot."
    )
    parser.add_argument(
        "--eqtl_bed",
        nargs=3,
        help="List of three BED files (e.g., eqtl_pleio_1_flank_10_hg38_bed, eqtl_pleio_2_flank_10_hg38_bed, eqtl_pleio_3_flank_10_hg38_bed)"
    )
    parser.add_argument(
        "--annot_bed",
        help="Path to the remap CRM file."
    )
    parser.add_argument(
        "--count_tsv",
        help="Path to the remap count TSV file."
    )
    parser.add_argument(
        "--barplot_png",
        help="Path to save the barplot remap CRM PNG file."
    )

    args = parser.parse_args()

    return args

#%%

args = parse_arguments()
title = args.title
#eqtl_pleio_1_flank_10_hg38_bed, eqtl_pleio_2_flank_10_hg38_bed, eqtl_pleio_3_flank_10_hg38_bed = args.bed_files
annot_bed_path = args.annot_bed
count_tsv_path = args.count_tsv
barplot_png_path = args.barplot_png

#%%

outdir_path = os.path.join(os.path.dirname(barplot_png_path))
pathlib.Path(outdir_path).mkdir(parents=True, exist_ok=True)

# %% input dir cmpt_count_per_rsid
#indir_path = os.path.dirname(eqtl_pleio_1_flank_10_hg38_bed)

out_df_columns = ['pleio_count', 'pleio_n_crm_count', 'pleio_1_crm_count', 'pleio_n_nocrm_count', 'pleio_1_nocrm_count',
                  'oddsr', 'p']
out_df = pandas.DataFrame(columns=out_df_columns)

# %% bedtools intersect
# flank = 10
# for count_pleio in range(1, 99):
#for count_pleio in range(1, pleio_high_cutoff + 1):
for count_pleio, eqtl_bed_path in enumerate(args.eqtl_bed, start=1):
    #eqtl_bed_path = os.path.join(indir_path, "eqtl_pleio_{}_flank_{}_hg38.bed".format(count_pleio, flank))
    if not os.path.isfile(eqtl_bed_path):
        break
    intersect_bed_path = os.path.join(outdir_path, "intesected_pleio_{}.bed".format(count_pleio))
    cmd_stf = "bedtools intersect -sorted -a {eqtl_bed_path} -b {annot_bed_path} -loj -wb -f 1"
    cmd = cmd_stf.format(
        **{'eqtl_bed_path': eqtl_bed_path, 'annot_bed_path': annot_bed_path, 'output_bed': intersect_bed_path})
    Logger.info(cmd)
    with open(intersect_bed_path, 'w') as fout:
        result = subprocess.run(shlex.split(cmd), stdout=fout)
    crm_pleio_df = pandas.read_csv(intersect_bed_path, sep="\t", header=None, usecols=[0, 1, 2, 3, 4, 5, 6, 7],
                                   names=['chrom', 'eqtl_start', 'eqtl_end', 'eqtl_rsid', 'gwas_category_count',
                                          'crm_chrom', 'crm_start', 'crm_end'])

    pleio_n_nocrm_count = (crm_pleio_df['crm_chrom'] == '.').sum()
    pleio_n_crm_count = (crm_pleio_df['crm_chrom'] != '.').sum()

    if count_pleio == 1:  # fisher test reference

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

    # Create the DataFrame to be added
    new_df = pandas.DataFrame(dict(zip(out_df_columns, out_row_lst)), index=[count_pleio])

    # Filter out columns that are empty or all-NA before concatenating
    new_df = new_df.dropna(axis=1, how='all')

    out_df = pandas.concat([out_df, new_df], axis=0)

out_df.to_csv(count_tsv_path, sep="\t", index=False)

####################################################################

# %% set signif symbols
out_df['signif'] = "ns"
out_df.loc[out_df['p'] <= 5.00e-02, 'signif'] = '*'
out_df.loc[out_df['p'] <= 1.00e-02, 'signif'] = '**'
out_df.loc[out_df['p'] <= 1.00e-03, 'signif'] = '***'
out_df.loc[out_df['p'] <= 1.00e-04, 'signif'] = '****'
out_df.rename({'pleio_count': 'gwas_category_count'}, axis=1, inplace=1)

# %%
out_df['gwas_category_count'] = [str(i) for i in out_df['gwas_category_count']]
order = out_df['gwas_category_count'].tolist()
xticklabels = order.copy()
#title = "CRM annotation enrichm."
xlabel = "Trait category count"
ylabel = "Odds ratio"
y = "oddsr"
x = "gwas_category_count"

# %%
pairs = [('1', x) for x in out_df['gwas_category_count'] if x != "1"]
formatted_pvalues = out_df['signif'].tolist()[1:]

# import pdb; pdb.set_trace()

# %% barplot
# ax = seaborn.barplot(x=x, y=y, data=out_df, order=order, palette=palette)
ax = seaborn.barplot(x=x, y=y, data=out_df, order=order, palette=palette, hue=x, legend=False)
#import pdb; pdb.set_trace()
annotator = Annotator(ax, pairs, data=out_df, x=x, y=y, order=order, size=label_fontsize)
annotator.set_custom_annotations(formatted_pvalues)
annotator.configure(**annotator_config_dic)
annotator.annotate()

ax.set_xticklabels(xticklabels)
plt.title(title, fontsize=label_fontsize)
plt.xlabel(xlabel, fontsize=label_fontsize)
# plt.xticks(fontsize=tick_fontsize, rotation=0)
# xticks_labels = [str(x) for x in (plt.xticks()[0] + 1)]
xticks_labels = [str(x + 1) for x in plt.xticks()[0]]
xticks_labels[-1] = '≥' + str(xticks_labels[-1])
plt.xticks(ticks=(plt.xticks()[0]), labels=xticks_labels, fontsize=tick_fontsize, rotation=0)
plt.ylabel(ylabel, fontsize=label_fontsize)
plt.yticks(fontsize=tick_fontsize)

plt.tight_layout()
plt.savefig(barplot_png_path, dpi=dpi)
plt.close()
