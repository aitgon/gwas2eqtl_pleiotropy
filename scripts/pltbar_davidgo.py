# import ssl
# ssl._create_default_https_context = ssl._create_unverified_context
import os
import pathlib

import matplotlib
import numpy
import pandas
import sys

import logging
import traceback as tb

import pandas
import seaborn
import suds.metrics as metrics
# from tests import *
from suds import *
from suds.client import Client
from datetime import datetime
from matplotlib import pyplot as plt

plt.rcParams["figure.figsize"] = (20, 4)


#%%
from gwas2eqtl_pleiotropy.Logger import Logger
from gwas2eqtl_pleiotropy.constants import label_fontsize, tick_fontsize, dpi

help_cmd_str = "todo"
try:
    count_per_rsid_gwas_egene_etissue_ods_path = sys.argv[1]
    davidgo_tsv_path = sys.argv[2]
    davidgo_png_path = sys.argv[3]
    if len(sys.argv) > 4:
        print("""Two many arguments!
        {}""".format(help_cmd_str))
        sys.exit(1)
except IndexError:
    print("""Argument missing!
    {}""".format(help_cmd_str))
    sys.exit(1)

outdir_path = os.path.dirname(davidgo_png_path)
pathlib.Path(outdir_path).mkdir(parents=True, exist_ok=True)

fin_df = pandas.read_excel(count_per_rsid_gwas_egene_etissue_ods_path, engine='odf')

#%%
max_gwas_category_count = fin_df['gwas_category_count'].max()
for pleio_i in range(2, max_gwas_category_count + 1):
    Logger.info(pleio_i)
    davidgo_pleio_tsv_path = davidgo_tsv_path + "_{}.tsv".format(pleio_i)
    fin_df = pandas.read_csv(davidgo_pleio_tsv_path, sep="\t", header=0)
    fin_df = fin_df.loc[fin_df['Bonferroni'] <= 0.05]

    if fin_df.shape[0] == 0:
        continue
    fin_df['Term'] = fin_df['Term'].str.split('~', expand=True)[1]
    fin_df = fin_df[['Term', 'Bonferroni']]
    fin_df['Bonferroni'] = -numpy.log10(fin_df['Bonferroni'])

    ax = seaborn.barplot(y="Term", x="Bonferroni", data=fin_df, orient='h')

    # ax.set(xticklabels=[])
    plt.grid(axis='y')
    plt.legend(fontsize=16)  # using a size in points
    plt.title("GWAS cat. count {}".format(pleio_i), fontsize=label_fontsize)
    plt.xlabel("Neg. log10 p-val", fontsize=label_fontsize)
    # plt.xticks(fontsize=tick_fontsize, rotation=0)
    plt.ylabel("GO term", fontsize=label_fontsize)
    plt.yticks(fontsize=tick_fontsize)

    plt.tight_layout()
    png_path = davidgo_png_path + "_{}.png".format(pleio_i)
    Logger.info(png_path)
    plt.savefig(png_path, dpi=dpi)
    plt.close()
