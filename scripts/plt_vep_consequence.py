import os
import pathlib
import sys
import pandas
import seaborn

from eqtl2gwas_pleiotropy.constants import tick_fontsize, label_fontsize, dpi
from matplotlib import pyplot as plt
from statannotations.Annotator import Annotator

consequence_tsv_path = 'out/gwas413/genome/5e-08/1000000/cmpt_vep_consequence_fisher.py/cons_stat.tsv'
consequence_png_path = 'out/gwas413/genome/5e-08/1000000/plt_vep_consequence.py/bar.png'

#%%
help_cmd_str = "todo"
try:
    consequence_tsv_path = sys.argv[1]
    consequence_png_path = sys.argv[2]
    if len(sys.argv) > 3:
        print("""Two many arguments!
        {}""".format(help_cmd_str))
        sys.exit(1)
except IndexError:
    print("""Argument missing!
    {}""".format(help_cmd_str))
    sys.exit(1)

#%%
outdir_path = os.path.dirname(consequence_png_path)
pathlib.Path(outdir_path).mkdir(exist_ok=True, parents=True)

#%%
in_df = pandas.read_csv(consequence_tsv_path, sep="\t", header=0)
in_df['gwas_category_count'] = in_df['gwas_category_count'].astype(int).astype(str)

#%% set signif symbols
in_df['signif'] = "ns"
in_df.loc[in_df['p'] <= 5.00e-02, 'signif'] = '*'
in_df.loc[in_df['p'] <= 1.00e-02, 'signif'] = '**'
in_df.loc[in_df['p'] <= 1.00e-03, 'signif'] = '***'
in_df.loc[in_df['p'] <= 1.00e-04, 'signif'] = '****'

#%%
consequence = 'missense_variant'

#%%
for consequence in in_df['consequence'].unique():

    #%%
    plt_df = in_df.loc[in_df['consequence'] == consequence, ['gwas_category_count', 'oddsr', 'p', 'signif']]

    #%%
    order = sorted(plt_df['gwas_category_count'].unique())
    seaborn.set_theme(style="whitegrid")
    xticklabels = order.copy()
    title = consequence
    xlabel = "GWAS category count"
    ylabel = "Odds ratio"
    ylim = [0, 9]
    y = "oddsr"
    x = "gwas_category_count"

    #%%
    ax = seaborn.barplot(x=x, y=y, data=plt_df, order=order, palette="rocket_r", ci=None)

    #%%
    formatted_pvalues = plt_df['signif'].tolist()

    #%%
    pairs = [(x, x) for x in plt_df['gwas_category_count']]
    annotator = Annotator(ax, pairs, data=plt_df, x=x, y=y, order=order, size=label_fontsize)
    annotator.set_custom_annotations(formatted_pvalues)
    annotator.configure(**{'fontsize': 16})
    annotator.annotate()


    plt.title(title, fontsize=label_fontsize)
    plt.xlabel(xlabel, fontsize=label_fontsize)
    plt.xticks(fontsize=tick_fontsize, rotation=0)
    plt.ylabel(ylabel, fontsize=label_fontsize)
    plt.ylim(ylim)
    plt.yticks(fontsize=tick_fontsize)
    ax.set_xticklabels(xticklabels)

    plt.tight_layout()
    plt_path = os.path.join(outdir_path, consequence + ".png")
    plt.savefig(plt_path)
    plt.close()
