import os
import pathlib
import sys
import pandas
import seaborn

from gwas2eqtl_pleiotropy.constants import tick_fontsize, label_fontsize, dpi, annotator_config_dic, seaborn_theme_dic
from matplotlib import pyplot as plt
from statannotations.Annotator import Annotator

seaborn.set_theme(**seaborn_theme_dic)

#%%
help_cmd_str = "todo"
try:
    consequence_tsv_path = sys.argv[1]
    upper_var_gwas_cat_count = sys.argv[2]
    consequence_png_path = sys.argv[3]
    if len(sys.argv) > 4:
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
in_df.loc[in_df['gwas_category_count'] >= upper_var_gwas_cat_count, "gwas_category_count"] = upper_var_gwas_cat_count

#%% set signif symbols
in_df['signif'] = "ns"
in_df.loc[in_df['p'] <= 5.00e-02, 'signif'] = '*'
in_df.loc[in_df['p'] <= 1.00e-02, 'signif'] = '**'
in_df.loc[in_df['p'] <= 1.00e-03, 'signif'] = '***'
in_df.loc[in_df['p'] <= 1.00e-04, 'signif'] = '****'

#%%
# consequence = 'missense_variant'

#%%
for consequence in in_df['consequence'].unique():

    #%%
    plt_df = in_df.loc[in_df['consequence'] == consequence, ['gwas_category_count', 'oddsr', 'p', 'signif']]

    #%%
    order = sorted(plt_df['gwas_category_count'].unique())
    xticklabels = order.copy()
    xticklabels[-1] = '≥{}'.format(order[-1])
    title = consequence
    xlabel = "GWAS category count"
    ylabel = "Odds ratio"
    # ylim = [0, 9]
    y = "oddsr"
    x = "gwas_category_count"

    #%%
    ax = seaborn.barplot(x=x, y=y, data=plt_df, order=order, palette="rocket_r")

    #%%
    formatted_pvalues = plt_df['signif'].tolist()

    #%%
    pairs = [(x, x) for x in plt_df['gwas_category_count']]
    annotator = Annotator(ax, pairs, data=plt_df, x=x, y=y, order=order, size=label_fontsize)
    annotator.set_custom_annotations(formatted_pvalues)
    annotator.configure(**annotator_config_dic)
    annotator.annotate()

    plt.title(title, fontsize=label_fontsize)
    plt.xlabel(xlabel, fontsize=label_fontsize)
    plt.xticks(fontsize=tick_fontsize, rotation=0)
    plt.ylabel(ylabel, fontsize=label_fontsize)
    # plt.ylim(ylim)
    plt.yticks(fontsize=tick_fontsize)
    ax.set_xticklabels(xticklabels)

    plt.tight_layout()
    plt_path = os.path.join(outdir_path, consequence + ".png")
    plt.savefig(plt_path, dpi=dpi)
    plt.close()
