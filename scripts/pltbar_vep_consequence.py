import os
import pathlib
import sys
import pandas
import seaborn
from statannotations.Annotator import Annotator

from gwas2eqtl_pleiotropy.constants import dpi, seaborn_theme_dic, label_fontsize, tick_fontsize, palette, \
    annotator_config_dic
from matplotlib import pyplot as plt

seaborn.set_theme(**seaborn_theme_dic)


#%%
help_cmd_str = "todo"
try:
    consequence_tsv_path = sys.argv[1]
    outdir_path = sys.argv[2]
    if len(sys.argv) > 3:
        print("""Two many arguments!
        {}""".format(help_cmd_str))
        sys.exit(1)
except IndexError:
    print("""Argument missing!
    {}""".format(help_cmd_str))
    sys.exit(1)

#%%
# outdir_path = os.path.dirname(consequence_png_path)
pathlib.Path(outdir_path).mkdir(exist_ok=True, parents=True)

#%%
df1 = pandas.read_csv(consequence_tsv_path, sep="\t", header=0)
df1['gwas_category_count'] = df1['gwas_category_count'].astype(int).astype(str)
consequence_signif_lst = df1.loc[df1['pfdr5perc'] < 0.05, 'consequence'].drop_duplicates().tolist()

for consequence in consequence_signif_lst:

    print(consequence)
    df2 = df1.loc[df1['consequence'] == consequence]
    nonzero_count_x_lst = df2.loc[df2['a_pleio_x_with_consequence'] > 0, 'consequence'].drop_duplicates().tolist()
    df2 = df2.loc[df2['consequence'].isin(nonzero_count_x_lst)]
    nonzero_count_1_lst = df2.loc[df2['b_pleio_1_with_consequence'] > 0, 'consequence'].drop_duplicates().tolist()
    df2 = df2.loc[df2['consequence'].isin(nonzero_count_1_lst)]

    #%% set signif symbols
    df2['signif'] = "ns"
    df2.loc[df2['pfdr5perc'] <= 5.00e-02, 'signif'] = '*'
    df2.loc[df2['pfdr5perc'] <= 1.00e-02, 'signif'] = '**'
    df2.loc[df2['pfdr5perc'] <= 1.00e-03, 'signif'] = '***'
    df2.loc[df2['pfdr5perc'] <= 1.00e-04, 'signif'] = '****'

    df2 = pandas.concat([df2, pandas.DataFrame(
        {'consequence': consequence, 'gwas_category_count': "1", 'oddsr': 1, 'signif': 'ns'}, index=[0])])
    df2.sort_values('gwas_category_count', inplace=True)

    title = consequence
    xlabel = "Trait category count"
    ylabel = "Odds ratio"
    y = "oddsr"
    x = "gwas_category_count"
    order = sorted(df2['gwas_category_count'].tolist())
    pairs = [('1', x) for x in df2['gwas_category_count'] if x != "1"]
    formatted_pvalues = df2['signif'].tolist()[1:]
    xticklabels = order.copy()


    # %% barplot
    ax = seaborn.barplot(x=x, y=y, data=df2, order=order, palette=palette)

    annotator = Annotator(ax, pairs, data=df2, x=x, y=y, order=order, size=label_fontsize)
    annotator.set_custom_annotations(formatted_pvalues)
    annotator.configure(**annotator_config_dic)
    annotator.annotate()

    ax.set_xticklabels(xticklabels)
    plt.title(title, fontsize=label_fontsize)
    plt.xlabel(xlabel, fontsize=label_fontsize)
    # plt.xticks(fontsize=tick_fontsize, rotation=0)
    xticks_labels = [str(x) for x in (plt.xticks()[0] + 1)]
    xticks_labels[-1] = '≥' + str(xticks_labels[-1])
    plt.xticks(ticks=(plt.xticks()[0]), labels=xticks_labels, fontsize=tick_fontsize, rotation=0)
    plt.ylabel(ylabel, fontsize=label_fontsize)
    plt.yticks(fontsize=tick_fontsize)
    plt.ylim([0, df2['oddsr'].max()*1.5])

    plt.tight_layout()
    png_path =os.path.join(outdir_path, "{}.png".format(consequence))
    plt.savefig(png_path, dpi=dpi)
    plt.close()
