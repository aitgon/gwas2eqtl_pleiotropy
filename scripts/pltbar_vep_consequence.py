import os
import pathlib
import sys
import pandas
import seaborn

from gwas2eqtl_pleiotropy.constants import dpi, seaborn_theme_dic, label_fontsize, tick_fontsize, palette
from matplotlib import pyplot as plt

seaborn.set_theme(**seaborn_theme_dic)


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
# import pdb; pdb.set_trace()
consequence_signif_lst = in_df.loc[in_df['pfdr5perc'] < 0.05, 'consequence'].drop_duplicates().tolist()
in_df = in_df.loc[in_df['consequence'].isin(consequence_signif_lst)]
nonzero_count_x_lst = in_df.loc[in_df['a_pleio_x_with_consequence'] > 0, 'consequence'].drop_duplicates().tolist()
in_df = in_df.loc[in_df['consequence'].isin(nonzero_count_x_lst)]
nonzero_count_1_lst = in_df.loc[in_df['b_pleio_1_with_consequence'] > 0, 'consequence'].drop_duplicates().tolist()
in_df = in_df.loc[in_df['consequence'].isin(nonzero_count_1_lst)]

#%% set signif symbols
in_df['signif'] = "ns"
in_df.loc[in_df['pfdr5perc'] <= 5.00e-02, 'signif'] = '*'
in_df.loc[in_df['pfdr5perc'] <= 1.00e-02, 'signif'] = '**'
in_df.loc[in_df['pfdr5perc'] <= 1.00e-03, 'signif'] = '***'
in_df.loc[in_df['pfdr5perc'] <= 1.00e-04, 'signif'] = '****'
# in_df.loc[in_df['gwas_category_count'] == max_gwas_category_count, 'gwas_category_count'] = '≥{}'.format(max_gwas_category_count)

################################################################################
# Draw a nested barplot by species and sex
in_df = in_df.sort_values(['consequence', 'gwas_category_count'])
in_df.loc[in_df['gwas_category_count'] == max(in_df['gwas_category_count']), 'gwas_category_count'] = '≥' + max(in_df['gwas_category_count'])

order = in_df['consequence'].unique()
hue_order = in_df['gwas_category_count'].unique()
g = seaborn.barplot(data=in_df, x="consequence", y="oddsr", hue="gwas_category_count", palette=palette, orient='v', order=order, hue_order=hue_order)
for c in g.containers:
    gwas_category_count = c.get_label()
    labels = in_df.loc[in_df['gwas_category_count'] == gwas_category_count, 'signif'].tolist()
    g.bar_label(c, labels=labels, label_type='edge', padding=5, fontsize=label_fontsize)

# legend_labels = hue_order
# legend_labels[-1] = '≥' + hue_order[-1]
# g.legend(labels=legend_labels, title='Category cnt.', loc='best', title_fontsize='xx-large', fontsize='18')
g.legend(title='Category cnt.', loc='best', title_fontsize='xx-large', fontsize='18')
plt.title("Variant Effect Predictor", fontsize=label_fontsize)
plt.ylabel("Odds vs cat. count 1", fontsize=tick_fontsize)
plt.yticks(fontsize=12, rotation=0)
plt.xlabel(None)
plt.ylim([0, 11])
plt.xticks(fontsize=label_fontsize)

plt.tight_layout()
plt.savefig(consequence_png_path, dpi=dpi)
plt.close()
