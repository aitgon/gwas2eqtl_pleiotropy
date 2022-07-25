import os
import sys

import pandas
import seaborn
from matplotlib import pyplot as plt

#%%
from eqtl2gwas_pleiotropy.constants import tick_fontsize, label_fontsize, dpi

upper_var_gwas_cat_count = 5
vep1_path = "out/gwas413/genome/5e-08/1000000/cmpt_vep_consequence.py/vep_pleio1.tsv"

#%%
help_cmd_str = "todo"
try:
    count_per_rsid_gwas_tsv_path = sys.argv[1]
    upper_var_gwas_cat_count = int(sys.argv[2])
    vep2_path = sys.argv[3]
    png_path = sys.argv[4]
    if len(sys.argv) > 5:
        print("""Two many arguments!
        {}""".format(help_cmd_str))
        sys.exit(1)
except IndexError:
    print("""Argument missing!
    {}""".format(help_cmd_str))
    sys.exit(1)


#%%
df = pandas.read_csv(count_per_rsid_gwas_tsv_path, sep="\t")
df.loc[df['gwas_category_count'] >= upper_var_gwas_cat_count, 'gwas_category_count'] = upper_var_gwas_cat_count
cat_df = pandas.DataFrame({'percentage': None, 'consequence': None, 'gwas_category_count': None}, index=[])
# import pdb; pdb.set_trace()
for pleio in range(1, upper_var_gwas_cat_count+1):
    vep_path = os.path.join(os.path.dirname(vep1_path), os.path.basename(vep1_path).replace("1", str(pleio)))
    # vep_path = vep0_path.format(pleio)
    vdf = pandas.read_csv(vep_path, sep="\t", comment='#', header=None)
    vdf = vdf.loc[vdf[0].drop_duplicates(keep='first').index]
    vdf = vdf[[0, 6]]
    vdf[6] = vdf[6].str.split(',', expand=True)[0]
    variant_count = vdf.shape[0]
    consequence_lst = vdf.groupby(6).count().index.tolist()
    percentage_lst = round(vdf.groupby(6).count() / variant_count * 100)[0].tolist()
    conseq_perc_df = pandas.DataFrame({'consequence': consequence_lst, 'percentage': percentage_lst}, index=[*range(len(consequence_lst))])
    conseq_perc_df['gwas_category_count'] = pleio
    cat_df = pandas.concat([cat_df, conseq_perc_df], axis=0)

cat_wide_df = pandas.pivot_table(cat_df, values='percentage', index='consequence', columns='gwas_category_count', fill_value=0)
cat_df = cat_df.loc[cat_df['consequence'].isin(cat_wide_df.loc[cat_wide_df.sum(axis=1) >= 10].index)]

order = ['intergenic_variant', 'upstream_gene_variant', 'intron_variant', 'missense_variant', '3_prime_UTR_variant', 'downstream_gene_variant']
xticklabels = [x[:-8] for x in order]
cat_df.loc[cat_df['gwas_category_count'] == upper_var_gwas_cat_count, 'gwas_category_count'] = "≥5"
ax = seaborn.catplot(x="percentage", y="consequence", hue="gwas_category_count", kind="bar", data=cat_df, order=order, legend=False, palette="rocket_r", orient="h")

title = "Variant Effect Predictor"
xlabel = "Variant consequence"
ylabel = "Variants [%]"
legend_loc = 'upper right'
label_fontsize = 16
tick_fontsize = 12
ax.set_yticklabels(xticklabels)

plt.grid(axis='x')
plt.legend(loc=legend_loc, title="GWAS cat. count")
plt.title(title, fontsize=label_fontsize)
plt.xlabel(xlabel, fontsize=label_fontsize)
plt.xticks(fontsize=tick_fontsize, rotation=0)
plt.yticks(fontsize=tick_fontsize, rotation=0)
plt.ylabel(ylabel, fontsize=label_fontsize)

plt.tight_layout()
plt.savefig(png_path, dpi=dpi)
plt.close()