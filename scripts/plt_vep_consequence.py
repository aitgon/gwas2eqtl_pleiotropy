import os
import pandas
import seaborn
from matplotlib import pyplot as plt

#%%
count_per_rsid_gwas_tsv_path = "out/gwas413/genome/5e-08/1000000/cmpt_count_per_rsid.py/count_per_rsid_gwas.tsv"
vep0_path = "out/gwas413/genome/5e-08/1000000/cmpt_vep_consequence.py/vep_pleio{}.tsv"
upper_var_gwas_cat_count = 5


#%%
df = pandas.read_csv(count_per_rsid_gwas_tsv_path, sep="\t")
df.loc[df['gwas_category_count'] >= upper_var_gwas_cat_count, 'gwas_category_count'] = upper_var_gwas_cat_count
cat_df = pandas.DataFrame({'percentage': None, 'consequence': None, 'gwas_category_count': None}, index=[])

for pleio in range(1, upper_var_gwas_cat_count+1):
    vep_path = vep0_path.format(pleio)
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

cat_df.loc[cat_df['gwas_category_count'] == upper_var_gwas_cat_count, 'gwas_category_count'] = "≥5"
seaborn.catplot(x="consequence", y="percentage", hue="gwas_category_count", kind="bar", data=cat_df, order=order, legend=False)
plt.xticks(fontsize=12, rotation=45)
plt.grid(axis='y')
plt.legend(loc='upper right')

plt.tight_layout()
plt.savefig('t.png')
plt.close()
