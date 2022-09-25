import os
import pandas
import pathlib
import sys
import seaborn
from gwas2eqtl_pleiotropy.Logger import Logger
from gwas2eqtl_pleiotropy.constants import seaborn_theme_dic, dpi
from matplotlib import pyplot as plt

#%%
seaborn.set_theme(**seaborn_theme_dic)

#%%
help_cmd_str = "todo"
try:
    disease_corr_tsv_path = sys.argv[1]
    gwas_cat_ods_path = sys.argv[2]
    htmp_disease_corr_png_path = sys.argv[3]
    if len(sys.argv) > 4:
        print("""Two many arguments!
        {}""".format(help_cmd_str))
        sys.exit(1)
except IndexError:
    print("""Argument missing!
    {}""".format(help_cmd_str))
    sys.exit(1)

outdir_path = os.path.dirname(htmp_disease_corr_png_path)
pathlib.Path(outdir_path).mkdir(parents=True, exist_ok=True)

#%%
Logger.info("Reading {}".format(disease_corr_tsv_path))
corr_df = pandas.read_csv(disease_corr_tsv_path, sep="\t")

#%%
gwas_cat_df = pandas.read_excel(gwas_cat_ods_path, engine="odf")
gwas_id_df = gwas_cat_df[['id', 'trait', 'category', 'code']].drop_duplicates()
gwas_category_df = gwas_cat_df[['category.1', 'code.1', 'category2', 'code2']].dropna(axis=0).drop_duplicates()

m_df = gwas_id_df.merge(gwas_category_df, left_on='code', right_on='code.1', how='outer')
# import pdb; pdb.set_trace()

gwas_id_df = gwas_id_df.merge(gwas_category_df, left_on='code', right_on='code.1', how='outer')

gwas_id_df.rename({'id': 'gwas_id'}, axis=1, inplace=True)
gwas_id_df = gwas_id_df[['gwas_id', 'trait', 'category2']]


#%%
corr_df = corr_df.merge(gwas_id_df, on='gwas_id')
corr_df.set_index(['gwas_id', 'trait', 'category2'], inplace=True)

#%%
# Prepare a vector of color mapped to the 'cyl' column
# import pdb; pdb.set_trace()
category2_lst = corr_df.index.get_level_values('category2').tolist()

# lut2 = dict(zip(category2_lst, seaborn.color_palette(palette="deep", n_colors=len(category2_lst), as_cmap=True)))
# row_colors2 = corr_df.index.get_level_values('category2').map(lut2)
lut = dict(zip(set(category2_lst), seaborn.hls_palette(len(set(category2_lst)), l=0.5, s=0.8)))
row_colors = pandas.DataFrame(category2_lst)[0].map(lut)

#%%
seaborn.set(font_scale=0.3)
# import pdb; pdb.set_trace()
g = seaborn.clustermap(corr_df, cmap='coolwarm', row_colors=[row_colors])

plt.tight_layout()
png_path = os.path.join(outdir_path, "blood.png")
plt.savefig(png_path, dpi=600)
plt.clf()
plt.close()

#%%
clustered_df = g.__dict__['data2d']
clustered_tsv_path = os.path.join(outdir_path, 'clustered.tsv')
clustered_df.to_csv(clustered_tsv_path, sep='\t', index=True, header=True)
