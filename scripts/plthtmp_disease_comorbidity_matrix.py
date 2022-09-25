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
gwas_cat_df = pandas.read_excel(gwas_cat_ods_path, engine="odf", usecols=['id', 'trait', 'manual_category'])
gwas_cat_df.rename({'id': 'gwas_id', 'manual_category': 'category'}, axis=1, inplace=True)

#%%
# import pdb; pdb.set_trace()

# iris = seaborn.load_dataset("iris")

# species = iris.pop("species")
# import pdb; pdb.set_trace()
corr_df = corr_df.merge(gwas_cat_df, on='gwas_id')
corr_df.set_index(['gwas_id', 'trait', 'category'], inplace=True)

#%%
# Prepare a vector of color mapped to the 'cyl' column
# lut2 = dict(zip(corr_df.index.get_level_values('category').unique(), "tab10"))
# row_colors2 = corr_df.index.get_level_values('category').map(lut2)

#%%
seaborn.set(font_scale=0.3)
g = seaborn.clustermap(corr_df, cmap='coolwarm')

plt.tight_layout()
png_path = os.path.join(outdir_path, "blood.png")
plt.savefig(png_path, dpi=600)
plt.clf()
plt.close()

#%%
clustered_df = g.__dict__['data2d']
clustered_tsv_path = os.path.join(outdir_path, 'clustered.tsv')
clustered_df.to_csv(clustered_tsv_path, sep='\t', index=True, header=True)
