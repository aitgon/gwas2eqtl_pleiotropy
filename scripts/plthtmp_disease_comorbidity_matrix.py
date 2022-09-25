import os
import pandas
import pathlib
import sys
import seaborn
from gwas2eqtl_pleiotropy.Logger import Logger
from gwas2eqtl_pleiotropy.constants import seaborn_theme_dic, dpi
from matplotlib import pyplot as plt
import scipy.spatial as sp, scipy.cluster.hierarchy as hc

from matplotlib.patches import Patch

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
corr_df = pandas.read_csv(disease_corr_tsv_path, sep="\t", index_col='gwas_id')

corr_df = corr_df.loc[corr_df.sum(axis=1) > 0, ]
corr_df = corr_df[corr_df.columns[corr_df.sum(axis=0) > 0]]

dis_df = 1 - corr_df   # distance matrix

#%%
gwas_cat_df = pandas.read_excel(gwas_cat_ods_path, engine="odf")
gwas_id_df = gwas_cat_df[['id', 'trait', 'category', 'code']].drop_duplicates()
gwas_id_df.set_index('id', inplace=True, verify_integrity=True)

gwas_category_df = gwas_cat_df[['category.1', 'code.1', 'category4', 'code2']].dropna(axis=0).drop_duplicates()
gwas_category_df.rename({'category.1': 'category', 'code.1': 'code'}, axis=1, inplace=True)
gwas_category_df.set_index('code', inplace=True, verify_integrity=True)
gwas_id_df = gwas_id_df.merge(gwas_category_df, left_on='code', right_index=True, how='inner')
gwas_id_df = gwas_id_df[['trait', 'category4']]

# import pdb; pdb.set_trace()
#%%
# dis_df = dis_df.merge(gwas_id_df, on='gwas_id')
# dis_df.set_index(['gwas_id', 'trait', 'category4'], inplace=True)

#%%
# Prepare a vector of color mapped to the 'cyl' column
# import pdb; pdb.set_trace()
category3_lst = dis_df.merge(gwas_id_df, left_index=True, right_index=True, how='left')['category4'].tolist()
# lut2 = dict(zip(category3_lst, seaborn.color_palette(palette="deep", n_colors=len(category3_lst), as_cmap=True)))
# row_colors2 = corr_df.index.get_level_values('category4').map(lut2)
# import pdb; pdb.set_trace()
lut = dict(zip(category3_lst, seaborn.color_palette(palette='Paired', n_colors=len(category3_lst))))
row_colors = pandas.DataFrame(category3_lst)[0].map(lut)
# print(zip(set(category3_lst), seaborn.color_palette(palette='bright', n_colors=len(category3_lst))))
# print(lut)
# print(row_colors)
# import pdb; pdb.set_trace()

#%%
# import pdb; pdb.set_trace()
seaborn.set(font_scale=0.2)

linkage = hc.linkage(sp.distance.squareform(dis_df), method='average')
# import pdb; pdb.set_trace()

# clustermap_args_dic = {'dendrogram_ratio': 0.1, 'colors_ratio': 0.1}
clustermap_args_dic = {}
clustermap_args_dic['cmap'] = 'mako'
clustermap_args_dic['row_colors'] = [row_colors]
clustermap_args_dic['col_colors'] = [row_colors]
clustermap_args_dic['row_linkage'] = linkage
clustermap_args_dic['col_linkage'] = linkage
clustermap_args_dic['row_cluster'] = False
clustermap_args_dic['col_cluster'] = False
# clustermap_args_dic['cbar_pos'] = None
# clustermap_args_dic['dendrogram_ratio'] = (0.05, 0.05)
g = seaborn.clustermap(dis_df, **clustermap_args_dic)
# g = seaborn.clustermap(dis_df, cmap='coolwarm', row_linkage=linkage, col_linkage=linkage)

handles = [Patch(facecolor=lut[name]) for name in lut]
plt.legend(handles, lut, title='Category3', bbox_transform=plt.gcf().transFigure, loc='best')

plt.subplots_adjust(top=0.95, left=0.05, right=0.95, bottom=0.05)
# plt.tight_layout()
# png_path = os.path.join(outdir_path, "blood.png")
plt.savefig(htmp_disease_corr_png_path, dpi=600)
plt.clf()
plt.close()

#%%
clustered_df = g.__dict__['data2d']
clustered_tsv_path = os.path.join(outdir_path, 'clustered.tsv')
clustered_df.to_csv(clustered_tsv_path, sep='\t', index=True, header=True)
