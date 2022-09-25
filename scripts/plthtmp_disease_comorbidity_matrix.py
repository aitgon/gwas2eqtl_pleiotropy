import os
import pandas
import pathlib
import sys
import seaborn
from gwas2eqtl_pleiotropy.Logger import Logger
from gwas2eqtl_pleiotropy.constants import seaborn_theme_dic, dpi
from matplotlib import pyplot as plt
import scipy.spatial as sp, scipy.cluster.hierarchy as hc
from matplotlib.pyplot import gcf
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

gwas_category_df = gwas_cat_df[['category.1', 'code.1', 'category3', 'code2']].dropna(axis=0).drop_duplicates()
gwas_category_df.rename({'category.1': 'category', 'code.1': 'code'}, axis=1, inplace=True)
gwas_category_df.set_index('code', inplace=True, verify_integrity=True)
gwas_id_df = gwas_id_df.merge(gwas_category_df, left_on='code', right_index=True, how='inner')
gwas_id_df = gwas_id_df.merge(dis_df, left_index=True, right_index=True)[['trait', 'category3']]
# import pdb; pdb.set_trace()
#%%
category_lst = gwas_id_df['category3'].unique()
category_lst_len = len(category_lst)
nb_lst = int(category_lst_len/2)
category1_lst = gwas_id_df['category3'].unique()[:nb_lst]
category2_lst = gwas_id_df['category3'].unique()[nb_lst:]
# import pdb; pdb.set_trace()
#%%
gwas_id_df['subset1'] = 'Other'
mask1 = gwas_id_df['category3'].isin(category1_lst)
gwas_id_df.loc[mask1, 'subset1'] = gwas_id_df.loc[mask1, 'category3']

gwas_id_df['subset2'] = 'Other'
gwas_id_df.loc[~mask1, 'subset2'] = gwas_id_df.loc[~mask1, 'category3']
# import pdb; pdb.set_trace()
#%%
annotation_df = dis_df.merge(gwas_id_df, left_index=True, right_index=True, how='left')[['subset1', 'subset2']]

# Label 1
subset1_labels = annotation_df["subset1"]
subset1_pal = seaborn.color_palette(palette='bright', n_colors=subset1_labels.unique().size)
subset1_lut = dict(zip(map(str, sorted(subset1_labels.unique())), subset1_pal))
subset1_colors = pandas.Series(subset1_labels, index=annotation_df.index).map(subset1_lut)

# Label 2
subset2_labels = annotation_df["subset2"]
subset2_pal = seaborn.color_palette(palette='tab10', n_colors=subset2_labels.unique().size)
subset2_lut = dict(zip(map(str, sorted(subset2_labels.unique())), subset2_pal))
subset2_colors = pandas.Series(subset2_labels, index=annotation_df.index).map(subset2_lut)

# import pdb; pdb.set_trace()
network_node_colors = pandas.DataFrame(subset1_colors).join(pandas.DataFrame(subset2_colors))

#%%
# category_annot_lst = dis_df.merge(gwas_id_df, left_index=True, right_index=True, how='left')['subset1'].tolist()
# lut = dict(zip(category_annot_lst, seaborn.color_palette(palette='bright', n_colors=len(category_annot_lst))))
# row_colors = pandas.DataFrame(category_annot_lst)[0].map(lut)
#
# #%%
# category_annot_lst = dis_df.merge(gwas_id_df, left_index=True, right_index=True, how='left')['subset2'].tolist()
# lut2 = dict(zip(category_annot_lst, seaborn.color_palette(palette='bright', n_colors=len(category_annot_lst))))
# row_colors2 = pandas.DataFrame(category_annot_lst)[0].map(lut2)

#%%
# seaborn.set(font_scale=0.2)

linkage = hc.linkage(sp.distance.squareform(dis_df), method='average')

clustermap_args_dic = {}
clustermap_args_dic['cmap'] = 'mako'
clustermap_args_dic['row_colors'] = network_node_colors
clustermap_args_dic['col_colors'] = network_node_colors
clustermap_args_dic['row_linkage'] = linkage
clustermap_args_dic['col_linkage'] = linkage
clustermap_args_dic['row_cluster'] = False
clustermap_args_dic['col_cluster'] = False
clustermap_args_dic['xticklabels'] = False
clustermap_args_dic['yticklabels'] = False
# clustermap_args_dic['cbar_pos'] = None
# clustermap_args_dic['dendrogram_ratio'] = (0.05, 0.05)
g = seaborn.clustermap(dis_df, **clustermap_args_dic)
# add legends
g.cax.set_visible(False)

for label in subset1_labels.unique():
    g.ax_col_dendrogram.bar(0, 0, color=subset1_lut[label], label=label, linewidth=0);
l1 = g.ax_col_dendrogram.legend(title='subset1', loc="upper left", bbox_to_anchor=(0.0, 1), ncol=3, bbox_transform=gcf().transFigure)

for label in subset2_labels.unique():
    g.ax_row_dendrogram.bar(0, 0, color=subset2_lut[label], label=label, linewidth=0);
l2 = g.ax_row_dendrogram.legend(title='subset2', loc="upper left", bbox_to_anchor=(0.0, 0.9), ncol=3, bbox_transform=gcf().transFigure)

plt.subplots_adjust(top=0.95, left=0.00, right=0.9, bottom=0.05)
# g.fig.subplots_adjust(right=0.7)
# plt.tight_layout()
# png_path = os.path.join(outdir_path, "blood.png")
plt.savefig(htmp_disease_corr_png_path, dpi=600)
plt.clf()
plt.close()

#%%
clustered_df = g.__dict__['data2d']
clustered_tsv_path = os.path.join(outdir_path, 'clustered.tsv')
clustered_df.to_csv(clustered_tsv_path, sep='\t', index=True, header=True)
