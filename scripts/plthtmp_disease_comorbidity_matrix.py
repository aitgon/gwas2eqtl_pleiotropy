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

#%%
seaborn.set_theme(**seaborn_theme_dic)



#%%
help_cmd_str = "todo"
try:
    disease_corr_tsv_path = sys.argv[1]
    gwas_metadata_ods_path = sys.argv[2]
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

#%% filter row with at least a number of correlations >0.1 larger than 30
corr_count_ser = (corr_df.abs() >= 0.05).sum(axis=1)
# max number of correlation 80
corr_count_ser = corr_count_ser.sort_values(ascending=False)[0:50]
# mask = (corr_df >= 0.04).sum(axis=1) > 10
mask = corr_df.index.isin(corr_count_ser.index)
corr_df = corr_df.loc[mask]
corr_df = corr_df[corr_df.columns[mask]]

#%% create distance matrix
dis_df = 1 - corr_df   # distance matrix

#%%
gwas_metadata_df = pandas.read_excel(gwas_metadata_ods_path, engine="odf", usecols=['gwas_id', 'trait', 'class_pleiotropy'])

#%%
gwas_metadata_df.set_index('gwas_id', inplace=True, verify_integrity=True)
annotation_df = dis_df.merge(gwas_metadata_df, left_index=True, right_index=True, how='left')[['trait', 'class_pleiotropy']]

# Label 1
class_labels = annotation_df["class_pleiotropy"]
class_pal = seaborn.color_palette(palette='bright', n_colors=class_labels.unique().size)
subset1_lut = dict(zip(map(str, sorted(class_labels.unique())), class_pal))
subset1_colors = pandas.Series(class_labels, index=annotation_df.index).map(subset1_lut)
network_node_colors = pandas.DataFrame(subset1_colors)

#%%
# import pdb; pdb.set_trace()
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

clustermap_args_dic['yticklabels'] = [s[0:40] for s in annotation_df['trait'].tolist()]

g = seaborn.clustermap(dis_df, **clustermap_args_dic)

seaborn.set(font_scale=1)
g.cax.set_visible(False)
g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_ymajorticklabels(), fontsize=10)
g.ax_heatmap.set(ylabel='Trait')

for label in class_labels.unique():
    g.ax_col_dendrogram.bar(0, 0, color=subset1_lut[label], label=label, linewidth=0);
l1 = g.ax_col_dendrogram.legend(title='Class', loc="upper left", bbox_to_anchor=(0.05, 0.95), ncol=4, bbox_transform=gcf().transFigure)

plt.subplots_adjust(top=1., right=0.6, bottom=0.05, left=-0.15)
plt.savefig(htmp_disease_corr_png_path, dpi=600)
plt.clf()
plt.close()

#%%
clustered_df = g.__dict__['data2d']
clustered_tsv_path = os.path.join(outdir_path, 'clustered.tsv')
clustered_df.to_csv(clustered_tsv_path, sep='\t', index=True, header=True)
