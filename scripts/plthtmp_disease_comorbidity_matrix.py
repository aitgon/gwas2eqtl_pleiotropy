import math
import os
import pandas
import pathlib
import sys
import seaborn

from gwas2eqtl_pleiotropy.Logger import Logger
from gwas2eqtl_pleiotropy.constants import seaborn_theme_dic, dpi
from matplotlib import pyplot as plt
from matplotlib.pyplot import gcf
from sqlalchemy import create_engine

#%%
seaborn.set_theme(**seaborn_theme_dic)


#%%
help_cmd_str = "todo"
try:
    sa_url = sys.argv[1]  # annotate
    disease_corr_tsv_path = sys.argv[2]
    gwas_metadata_ods_path = sys.argv[3]
    htmp_disease_corr_png_path = sys.argv[4]
    if len(sys.argv) > 5:
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
corr_count_ser = (corr_df.abs() >= 0.005).sum(axis=1)
# max number of correlation
corr_count_ser = corr_count_ser.sort_values(ascending=False)[0:30]

mask = corr_df.index.isin(corr_count_ser.index)
corr_df = corr_df.loc[mask]
corr_df = corr_df[corr_df.columns[mask]]

#%% create distance matrix
dis_df = 1 - corr_df   # distance matrix
gwas_metadata_df = pandas.read_sql('select distinct gwas_id, gwas_trait_ontology_term, batch, pmid , gwas_category_ontology_term from colocpleio', con=create_engine(sa_url))

#%%
# gwas_metadata_df = pandas.read_excel(gwas_metadata_ods_path, engine="odf", usecols=['gwas_id', 'trait', 'gwas_category'])

#%%
# import pdb; pdb.set_trace()

#%%
gwas_metadata_df.set_index('gwas_id', inplace=True, verify_integrity=True)
# annotation_df = dis_df.merge(gwas_metadata_df, left_index=True, right_index=True, how='left')[['trait', 'gwas_category']]
annotation_df = dis_df.merge(gwas_metadata_df, left_index=True, right_index=True, how='left')[['batch', 'pmid', 'gwas_trait_ontology_term', 'gwas_category_ontology_term']]

# pmid_df = pandas.read_sql('select distinct gwas_id, pmid from colocpleio', con=create_engine(sa_url), index_col='gwas_id')
# pmid_df = pmid_df.loc[annotation_df.index.tolist(), ]
# import pdb; pdb.set_trace()
# annotation_df = annotation_df.merge(pmid_df, left_index=True, right_index=True)
annotation_df['pmid'] = annotation_df['pmid'].replace(math.nan, 0)
annotation_df['pmid'] = annotation_df['pmid'].astype(int)

#%%
annotation_df = annotation_df.drop_duplicates(['pmid', 'gwas_trait_ontology_term'])
dis_df = dis_df.loc[annotation_df.index]
dis_df = dis_df[annotation_df.index]

#%%
annotation_df['trait'] = annotation_df['batch'] + '_' + annotation_df['pmid'].astype(str) + '_' + annotation_df['gwas_trait_ontology_term']

# pmid_trait_dupli_mask = annotation_df['trait'].duplicated(keep='first')
# pmid_trait_uniq_mask = ~annotation_df['trait'].duplicated(keep=False)
# pmid_trait_mask = pmid_trait_dupli_mask + pmid_trait_uniq_mask
# dis_df = dis_df.loc[pmid_trait_mask, pmid_trait_mask]
# annotation_df = annotation_df.loc[pmid_trait_mask, ]

# dataset_a = (annotation_df.index).to_series().str.split('-', expand=True)[0]
# dataset_b = (annotation_df.index).to_series().str.split('-', expand=True)[1]
# annotation_df['trait'] = dataset_a + '_' + dataset_b + '_' + annotation_df['trait']
# import pdb; pdb.set_trace()
# Label 1
category_labels = annotation_df["gwas_category_ontology_term"]
category_pal = seaborn.color_palette(palette='bright', n_colors=category_labels.unique().size)
subset1_lut = dict(zip(map(str, sorted(category_labels.unique())), category_pal))
subset1_colors = pandas.Series(category_labels, index=annotation_df.index).map(subset1_lut)
network_node_colors = pandas.DataFrame(subset1_colors)

#%%
# linkage = hc.linkage(sp.distance.squareform(dis_df), method='average')

clustermap_args_dic = {}
clustermap_args_dic['cmap'] = 'mako'
clustermap_args_dic['row_colors'] = network_node_colors
clustermap_args_dic['col_colors'] = network_node_colors
clustermap_args_dic['xticklabels'] = False

# import pdb; pdb.set_trace()
clustermap_args_dic['yticklabels'] = [s[0:40] for s in annotation_df['trait'].tolist()]
g = seaborn.clustermap(dis_df, **clustermap_args_dic)

seaborn.set(font_scale=1)
g.cax.set_visible(False)
# import pdb; pdb.set_trace()
g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_ymajorticklabels(), fontsize=12)
g.ax_heatmap.set(ylabel='Trait')
#
for label in category_labels.unique():
    g.ax_col_dendrogram.bar(0, 0, color=subset1_lut[label], label=label, linewidth=0);
l1 = g.ax_col_dendrogram.legend(title='Class', loc="upper left", bbox_to_anchor=(0.05, 0.95), ncol=3, bbox_transform=gcf().transFigure)

plt.subplots_adjust(top=1., right=0.6, bottom=0.05, left=0.)
plt.savefig(htmp_disease_corr_png_path, dpi=600)
plt.clf()
plt.close()

#%%
clustered_df = g.__dict__['data2d']
clustered_tsv_path = os.path.join(outdir_path, 'clustered.tsv')
clustered_df.to_csv(clustered_tsv_path, sep='\t', index=True, header=True)
