import os
import pandas
import pathlib
import seaborn
import sys

import sqlalchemy
from matplotlib import pyplot as plt

from gwas2eqtl_pleiotropy.constants import label_fontsize

#%%
help_cmd_str = "todo"
try:
    snp_pp_h4 = float(sys.argv[1])
    manuscript_pleio_cutoff = int(sys.argv[2])
    db_url = sys.argv[3]
    count_per_rsid_gwas_egene_etissue_ods = sys.argv[4]
    count_per_rsid_gwas_egene_etissue_corr_png = sys.argv[5]
    if len(sys.argv) > 6:
        print("""Two many arguments!
        {}""".format(help_cmd_str))
        sys.exit(1)
except IndexError:
    print("""Argument missing!
    {}""".format(help_cmd_str))
    sys.exit(1)

outdir_path = os.path.dirname(count_per_rsid_gwas_egene_etissue_ods)
pathlib.Path(outdir_path).mkdir(parents=True, exist_ok=True)

columns = ['chrom', 'pos19', 'pos38', 'cytoband', 'rsid', 'ref', 'alt', 'gwas_trait_ontology_term',
           'gwas_category_ontology_term', 'eqtl_beta', 'eqtl_gene_id', 'eqtl_gene_symbol', 'gwas_id',
           'eqtl_id', 'etissue_category_term', 'pubmed_count', 'domains_watanabe2019']
cols_str = ','.join(columns)
sql = 'select distinct {} from colocpleio where snp_pp_h4>={}'.format(cols_str, snp_pp_h4)
engine = sqlalchemy.create_engine(db_url)
with engine.begin() as conn:
    coloc_df = pandas.read_sql(sqlalchemy.text(sql), con=conn)

#%% definition of variant for aggregation
variant_def_lst = ['chrom', 'cytoband', 'pos19', 'pos38', 'rsid', 'ref', 'alt']

#%% per rsid, aggregate gwas categories
gwas_df = coloc_df[variant_def_lst + ['gwas_category_ontology_term']].drop_duplicates()
agg_dic = {'gwas_category_ontology_term': lambda x: x.tolist()}
gwas_df = gwas_df.groupby(variant_def_lst).agg(agg_dic).reset_index()
gwas_df['gwas_category_count'] = gwas_df['gwas_category_ontology_term'].apply(len)
gwas_df.rename({'gwas_category_ontology_term': 'gwas_category_lst'}, axis=1, inplace=True)

#%% per rsid, aggregate gwas categories
gwas_ontology_df = coloc_df[variant_def_lst + ['gwas_trait_ontology_term']].drop_duplicates()
agg_dic = {'gwas_trait_ontology_term': lambda x: x.tolist()}
gwas_ontology_df = gwas_ontology_df.groupby(variant_def_lst).agg(agg_dic).reset_index()
gwas_ontology_df['gwas_trait_count'] = gwas_ontology_df['gwas_trait_ontology_term'].apply(len)

#%% per rsid, aggregate eqtl genes
egene_df = coloc_df[variant_def_lst + ['eqtl_gene_id', 'eqtl_gene_symbol']].drop_duplicates()
agg_dic = {'eqtl_gene_id': lambda x: x.tolist(), 'eqtl_gene_symbol': lambda x: x.tolist()}
egene_df = egene_df.groupby(variant_def_lst).agg(agg_dic).reset_index()
egene_df['eqtl_gene_id_count'] = egene_df['eqtl_gene_id'].apply(len)
egene_df.rename({'eqtl_gene_id': 'egene_lst', 'eqtl_gene_symbol': 'eqtl_gene_symbol_lst', 'eqtl_gene_id_count': 'egene_count'}, axis=1, inplace=True)

#%% per rsid, aggregate eqtl tissues
etissue_df = coloc_df[variant_def_lst + ['etissue_category_term']].drop_duplicates()
agg_dic = {'etissue_category_term': lambda x: x.tolist()}
etissue_df = etissue_df.groupby(variant_def_lst).agg(agg_dic).reset_index()
etissue_df['etissue_category_term_count'] = etissue_df['etissue_category_term'].apply(len)
etissue_df.rename({'etissue_category_term': 'etissue_category_term_lst'}, axis=1, inplace=True)

#%% per rsid, aggregate eqtl gene pubmed count
gene_pubmed_marker_df = coloc_df[variant_def_lst + ['eqtl_gene_id', 'eqtl_gene_symbol', 'pubmed_count', 'domains_watanabe2019']].drop_duplicates()
gene_pubmed_marker_df.rename({'eqtl_gene_id': 'eqtl_gene_marker_id', 'eqtl_gene_symbol': 'eqtl_gene_marker_symbol'}, axis=1, inplace=True)
gene_pubmed_marker_df['pubmed_count'] = gene_pubmed_marker_df['pubmed_count'].fillna(0)
gene_pubmed_marker_df.sort_values('pubmed_count', ascending=False, inplace=True)
gene_pubmed_marker_df.drop_duplicates(variant_def_lst, inplace=True)

#%% merge all aggregation columns
m_df = pandas.merge(gwas_df, gwas_ontology_df, on=variant_def_lst)
m_df = pandas.merge(m_df, egene_df, on=variant_def_lst)
m_df = pandas.merge(m_df, etissue_df, on=variant_def_lst)
m_df = pandas.merge(m_df, gene_pubmed_marker_df, on=variant_def_lst)
m_df.sort_values(['gwas_category_count', 'gwas_trait_count', 'chrom', 'pos38', 'rsid'], ascending=[False, False, True, True, True], inplace=True)
columns = ['chrom', 'pos19', 'pos38', 'cytoband', 'rsid', 'ref', 'alt', 'eqtl_gene_marker_symbol',
           'pubmed_count', 'gwas_category_count', 'gwas_trait_count',
           'gwas_category_lst', 'gwas_trait_ontology_term', 'egene_lst', 'eqtl_gene_symbol_lst', 'egene_count',
           'etissue_category_term_lst', 'etissue_category_term_count', 'eqtl_gene_marker_id', 'domains_watanabe2019']
m_df = m_df[columns]

m_df['gwas_category_lst'] = m_df['gwas_category_lst'].apply(lambda x: ";".join(sorted([i for i in x if not i is None])))
m_df['gwas_trait_ontology_term'] = m_df['gwas_trait_ontology_term'].apply(lambda x: ";".join(sorted([i for i in x if not i is None])))
m_df['egene_lst'] = m_df['egene_lst'].apply(lambda x: ";".join(sorted([i for i in x if not i is None])))
m_df['eqtl_gene_symbol_lst'] = m_df['eqtl_gene_symbol_lst'].apply(lambda x: ";".join(sorted([i for i in x if not i is None])))
m_df['etissue_category_term_lst'] = m_df['etissue_category_term_lst'].apply(lambda x: ";".join(sorted([i for i in x if not i is None])))

with pandas.ExcelWriter(count_per_rsid_gwas_egene_etissue_ods, engine="odf") as writer:
    m_df.to_excel(writer, index=False)

#%% gwas variant pleiotropy for MS
ms_df = m_df.copy()
ms_df = ms_df.drop_duplicates('cytoband', keep='first')
ms_df = ms_df.loc[ms_df['gwas_category_count'] >= manuscript_pleio_cutoff]
columns = ['chrom', 'pos38', 'cytoband', 'rsid', 'eqtl_gene_marker_symbol', 'gwas_category_lst']
ms_df = ms_df[columns]
# format output
# ms_df['gwas_category_lst'] = ms_df['gwas_category_lst'].apply(lambda x: '; '.join(x))
ms_df['pos38'] = ms_df['pos38'].apply(lambda x : '{0:,}'.format(x))
ms_df['rsid'] = 'rs' + ms_df['rsid'].astype(str)
ms_df['gwas_category_lst'] = ms_df['gwas_category_lst'].str.replace(';', '; ')
tsv_path = os.path.join(outdir_path, "count_per_rsid_gwas_ms.tsv")
ms_df.to_csv(tsv_path, sep="\t", index=False)

#%%
m2df = m_df[['gwas_category_count', 'etissue_category_term_count', 'egene_count']].copy()
m2df.rename({'gwas_category_count': 'GWAS. cat. cnt.', 'etissue_category_term_count': 'eTissue cat. cnt.', 'egene_count': 'eGene cnt.', }, axis=1, inplace=True)
corr = m2df.corr(method='spearman')
plt.subplots_adjust(left=0.2, right=0.8, top=0.9, bottom=0.)
ax = seaborn.heatmap(corr, annot=True, xticklabels=False)
ax.set_yticklabels(ax.get_ymajorticklabels(), fontsize=12, rotation=45)
plt.title("Spearman correlation", fontsize=16)
plt.savefig(count_per_rsid_gwas_egene_etissue_corr_png)
plt.clf()
plt.close()

#%########################################### watanabe 2019 category count
m2df = m_df[['rsid', 'ref', 'alt', 'gwas_category_count', 'domains_watanabe2019']].copy()

order = [*range(1, m2df['gwas_category_count'].max()+1)]
ax = seaborn.boxplot(x='gwas_category_count', y='domains_watanabe2019', data=m2df, order=order, palette="rocket_r")
plt.ylabel("Watanabe cat. cnt.", fontsize=label_fontsize)
plt.xlabel("GWAS category count", fontsize=label_fontsize)
png_path = os.path.join(outdir_path, 'watanabe_cat_count.png')
plt.savefig(png_path)
plt.clf()
plt.close()

#%########################################### watanabe 2019 category count
m2df_cat_count_df = m2df.groupby('gwas_category_count').agg({'domains_watanabe2019': len})
m2df_watanabe2019_cat_count_df = m2df.loc[~m2df['domains_watanabe2019'].isna()].groupby('gwas_category_count').agg({'domains_watanabe2019': len})
m2df_watanabe2019_cat_count_df = m2df_watanabe2019_cat_count_df.merge(m2df_cat_count_df, left_index=True, right_index=True)
m2df_watanabe2019_cat_count_df.columns = ['watanabe_count', 'count']
m2df_watanabe2019_cat_count_df['watanabe_perc'] = (m2df_watanabe2019_cat_count_df['watanabe_count'] / m2df_watanabe2019_cat_count_df['count'] * 100).astype(int)

order = [*range(1, m2df['gwas_category_count'].max()+1)]
ax = seaborn.barplot(x=m2df_watanabe2019_cat_count_df.index, y='watanabe_perc', data=m2df_watanabe2019_cat_count_df, order=order, palette="rocket_r")
plt.title("Part of Watanabe in our study", fontsize=label_fontsize)
plt.ylabel("Percentage", fontsize=label_fontsize)
plt.xlabel("GWAS category count", fontsize=label_fontsize)
png_path = os.path.join(outdir_path, 'watanabe_percentage.png')
plt.savefig(png_path)
plt.clf()
plt.close()


#%########################################### bed files, flanking=0
# bed files of variants splitted by gwas categories
flank = 10
bed_df = gwas_df.copy()
bed_df['chrom'] = 'chr' + bed_df['chrom'].astype(str)
bed_df['start'] = bed_df['pos38'] - 1 - flank
bed_df['end'] = bed_df['pos38'] + flank
bed_df['rsid'] = 'rs' + bed_df['rsid'].astype(str)
bed_df = bed_df[['chrom', 'start', 'end', 'rsid', 'gwas_category_count']]

for gwas_category_count in sorted(bed_df['gwas_category_count'].unique()):
    pleio_bed_path = os.path.join(outdir_path, "eqtl_pleio_{}_flank_{}_hg38.bed".format(gwas_category_count, flank))
    bed_pleio_df = bed_df.loc[bed_df['gwas_category_count'] == gwas_category_count, ]
    bed_pleio_df = bed_pleio_df.sort_values(by=['chrom', 'start', 'end'])
    bed_pleio_df.to_csv(pleio_bed_path, sep="\t", index=False, header=False)

#%########################################### bed files, flanking=0
# bed files of variants splitted by gwas categories
flank = 0
bed_df = gwas_df.copy()
bed_df['chrom'] = 'chr' + bed_df['chrom'].astype(str)
bed_df['start'] = bed_df['pos38'] - 1 - flank
bed_df['end'] = bed_df['pos38'] + flank
bed_df['rsid'] = 'rs' + bed_df['rsid'].astype(str)
bed_df = bed_df[['chrom', 'start', 'end', 'rsid', 'gwas_category_count']]

for gwas_category_count in sorted(bed_df['gwas_category_count'].unique()):
    pleio_bed_path = os.path.join(outdir_path, "eqtl_pleio_{}_flank_{}_hg38.bed".format(gwas_category_count, flank))
    bed_pleio_df = bed_df.loc[bed_df['gwas_category_count'] == gwas_category_count, ]
    bed_pleio_df = bed_pleio_df.sort_values(by=['chrom', 'start', 'end'])
    bed_pleio_df.to_csv(pleio_bed_path, sep="\t", index=False, header=False)
