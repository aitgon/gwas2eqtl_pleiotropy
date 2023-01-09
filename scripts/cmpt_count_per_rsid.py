import os
import pandas
import pathlib
import seaborn
import sys

from matplotlib import pyplot as plt


#%%
help_cmd_str = "todo"
try:
    snp_pp_h4 = float(sys.argv[1])
    max_gwas_category_count = int(sys.argv[2])
    manuscript_pleio_cutoff = int(sys.argv[3])
    url = sys.argv[4]
    count_per_rsid_gwas_tsv_path = sys.argv[5]
    count_per_rsid_gwas_egene_etissue_corr_png = sys.argv[6]
    if len(sys.argv) > 7:
        print("""Two many arguments!
        {}""".format(help_cmd_str))
        sys.exit(1)
except IndexError:
    print("""Argument missing!
    {}""".format(help_cmd_str))
    sys.exit(1)

outdir_path = os.path.dirname(count_per_rsid_gwas_tsv_path)
pathlib.Path(outdir_path).mkdir(parents=True, exist_ok=True)

sql = 'select * from colocpleio where snp_pp_h4>={}'.format(snp_pp_h4)
columns = ['chrom', 'pos38', 'cytoband', 'rsid', 'eqtl_beta', 'eqtl_gene_id', 'gwas_id', 'eqtl_id', 'etissue_category_term', 'pubmed_count']
coloc_df = pandas.read_sql(sql, con=url, columns=columns).drop_duplicates()

#%% definition of variant for aggregation
variant_def_lst = ['chrom', 'cytoband', 'pos38', 'rsid', 'alt']

#%% per rsid, aggregate gwas categories
gwas_df = coloc_df[variant_def_lst + ['gwas_category']].drop_duplicates()
agg_dic = {'gwas_category': lambda x: x.tolist()}
gwas_df = gwas_df.groupby(variant_def_lst).agg(agg_dic).reset_index()
gwas_df['gwas_category_count'] = gwas_df['gwas_category'].apply(len)
gwas_df['gwas_category_count'].unique()

#%% per rsid, aggregate eqtl genes
egene_df = coloc_df[variant_def_lst + ['eqtl_gene_id', 'eqtl_gene_symbol']].drop_duplicates()
agg_dic = {'eqtl_gene_id': ['size', (lambda x: (";".join(sorted(x))))], 'eqtl_gene_symbol': lambda x: ";".join(sorted([str(i) for i in x]))}
egene_df = egene_df.groupby(variant_def_lst).agg(agg_dic).reset_index()

#%% per rsid, aggregate eqtl gene pubmed count
gene_pubmed_df = coloc_df[variant_def_lst + ['eqtl_gene_id', 'eqtl_gene_symbol', 'pubmed_count']].drop_duplicates()
gene_pubmed_df['pubmed_count'] = gene_pubmed_df['pubmed_count'].fillna(0)
gene_pubmed_df.sort_values('pubmed_count', ascending=False, inplace=True)
gene_pubmed_df.drop_duplicates(variant_def_lst, inplace=True)

#%% per rsid, aggregate eqtl tissues
etissue_df = coloc_df[variant_def_lst + ['etissue_category_term']].drop_duplicates()
agg_dic = {'etissue_category_term': ['size', lambda x: ";".join(sorted(x))]}
etissue_df = etissue_df.groupby(variant_def_lst).agg(agg_dic).reset_index()

#%% merge all aggregation columns
m_df = pandas.merge(gwas_df, egene_df, on=variant_def_lst)
m_df = pandas.merge(m_df, gene_pubmed_df, on=variant_def_lst)
m_df = pandas.merge(m_df, etissue_df, on=variant_def_lst)
import pdb; pdb.set_trace()
m_df = m_df[['chrom', 'cytoband', 'pos38', 'rsid', 'gwas_category_count', 'gwas_category_lst', 'egene_count', 'eqtl_gene_symbol_lst', 'etissue_label_count', 'etissue_category_term_lst', 'egene_lst']]
m_df.sort_values(['gwas_category_count', 'chrom', 'pos38', 'rsid'], ascending=[False, True, True, True], inplace=True)
tsv_path = os.path.join(outdir_path, "count_per_rsid_gwas_egene_etissue.tsv")
m_df.to_csv(tsv_path, sep="\t", index=False)

#%%
etissue_df.columns = ['chrom', 'cytoband', 'pos38', 'rsid', 'etissue_label_count', 'etissue_category_term_lst']
egene_df.columns = ['chrom', 'cytoband', 'pos38', 'rsid', 'egene_count', 'egene_lst', 'eqtl_gene_symbol_lst']
gwas_df.columns = ['chrom', 'cytoband', 'pos38', 'rsid', 'gwas_category_count', 'gwas_category_lst']
gwas_df = gwas_df.sort_values(by=['gwas_category_count', 'chrom', 'pos38', 'gwas_category_lst', 'rsid'], ascending=[False, True, True, True, True])
egene_df = egene_df.sort_values(by=['egene_count', 'chrom', 'pos38', 'egene_lst', 'eqtl_gene_symbol_lst', 'rsid'], ascending=[False, True, True, True, True, True])
etissue_df = etissue_df.sort_values(by=['etissue_label_count', 'chrom', 'pos38', 'etissue_category_term_lst', 'rsid'], ascending=[False, True, True, True, True])
tsv_path = os.path.join(outdir_path, "count_per_rsid_egene.tsv")
egene_df.to_csv(tsv_path, sep="\t", index=False)
tsv_path = os.path.join(outdir_path, "count_per_rsid_etissue.tsv")
etissue_df.to_csv(tsv_path, sep="\t", index=False)
tsv_path = os.path.join(outdir_path, "count_per_rsid_gwas.tsv")
gwas_df.to_csv(tsv_path, sep="\t", index=False)

#%% gwas variant pleiotropy for MS
gwas_ms_df = gwas_df.copy()
gwas_ms_df = gwas_ms_df.drop_duplicates('cytoband', keep='first')
gwas_ms_df = gwas_ms_df.loc[gwas_ms_df['gwas_category_count'] >= manuscript_pleio_cutoff]
# format output
gwas_ms_df.drop(['gwas_category_count'], inplace=True, axis=1)
gwas_ms_df['gwas_category_lst'] = gwas_ms_df['gwas_category_lst'].str.replace(',', ', ')
gwas_ms_df['pos38']=gwas_ms_df['pos38'].apply(lambda x : '{0:,}'.format(x))
tsv_path = os.path.join(outdir_path, "count_per_rsid_gwas_ms.tsv")
gwas_ms_df.to_csv(tsv_path, sep="\t", index=False)



m2df = m_df[['gwas_category_count', 'etissue_label_count', 'egene_count']]
corr = m2df.corr(method='spearman')
plt.subplots_adjust(left=0.3, right=0.8, top=0.9, bottom=0.35)
seaborn.heatmap(corr, annot=True)
plt.title("Count corr. GWAS, gene, tissue", fontsize=16)
plt.savefig(count_per_rsid_gwas_egene_etissue_corr_png)
# plt.savefig("img.png")
plt.clf()
plt.close()

#%########################################### bed files, flanking=0
# bed files of variants splitted by gwas categories
flank = 10
variant_bed_df = gwas_df.copy()
variant_bed_df['chrom'] = 'chr' + variant_bed_df['chrom'].astype(str)
variant_bed_df['start'] = variant_bed_df['pos38'] - 1 - flank
variant_bed_df['end'] = variant_bed_df['pos38'] + flank
variant_bed_df = variant_bed_df[['chrom', 'start', 'end', 'rsid', 'gwas_category_count', 'gwas_category_lst']]

for count_pleio in range(1, 6):
    variant_pleio_i_bed_path = os.path.join(outdir_path, "variant_pleio_{}_flank_{}_hg38.bed".format(count_pleio, flank))
    if count_pleio == max_gwas_category_count:
        variant_pleio_i_bed_df = variant_bed_df.loc[variant_bed_df['gwas_category_count'] >= count_pleio, ]
    else:
        variant_pleio_i_bed_df = variant_bed_df.loc[variant_bed_df['gwas_category_count'] == count_pleio,]
    variant_pleio_i_bed_df = variant_pleio_i_bed_df.sort_values(by=['chrom', 'start', 'end'])
    variant_pleio_i_bed_df.to_csv(variant_pleio_i_bed_path, sep="\t", index=False, header=False)
