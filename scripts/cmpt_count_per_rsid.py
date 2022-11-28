import os
import pandas
import pathlib
import seaborn
import sys

from matplotlib import pyplot as plt


#%%
help_cmd_str = "todo"
try:
    annotated_tsv_path = sys.argv[1]
    max_gwas_class_count = int(sys.argv[2])
    count_per_rsid_gwas_tsv_path = sys.argv[3]
    count_per_rsid_gwas_egene_etissue_corr_png = sys.argv[4]
    if len(sys.argv) > 5:
        print("""Two many arguments!
        {}""".format(help_cmd_str))
        sys.exit(1)
except IndexError:
    print("""Argument missing!
    {}""".format(help_cmd_str))
    sys.exit(1)

outdir_path = os.path.dirname(count_per_rsid_gwas_tsv_path)
pathlib.Path(outdir_path).mkdir(parents=True, exist_ok=True)

#%%
coloc_df = pandas.read_csv(annotated_tsv_path, sep="\t")
coloc_df['chrom'] = coloc_df['chrom'].astype(int)

#%% etissue pleiotropy
etissue_df = coloc_df[['chrom', 'cytoband', 'pos', 'rsid', 'etissue_class']].drop_duplicates().groupby(['chrom', 'cytoband', 'pos', 'rsid']).agg({'etissue_class': ['size', (lambda x: (",".join(sorted(x))))]}).reset_index()
etissue_df.columns = ['chrom', 'cytoband', 'pos', 'rsid', 'etissue_label_count', 'etissue_class_lst']
etissue_df = etissue_df.sort_values(by=['etissue_label_count', 'chrom', 'pos', 'etissue_class_lst', 'rsid'], ascending=[False, True, True, True, True])
tsv_path = os.path.join(outdir_path, "count_per_rsid_etissue.tsv")
etissue_df.to_csv(tsv_path, sep="\t", index=False)

#%% egene pleiotropy
egene_df = coloc_df[['chrom', 'cytoband', 'pos', 'rsid', 'egene', 'egene_symbol']].drop_duplicates().groupby(['chrom', 'cytoband', 'pos', 'rsid']).agg({'egene': ['size', (lambda x: (",".join(sorted(x))))], 'egene_symbol': (lambda x: (",".join(sorted([str(i) for i in x]))))}).reset_index()
egene_df.columns = ['chrom', 'cytoband', 'pos', 'rsid', 'egene_count', 'egene_lst', 'egene_symbol_lst']
egene_df = egene_df.sort_values(by=['egene_count', 'chrom', 'pos', 'egene_lst', 'egene_symbol_lst', 'rsid'], ascending=[False, True, True, True, True, True])
tsv_path = os.path.join(outdir_path, "count_per_rsid_egene.tsv")
egene_df.to_csv(tsv_path, sep="\t", index=False)

#%% gwas pleiotropy
gwas_df = coloc_df[['chrom', 'cytoband', 'pos', 'rsid', 'gwas_class']].drop_duplicates().groupby(['chrom', 'cytoband', 'pos', 'rsid']).agg({'gwas_class': ['size', (lambda x: (",".join(sorted(x))))]}).reset_index()
gwas_df.columns = ['chrom', 'cytoband', 'pos', 'rsid', 'gwas_class_count', 'gwas_class_lst']
gwas_df = gwas_df.sort_values(by=['gwas_class_count', 'chrom', 'pos', 'gwas_class_lst', 'rsid'], ascending=[False, True, True, True, True])
tsv_path = os.path.join(outdir_path, "count_per_rsid_gwas.tsv")
gwas_df.to_csv(tsv_path, sep="\t", index=False)

#%% gwas variant pleiotropy for MS
gwas_ms_df = gwas_df.copy()
# gwas_ms_df['cytoband'] = gwas_ms_df['cytoband'].str.split('.', expand=True)[0]
gwas_ms_df = gwas_ms_df.drop_duplicates('cytoband', keep='first')
gwas_ms_df = gwas_ms_df.loc[gwas_ms_df['gwas_class_count'] >= 5]
# format output
gwas_ms_df.drop(['gwas_class_count'], inplace=True, axis=1)
gwas_ms_df['gwas_class_lst'] = gwas_ms_df['gwas_class_lst'].str.replace(',', ', ')
gwas_ms_df['pos']=gwas_ms_df['pos'].apply(lambda x : '{0:,}'.format(x))
tsv_path = os.path.join(outdir_path, "count_per_rsid_gwas_ms.tsv")
gwas_ms_df.to_csv(tsv_path, sep="\t", index=False)

#%% merge gwas, egene and etissue
m_df = pandas.merge(gwas_df, egene_df, on=['chrom', 'cytoband', 'pos', 'rsid'])
m_df = pandas.merge(m_df, etissue_df, on=['chrom', 'cytoband', 'pos', 'rsid'])
m_df = m_df[['chrom', 'cytoband', 'pos', 'rsid', 'gwas_class_count', 'gwas_class_lst', 'egene_count', 'egene_symbol_lst', 'etissue_label_count', 'etissue_class_lst', 'egene_lst']]
m_df.sort_values(['gwas_class_count', 'chrom', 'pos', 'rsid'], ascending=[False, True, True, True], inplace=True)
tsv_path = os.path.join(outdir_path, "count_per_rsid_gwas_egene_etissue.tsv")
m_df.to_csv(tsv_path, sep="\t", index=False)

m2df = m_df[['gwas_class_count', 'etissue_label_count', 'egene_count']]
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
variant_bed_df['start'] = variant_bed_df['pos'] - 1 - flank
variant_bed_df['end'] = variant_bed_df['pos'] + flank
variant_bed_df = variant_bed_df[['chrom', 'start', 'end', 'rsid', 'gwas_class_count', 'gwas_class_lst']]

for count_pleio in range(1, 6):
    variant_pleio_i_bed_path = os.path.join(outdir_path, "variant_pleio_{}_flank_{}_hg38.bed".format(count_pleio, flank))
    if count_pleio == max_gwas_class_count:
        variant_pleio_i_bed_df = variant_bed_df.loc[variant_bed_df['gwas_class_count'] >= count_pleio, ]
    else:
        variant_pleio_i_bed_df = variant_bed_df.loc[variant_bed_df['gwas_class_count'] == count_pleio,]
    variant_pleio_i_bed_df = variant_pleio_i_bed_df.sort_values(by=['chrom', 'start', 'end'])
    variant_pleio_i_bed_df.to_csv(variant_pleio_i_bed_path, sep="\t", index=False, header=False)
