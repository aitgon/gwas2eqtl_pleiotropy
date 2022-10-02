import os
import pandas
import pathlib
import sys


#%%
help_cmd_str = "todo"
try:
    annotated_tsv_path = sys.argv[1]
    upper_var_gwas_cat_count = int(sys.argv[2])
    count_per_rsid_gwas_tsv_path = sys.argv[3]
    if len(sys.argv) > 4:
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
etissue_df = coloc_df[['chrom', 'cytoband', 'pos', 'rsid', 'etissue_category']].drop_duplicates().groupby(['chrom', 'cytoband', 'pos', 'rsid']).agg({'etissue_category': ['size', (lambda x: (",".join(sorted(x))))]}).reset_index()
etissue_df.columns = ['chrom', 'cytoband', 'pos', 'rsid', 'etissue_label_count', 'etissue_subcategory_lst']
etissue_df = etissue_df.sort_values(by=['etissue_label_count', 'chrom', 'pos', 'etissue_subcategory_lst', 'rsid'], ascending=[False, True, True, True, True])
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
gwas_df.columns = ['chrom', 'cytoband', 'pos', 'rsid', 'gwas_category_count', 'gwas_category_lst']
gwas_df = gwas_df.sort_values(by=['gwas_category_count', 'chrom', 'pos', 'gwas_category_lst', 'rsid'], ascending=[False, True, True, True, True])
tsv_path = os.path.join(outdir_path, "count_per_rsid_gwas.tsv")
gwas_df.to_csv(tsv_path, sep="\t", index=False)

#%% gwas variant pleiotropy for MS
gwas_ms_df = gwas_df.copy()
# gwas_ms_df['cytoband'] = gwas_ms_df['cytoband'].str.split('.', expand=True)[0]
gwas_ms_df = gwas_ms_df.drop_duplicates('cytoband', keep='first')
gwas_ms_df = gwas_ms_df.loc[gwas_ms_df['gwas_category_count'] >= 5]
# format output
gwas_ms_df.drop(['gwas_category_count'], inplace=True, axis=1)
gwas_ms_df['gwas_category_lst'] = gwas_ms_df['gwas_category_lst'].str.replace(',', ', ')
gwas_ms_df['pos']=gwas_ms_df['pos'].apply(lambda x : '{0:,}'.format(x))
tsv_path = os.path.join(outdir_path, "count_per_rsid_gwas_ms.tsv")
gwas_ms_df.to_csv(tsv_path, sep="\t", index=False)

#%% merge gwas, egene and etissue
m_df = pandas.merge(gwas_df, egene_df, on=['chrom', 'cytoband', 'pos', 'rsid'])
m_df = pandas.merge(m_df, etissue_df, on=['chrom', 'cytoband', 'pos', 'rsid'])
m_df = m_df[['chrom', 'cytoband', 'pos', 'rsid', 'gwas_category_count', 'gwas_category_lst', 'egene_count', 'egene_symbol_lst', 'etissue_label_count', 'etissue_subcategory_lst', 'egene_lst']]
m_df.sort_values(['gwas_category_count', 'chrom', 'pos', 'rsid'], ascending=[False, True, True, True], inplace=True)
tsv_path = os.path.join(outdir_path, "count_per_rsid_gwas_egene_etissue.tsv")
m_df.to_csv(tsv_path, sep="\t", index=False)

#%########################################### bed files, flanking=0
# bed files of variants splitted by gwas categories
flank = 10
variant_bed_df = gwas_df.copy()
variant_bed_df['chrom'] = 'chr' + variant_bed_df['chrom'].astype(str)
variant_bed_df['start'] = variant_bed_df['pos'] - 1 - flank
variant_bed_df['end'] = variant_bed_df['pos'] + flank
variant_bed_df = variant_bed_df[['chrom', 'start', 'end', 'rsid', 'gwas_category_count', 'gwas_category_lst']]

for count_pleio in range(1, 6):
    variant_pleio_i_bed_path = os.path.join(outdir_path, "variant_pleio_{}_flank_{}_hg38.bed".format(count_pleio, flank))
    if count_pleio == upper_var_gwas_cat_count:
        variant_pleio_i_bed_df = variant_bed_df.loc[variant_bed_df['gwas_category_count'] >= count_pleio, ]
    else:
        variant_pleio_i_bed_df = variant_bed_df.loc[variant_bed_df['gwas_category_count'] == count_pleio,]
    variant_pleio_i_bed_df = variant_pleio_i_bed_df.sort_values(by=['chrom', 'start', 'end'])
    variant_pleio_i_bed_df.to_csv(variant_pleio_i_bed_path, sep="\t", index=False, header=False)
