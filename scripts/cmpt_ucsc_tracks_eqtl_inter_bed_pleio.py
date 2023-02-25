import os
import pandas
import pathlib
import sys
import sqlalchemy

#%%
help_cmd_str = "todo"
try:
    snp_pp_h4 = float(sys.argv[1])
    gwas_category_count = int(sys.argv[2])
    db_url = sys.argv[3]
    count_per_rsid_gwas_egene_etissue_ods = sys.argv[4]
    outdir_path = sys.argv[5]
    if len(sys.argv) > 6:
        print("""Two many arguments!
        {}""".format(help_cmd_str))
        sys.exit(1)
except IndexError:
    print("""Argument missing!
    {}""".format(help_cmd_str))
    sys.exit(1)

pathlib.Path(outdir_path).mkdir(parents=True, exist_ok=True)

#%%
sql = 'select * from colocpleio where snp_pp_h4>={}'.format(snp_pp_h4)
engine = sqlalchemy.create_engine(db_url)
with engine.begin() as conn:
    df = pandas.read_sql(sqlalchemy.text(sql), con=conn).drop_duplicates()

#%%
usecols = ['chrom', 'pos19', 'pos38', 'rsid', 'ref', 'alt', 'gwas_category_count']
pleio_df = pandas.read_excel(count_per_rsid_gwas_egene_etissue_ods, engine='odf', usecols=usecols)
df = df.merge(pleio_df, on=['chrom', 'pos19', 'pos38', 'rsid', 'ref', 'alt'])
df = df.loc[df['gwas_category_count'] == gwas_category_count]

#%%
df = df.loc[~df['eqtl_refseq_transcript_start38'].isna()]
df = df[['chrom', 'pos38', 'rsid', 'eqtl_gene_symbol', 'eqtl_beta', 'eqtl_id',
 'eqtl_refseq_transcript_start38', 'eqtl_refseq_transcript_end38', 'eqtl_refseq_transcript_strand',
 'etissue_category_term', 'gwas_trait', 'gwas_category_ontology_term', 'gwas_id', 'snp_pp_h4', 'gwas_category_count']].drop_duplicates()

# select longest transcript
df['transcript_length'] = df['eqtl_refseq_transcript_end38'] - df['eqtl_refseq_transcript_start38']
df = df.sort_values(['transcript_length'], ascending=[False])
df = df.drop_duplicates(['chrom', 'pos38', 'rsid', 'eqtl_gene_symbol', 'eqtl_beta', 'eqtl_id', 'eqtl_refseq_transcript_strand'], keep='first')

df2 = df.copy()
df2["sourceChrom"] = "chr" + df2['chrom'].astype(str)
df2["sourceStart"] = df2['pos38'].astype(int) - 1
df2["sourceEnd"] = df2['pos38'].astype(int)
df2["sourceName"] = 'rs' + df2['rsid'].astype(str) + '/' + df2['gwas_trait'].str.replace(" ", "_")
df2["sourceStrand"] = '.'

df2["targetChrom"] = "chr" + df2['chrom'].astype(str)
df2["targetStart"] = df2['eqtl_refseq_transcript_start38'].astype(int) - 1
df2["targetEnd"] = df2['eqtl_refseq_transcript_end38'].astype(int)
df2["targetName"] = df2['eqtl_gene_symbol'].astype(str)
df2["targetStrand"] = df2['eqtl_refseq_transcript_strand']

df2["#chrom"] = "chr" + df2['chrom'].astype(str)
df2['chromStart'] = df2[['sourceStart', 'sourceEnd', 'targetStart', 'targetEnd']].min(axis=1)
df2['chromEnd'] = df2[['sourceStart', 'sourceEnd', 'targetStart', 'targetEnd']].max(axis=1)

df2['score'] = (df2['snp_pp_h4'] * 1000).astype(int)

df2['exp'] = df2['eqtl_id'] + '/' + df2['gwas_id']
df2['value'] = df2['eqtl_beta']
df2['color'] = '#008000'
df2.loc[df['eqtl_beta'] < 0, 'color'] = '#FF0000'
df2['name'] = df2['sourceName'] + '/' + df2['targetName'] + '/' + df2['exp']

eqtl_tsv_path = 'https://raw.githubusercontent.com/eQTL-Catalogue/eQTL-Catalogue-resources/master/tabix/tabix_ftp_paths.tsv'
eqtl_df = pandas.read_csv(eqtl_tsv_path, sep="\t", usecols=[0, 6, 8])
eqtl_df = eqtl_df.loc[eqtl_df['ftp_path'].str.contains('/ge/|/microarray/',regex=True,na=False), ]
# keep urls containing "ge" or "microarray"
eqtl_df['index'] = (eqtl_df['ftp_path'].str.replace('.all.tsv.gz', '', regex=True)).str.split('/', expand=True)[10].tolist()
eqtl_df.set_index('index', drop=True, verify_integrity=True, inplace=True)
eqtl_id_lst = (eqtl_df.index).tolist()

for eqtl_id in sorted(eqtl_id_lst):
    # print(eqtl_id)
    cell_df = df2.loc[df2['eqtl_id'] == eqtl_id].copy()
    cell_df.drop(df.columns.tolist(), axis=1, inplace=True)
    cell_df['etissue_category_term'] = eqtl_id

    cell_track_tsv = os.path.join(outdir_path, "{}.inter.bed".format(eqtl_id.replace(' ', '_')))

    cell_df = cell_df[['#chrom', 'chromStart', 'chromEnd', 'name', 'score', 'value', 'exp',
               'color', 'sourceChrom', 'sourceStart', 'sourceEnd', 'sourceName', 'sourceStrand',
               'targetChrom', 'targetStart', 'targetEnd', 'targetName', 'targetStrand']]
    cell_df.sort_values(by=cell_df.columns.tolist(), inplace=True)
    cell_df.to_csv(cell_track_tsv, sep='\t', header=True, index=False)
