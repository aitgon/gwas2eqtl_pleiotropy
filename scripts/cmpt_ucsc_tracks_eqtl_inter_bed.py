import os
import pandas
import pathlib
import sys

import sqlalchemy


#%%
help_cmd_str = "todo"
try:
    snp_pp_h4 = float(sys.argv[1])
    sa_url = sys.argv[2]
    outdir_path = sys.argv[3]
    if len(sys.argv) > 4:
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
engine = sqlalchemy.create_engine(sa_url)
with engine.begin() as conn:
    df = pandas.read_sql(sqlalchemy.text(sql), con=conn).drop_duplicates()

df = df.loc[~df['eqtl_refseq_transcript_start38'].isna()]
df = df[['chrom', 'pos38', 'rsid', 'eqtl_gene_symbol', 'eqtl_beta', 'eqtl_id',
         'eqtl_refseq_transcript_start38', 'eqtl_refseq_transcript_end38', 'eqtl_refseq_transcript_strand',
         'etissue_category_term', 'gwas_trait', 'gwas_category_ontology_term', 'gwas_id']].drop_duplicates()

# select longest transcript
df['transcript_length'] = df['eqtl_refseq_transcript_end38'] - df['eqtl_refseq_transcript_start38']
df = df.sort_values(['transcript_length'], ascending=[False])
df = df.drop_duplicates(['chrom', 'pos38', 'rsid', 'eqtl_gene_symbol', 'eqtl_beta', 'eqtl_id', 'eqtl_refseq_transcript_strand'], keep='first')

df2 = df.copy()
df2["sourceChrom"] = "chr" + df2['chrom'].astype(str)
df2["sourceStart"] = df2['pos38'].astype(int) - 1
df2["sourceEnd"] = df2['pos38'].astype(int)
df2["sourceName"] = 'rs' + df2['rsid'].astype(str) + '/' + df2['gwas_category_ontology_term'].str.replace(" ", "")
df2["sourceStrand"] = '.'

df2["targetChrom"] = "chr" + df2['chrom'].astype(str)
df2["targetStart"] = df2['eqtl_refseq_transcript_start38'].astype(int) - 1
df2["targetEnd"] = df2['eqtl_refseq_transcript_end38'].astype(int)
df2["targetName"] = df2['eqtl_gene_symbol'].astype(str)
df2["targetStrand"] = df2['eqtl_refseq_transcript_strand']

df2["#chrom"] = "chr" + df2['chrom'].astype(str)
df2['chromStart'] = df2[['sourceStart', 'sourceEnd', 'targetStart', 'targetEnd']].min(axis=1)
df2['chromEnd'] = df2[['sourceStart', 'sourceEnd', 'targetStart', 'targetEnd']].max(axis=1)

df2['score'] = 1000

df2['exp'] = df['eqtl_id'] + '/' + df['gwas_trait'].str.replace(' ', '_') + '/' + df['gwas_id']
df2['value'] = df['eqtl_beta'].abs()
df2['color'] = '#008000'
df2.loc[df['eqtl_beta'] < 0, 'color'] = '#FF0000'
df2['name'] = df2['sourceName'] + '/' + df2['targetName'] + '/' + df2['exp']

for eqtl_id in sorted(df2['eqtl_id'].unique()):
    print(eqtl_id)
    cell_df = df2.loc[df2['eqtl_id'] == eqtl_id].copy()
    cell_df.drop(df.columns.tolist(), axis=1, inplace=True)
    cell_df['etissue_category_term'] = eqtl_id

#     track_config_str = """browser position chr5:132239646-132497907
# browser pack GWAS/eQTL colocalization
# track type=interact name="coloc_gwas_eqtl_{eqtl_id}" description="An interact file" interactDirectional=true maxHeightPixels=200:100:50 visibility=full
# """.format(eqtl_id=eqtl_id)
    cell_track_tsv = os.path.join(outdir_path, "{}.inter.bed".format(eqtl_id.replace(' ', '_')))
    # with open(cell_track_tsv, 'w') as fout:
    #     fout.write(track_config_str)

    cell_df = cell_df[['#chrom', 'chromStart', 'chromEnd', 'name', 'score', 'value', 'exp',
               'color', 'sourceChrom', 'sourceStart', 'sourceEnd', 'sourceName', 'sourceStrand',
               'targetChrom', 'targetStart', 'targetEnd', 'targetName', 'targetStrand']]
    cell_df.sort_values(by=cell_df.columns.tolist(), inplace=True)
    # import pdb; pdb.set_trace()
    cell_df.to_csv(cell_track_tsv, sep='\t', header=True, index=False)
