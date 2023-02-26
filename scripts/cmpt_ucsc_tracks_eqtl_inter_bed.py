import os

import numpy
import pandas
import pathlib
import sys

import sqlalchemy


#%%
help_cmd_str = "todo"
try:
    snp_pp_h4 = float(sys.argv[1])
    db_url = sys.argv[2]
    bed_outdir_path = sys.argv[3]
    if len(sys.argv) > 4:
        print("""Two many arguments!
        {}""".format(help_cmd_str))
        sys.exit(1)
except IndexError:
    print("""Argument missing!
    {}""".format(help_cmd_str))
    sys.exit(1)

pathlib.Path(bed_outdir_path).mkdir(parents=True, exist_ok=True)

#%%
sql = 'select * from colocpleio where snp_pp_h4>={}'.format(snp_pp_h4)
engine = sqlalchemy.create_engine(db_url)
with engine.begin() as conn:
    df = pandas.read_sql(sqlalchemy.text(sql), con=conn).drop_duplicates()

df = df.loc[~df['eqtl_refseq_transcript_start38'].isna()]
df = df[['chrom', 'pos38', 'rsid', 'eqtl_gene_symbol', 'eqtl_beta', 'eqtl_id',
 'eqtl_refseq_transcript_start38', 'eqtl_refseq_transcript_end38', 'eqtl_refseq_transcript_strand',
 'etissue_category_term', 'gwas_trait', 'gwas_category_ontology_term', 'gwas_beta', 'gwas_id', 'snp_pp_h4']].drop_duplicates()

# select longest transcript
df['transcript_length'] = df['eqtl_refseq_transcript_end38'] - df['eqtl_refseq_transcript_start38']
df = df.sort_values(['transcript_length'], ascending=[False])
df = df.drop_duplicates(['chrom', 'pos38', 'rsid', 'eqtl_gene_symbol', 'eqtl_beta', 'eqtl_id', 'eqtl_refseq_transcript_strand'], keep='first')

df2 = df.copy()
df2["sourceChrom"] = "chr" + df2['chrom'].astype(str)
df2["sourceStart"] = df2['pos38'].astype(int) - 1
df2["sourceEnd"] = df2['pos38'].astype(int)
df2["sourceName"] = 'rs' + df2['rsid'].astype(str) + '/' + df2['gwas_trait'].str.replace(" ", "_") + '/' + df2['gwas_beta'].astype(str)
df2["sourceStrand"] = '.'

df2["targetChrom"] = "chr" + df2['chrom'].astype(str)
df2["targetStart"] = df2['eqtl_refseq_transcript_start38'].astype(int) - 1
df2["targetEnd"] = df2['eqtl_refseq_transcript_end38'].astype(int)
df2["targetName"] = df2['eqtl_gene_symbol'].astype(str) + '/' + df2['eqtl_beta'].astype(str)
df2["targetStrand"] = df2['eqtl_refseq_transcript_strand']

df2["#chrom"] = "chr" + df2['chrom'].astype(str)
df2['chromStart'] = df2[['sourceStart', 'sourceEnd', 'targetStart', 'targetEnd']].min(axis=1)
df2['chromEnd'] = df2[['sourceStart', 'sourceEnd', 'targetStart', 'targetEnd']].max(axis=1)

df2['score'] = (df2['snp_pp_h4'] * 1000).astype(int)

df2['exp'] = df['eqtl_id'] + '/' + df['gwas_id']
df2['value'] = df['eqtl_beta']

df2['color'] = '#FF0000'
df2.loc[df['eqtl_beta'] < 0, 'color'] = '#0000FF'
df2['name'] = df2['sourceName'] + '/' + df2['targetName'] + '/' + df2['exp']

df2['beta_equal'] = numpy.sign(df2['gwas_beta'])*numpy.sign(df2['eqtl_beta'])

for eqtl_id in sorted(df2['eqtl_id'].unique()):

    cell_df = df2.loc[df2['eqtl_id'] == eqtl_id].copy()
    cell_df.drop(df.columns.tolist(), axis=1, inplace=True)


    cell_beta_equal_mask = cell_df['beta_equal'] == 1
    cell_df = cell_df[['#chrom', 'chromStart', 'chromEnd', 'name', 'score', 'value', 'exp',
               'color', 'sourceChrom', 'sourceStart', 'sourceEnd', 'sourceName', 'sourceStrand',
               'targetChrom', 'targetStart', 'targetEnd', 'targetName', 'targetStrand']]

    cell_beta_equal_df = cell_df.loc[cell_beta_equal_mask].copy()
    cell_beta_equal_df.sort_values(by=cell_df.columns.tolist(), inplace=True)
    cell_track_beta_equal_tsv = os.path.join(bed_outdir_path, 'beta_equal', "{}.inter.bed".format(eqtl_id.replace(' ', '_')))
    pathlib.Path(os.path.dirname(cell_track_beta_equal_tsv)).mkdir(parents=True, exist_ok=True)
    cell_beta_equal_df.to_csv(cell_track_beta_equal_tsv, sep='\t', header=True, index=False)

    cell_beta_unequal_df = cell_df.loc[~cell_beta_equal_mask].copy()
    cell_beta_unequal_df.sort_values(by=cell_df.columns.tolist(), inplace=True)
    cell_track_beta_unequal_tsv = os.path.join(bed_outdir_path, 'beta_unequal', "{}.inter.bed".format(eqtl_id.replace(' ', '_')))
    pathlib.Path(os.path.dirname(cell_track_beta_unequal_tsv)).mkdir(parents=True, exist_ok=True)
    cell_beta_unequal_df.to_csv(cell_track_beta_unequal_tsv, sep='\t', header=True, index=False)
