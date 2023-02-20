import os
import pandas
import pathlib
import seaborn
import sys

import sqlalchemy
from matplotlib import pyplot as plt


#%%
help_cmd_str = "todo"
try:
    snp_pp_h4 = float(sys.argv[1])
    sa_url = sys.argv[2]
    track_tsv = sys.argv[3]
    if len(sys.argv) > 4:
        print("""Two many arguments!
        {}""".format(help_cmd_str))
        sys.exit(1)
except IndexError:
    print("""Argument missing!
    {}""".format(help_cmd_str))
    sys.exit(1)

outdir_path = os.path.dirname(track_tsv)
pathlib.Path(outdir_path).mkdir(parents=True, exist_ok=True)

#%%
sql = 'select * from colocpleio where snp_pp_h4>={}'.format(snp_pp_h4)
engine = sqlalchemy.create_engine(sa_url)
with engine.begin() as conn:
    df = pandas.read_sql(sqlalchemy.text(sql), con=conn).drop_duplicates()

df = df.loc[~df['eqtl_refseq_transcript_start38'].isna()]

df = df[['chrom', 'pos38', 'rsid', 'eqtl_gene_symbol', 'eqtl_beta', 'eqtl_id',
         'eqtl_refseq_transcript_start38', 'eqtl_refseq_transcript_end38', 'eqtl_refseq_transcript_strand',
         'etissue_category_term', 'gwas_trait', 'gwas_id']].drop_duplicates()

# select longest transcript
df['transcript_length'] = df['eqtl_refseq_transcript_end38'] - df['eqtl_refseq_transcript_start38']
df = df.sort_values(['transcript_length'], ascending=[False])
df = df.drop_duplicates(['chrom', 'pos38', 'rsid', 'eqtl_gene_symbol', 'eqtl_beta', 'eqtl_id', 'eqtl_refseq_transcript_strand'], keep='first')

df2 = df.copy()
df2["sourceChrom"] = "chr" + df2['chrom'].astype(str)
df2["sourceStart"] = df2['pos38'].astype(int) - 1
df2["sourceEnd"] = df2['pos38'].astype(int)
df2["sourceName"] = 'rs' + df2['rsid'].astype(str)
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

df2['exp'] = df['eqtl_id'] + '/' + df['gwas_trait'] + '/' + df['gwas_id']
df2['value'] = df['eqtl_beta'].abs()
df2['color'] = '#008000'
df2.loc[df['eqtl_beta'] < 0, 'color'] = '#FF0000'
df2['name'] = df2['sourceName'] + '/' + df2['targetName'] + '/' + df2['exp']

for etissue_category_term in sorted(df2['etissue_category_term'].unique()):
    etissue_df = df2.loc[df2['etissue_category_term'] == etissue_category_term]
    etissue_df.drop(df.columns.tolist(), axis=1, inplace=True)
    etissue_df['etissue_category_term'] = etissue_category_term

    track_config_str = """browser position chr5:132239646-132497907
    browser pack Pleiotropic_eQTLs
    track type=interact name="gwas2eqtl_{etissue_category_term}" description="An interact file" interactDirectional=true maxHeightPixels=200:100:50 visibility=full
    """.format(etissue_category_term=etissue_category_term)
    with open(track_tsv, 'w') as fout:
        fout.write(track_config_str)

    df2 = df2[['#chrom', 'chromStart', 'chromEnd', 'name', 'score', 'value', 'exp',
               'color', 'sourceChrom', 'sourceStart', 'sourceEnd', 'sourceName', 'sourceStrand',
               'targetChrom', 'targetStart', 'targetEnd', 'targetName', 'targetStrand']]
    df2.to_csv(track_tsv, sep='\t', header=True, mode='a', index=False)
