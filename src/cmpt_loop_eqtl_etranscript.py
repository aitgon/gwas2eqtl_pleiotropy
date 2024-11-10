import os
import pandas
import pathlib
import seaborn
import sys
import sqlalchemy

from matplotlib import pyplot as plt
from statannotations.Annotator import Annotator

from constants import label_fontsize, tick_fontsize, seaborn_theme_dic, annotator_config_dic

#%%
help_cmd_str = "todo"
try:
    snp_pp_h4 = float(sys.argv[1])
    pleio_high_cutoff = int(sys.argv[2])
    db_url = sys.argv[3]
    loop_eqtl_etranscript_pleio_1_hg38_bed = sys.argv[4]
    if len(sys.argv) > 5:
        print("""Two many arguments!
        {}""".format(help_cmd_str))
        sys.exit(1)
except IndexError:
    print("""Argument missing!
    {}""".format(help_cmd_str))
    sys.exit(1)

seaborn.set_theme(**seaborn_theme_dic)
outdir_path = os.path.dirname(loop_eqtl_etranscript_pleio_1_hg38_bed)
pathlib.Path(outdir_path).mkdir(parents=True, exist_ok=True)

columns = ['chrom', 'pos38', 'rsid', 'ref', 'alt', 'gwas_category_ontology_term', 'eqtl_gene_id', 'eqtl_refseq_transcript_id', 'eqtl_refseq_transcript_start38']
cols_str = ','.join(columns)
sql = 'select distinct {} from colocpleio where snp_pp_h4>={}'.format(cols_str, snp_pp_h4)
engine = sqlalchemy.create_engine(db_url)
with engine.begin() as conn:
    coloc_df = pandas.read_sql(sqlalchemy.text(sql), con=conn)

#%% definition of variant for aggregation
variant_def_lst = ['chrom', 'pos38', 'rsid', 'ref', 'alt']

#%% per rsid, aggregate gwas categories
gwas_df = coloc_df[variant_def_lst + ['gwas_category_ontology_term']].drop_duplicates()
agg_dic = {'gwas_category_ontology_term': lambda x: x.tolist()}
gwas_df = gwas_df.groupby(variant_def_lst).agg(agg_dic).reset_index()
gwas_df['gwas_category_count'] = gwas_df['gwas_category_ontology_term'].apply(len)
gwas_df.rename({'gwas_category_ontology_term': 'gwas_category_lst'}, axis=1, inplace=True)

#%% per rsid, aggregate gwas categories
gwas_ontology_df = coloc_df[variant_def_lst + ['gwas_category_ontology_term']].drop_duplicates()
agg_dic = {'gwas_category_ontology_term': lambda x: x.tolist()}
gwas_ontology_df = gwas_ontology_df.groupby(variant_def_lst).agg(agg_dic).reset_index()
gwas_ontology_df['gwas_category_count'] = gwas_ontology_df['gwas_category_ontology_term'].apply(len)

#%% merge all aggregation columns
m_df = pandas.merge(coloc_df, gwas_df, on=variant_def_lst)

m_df = m_df[['chrom', 'pos38', 'eqtl_refseq_transcript_start38', 'gwas_category_count']]
bed_df = m_df[['chrom', 'pos38', 'eqtl_refseq_transcript_start38', 'gwas_category_count']].copy()
bed_df.loc[m_df['pos38'] > m_df['eqtl_refseq_transcript_start38'], 'eqtl_refseq_transcript_start38'] = m_df.loc[m_df['pos38'] > m_df['eqtl_refseq_transcript_start38'], 'pos38']
bed_df.loc[m_df['pos38'] > m_df['eqtl_refseq_transcript_start38'], 'pos38'] = m_df.loc[m_df['pos38'] > m_df['eqtl_refseq_transcript_start38'], 'eqtl_refseq_transcript_start38']
bed_df.rename({'pos38': 'start', 'eqtl_refseq_transcript_start38': 'end'}, axis=1, inplace=True)
bed_df = bed_df.dropna(how='any')
bed_df['chrom'] = 'chr' + bed_df['chrom'].astype(str)
#import pdb; pdb.set_trace()
bed_df['start']= bed_df['start'].astype(int)
bed_df['end'] = bed_df['end'].astype(int)
bed_df['name'] = bed_df['chrom'] + "_" + bed_df['start'].astype(str) + "_" + bed_df['end'].astype(str)
bed_df['start'] = bed_df['start'] - 1
bed_df = bed_df[['chrom', 'start', 'end', 'name', 'gwas_category_count']].drop_duplicates()

for gwas_category_count in sorted(bed_df['gwas_category_count'].unique()):
    pleio_bed_path = loop_eqtl_etranscript_pleio_1_hg38_bed.replace('pleio_1', 'pleio_{}'.format(gwas_category_count))
    if gwas_category_count == pleio_high_cutoff:
        bed_pleio_df = bed_df.loc[bed_df['gwas_category_count'] >= pleio_high_cutoff,]
    else:
        bed_pleio_df = bed_df.loc[bed_df['gwas_category_count'] == gwas_category_count, ]
    bed_pleio_df = bed_pleio_df.sort_values(by=['chrom', 'start', 'end'])
    bed_pleio_df.to_csv(pleio_bed_path, sep="\t", index=False, header=False)
    if gwas_category_count == pleio_high_cutoff:
        break

