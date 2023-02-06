import os
import pathlib

import pandas
import sys

from gwas2eqtl_pleiotropy.Logger import Logger

#%%
help_cmd_str = "todo"
try:
    snp_pp_h4 = float(sys.argv[1])
    url = sys.argv[2]
    loci_explained_perc_tsv = sys.argv[3]
    if len(sys.argv) > 4:
        print("""Two many arguments!
        {}""".format(help_cmd_str))
        sys.exit(1)
except IndexError:
    print("""Argument missing!
    {}""".format(help_cmd_str))
    sys.exit(1)

pathlib.Path(os.path.dirname(loci_explained_perc_tsv)).mkdir(exist_ok=True, parents=True)

sql = 'select distinct chrom, pos, ea, tophits.gwas_id, trait as gwas_trait, op.consortium, op.pmid from tophits, open_gwas_info op where op.gwas_id=tophits.gwas_id'
tophits_df = pandas.read_sql(sql, con=url)
tophits_df.rename({'ea': 'alt', 'pos': 'pos38'}, inplace=True, axis=1)

concat_df = pandas.DataFrame()

for gwas_id in sorted(tophits_df['gwas_id'].unique()):
    Logger.info(gwas_id)
    sql = "select distinct chrom, pos38, alt, gwas_id, gwas_trait from colocpleio where colocpleio.gwas_id='{}' and snp_pp_h4>={}".format(gwas_id, snp_pp_h4)
    coloc_df = pandas.read_sql(sql, con=url)
    tophits_gwas_df = tophits_df.query("gwas_id=='{gwas_id}'".format(gwas_id=gwas_id))
    tophits_gwas_coloc_df = coloc_df.query("gwas_id=='{gwas_id}'".format(gwas_id=gwas_id))
    merge_col_lst = ['chrom', 'pos38', 'alt', 'gwas_id', 'gwas_trait']
    tophits_gwas_coloc_df = (tophits_gwas_df[merge_col_lst]).merge(tophits_gwas_coloc_df[merge_col_lst], on=merge_col_lst)
    tophits_gwas_df = tophits_gwas_df.merge(tophits_gwas_coloc_df, on=['chrom', 'pos38', 'alt', 'gwas_id', 'gwas_trait'], how='left', indicator=True)
    concat_df = pandas.concat([concat_df, tophits_gwas_df], axis=0)

# Percentage explained all
out_df = concat_df.groupby(['gwas_id', 'gwas_trait', '_merge']).size().reset_index()
out_df = out_df.loc[out_df['_merge'] != 'right_only']
out_df = tophits_df[['gwas_id', 'gwas_trait']].merge(out_df, on=['gwas_id', 'gwas_trait']).drop_duplicates()
out_df = out_df.pivot_table(index=['gwas_id', 'gwas_trait'], columns=['_merge'], values=0)

out_df['loci_count'] = out_df.apply(sum, axis=1)
out_df.rename({'both': 'loci_coloc_count', 'left_only': 'loci_noncoloc_count'}, axis=1, inplace=True)
out_df['loci_coloc_perc'] = out_df['loci_coloc_count'] / out_df['loci_count'] * 100
out_df['loci_coloc_perc'] = out_df['loci_coloc_perc'].apply(int)

# Percentage explained non mch
concat_nomhc_df = concat_df.query('not (chrom==6 & pos38>=25000000 & pos38<=35000000)')
out_nomhc_df = concat_nomhc_df.groupby(['gwas_id', 'gwas_trait', '_merge']).size().reset_index()
out_nomhc_df = out_nomhc_df.loc[out_nomhc_df['_merge'] != 'right_only']
out_nomhc_df = tophits_df[['gwas_id', 'gwas_trait']].merge(out_nomhc_df, on=['gwas_id', 'gwas_trait']).drop_duplicates()
out_nomhc_df = out_nomhc_df.pivot_table(index=['gwas_id', 'gwas_trait'], columns=['_merge'], values=0)

out_nomhc_df['loci_nomhc_count'] = out_nomhc_df.apply(sum, axis=1)
out_nomhc_df.rename({'both': 'loci_nomhc_coloc_count', 'left_only': 'loci_nomhc_noncoloc_count'}, axis=1, inplace=True)
out_nomhc_df['loci_nomhc_coloc_perc'] = out_nomhc_df['loci_nomhc_coloc_count'] / out_nomhc_df['loci_nomhc_count'] * 100
out_nomhc_df['loci_nomhc_coloc_perc'] = out_nomhc_df['loci_nomhc_coloc_perc'].apply(int)
out_df = (out_df.merge(out_nomhc_df[['loci_nomhc_count', 'loci_nomhc_coloc_count', 'loci_nomhc_coloc_perc']], left_index=True, right_index=True))

out_df = out_df.merge(tophits_df[['gwas_id', 'gwas_trait', 'consortium', 'pmid']].drop_duplicates(), left_index=True, right_on=['gwas_id', 'gwas_trait'])
out_df.set_index(['gwas_id', 'gwas_trait'], verify_integrity=True, inplace=True)
out_df = out_df[['consortium', 'pmid', 'loci_count', 'loci_coloc_count', 'loci_coloc_perc', 'loci_nomhc_count', 'loci_nomhc_coloc_count', 'loci_nomhc_coloc_perc']]
out_df.to_csv(loci_explained_perc_tsv, sep='\t', index=True)
