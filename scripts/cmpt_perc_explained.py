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
    perc_explained_tsv_path = sys.argv[3]
    if len(sys.argv) > 4:
        print("""Two many arguments!
        {}""".format(help_cmd_str))
        sys.exit(1)
except IndexError:
    print("""Argument missing!
    {}""".format(help_cmd_str))
    sys.exit(1)

pathlib.Path(os.path.dirname(perc_explained_tsv_path)).mkdir(exist_ok=True, parents=True)

sql = 'select distinct chrom,pos,ea, tophits.gwas_id,gwas_trait from tophits,gwas_annot where gwas_annot.gwas_id=tophits.gwas_id'
tophits_df = pandas.read_sql(sql, con=url)
tophits_df.rename({'ea': 'alt', 'pos': 'pos38'}, inplace=True, axis=1)

concat_df = pandas.DataFrame()

for gwas_id in sorted(tophits_df['gwas_id'].unique()):
    Logger.info(gwas_id)
    sql = "select distinct chrom,pos38,alt,gwas_id,gwas_trait from colocpleio where colocpleio.gwas_id='{}' and snp_pp_h4>={}".format(gwas_id, snp_pp_h4)
    coloc_df = pandas.read_sql(sql, con=url)
    tophits_gwas_df = tophits_df.query("gwas_id=='{gwas_id}'".format(gwas_id=gwas_id))
    tophits_gwas_coloc_df = coloc_df.query("gwas_id=='{gwas_id}'".format(gwas_id=gwas_id))
    tophits_gwas_coloc_df = tophits_gwas_df.merge(tophits_gwas_coloc_df, on=['chrom', 'pos38', 'alt', 'gwas_id', 'gwas_trait'])
    tophits_gwas_df = tophits_gwas_df.merge(tophits_gwas_coloc_df, on=['chrom', 'pos38', 'alt', 'gwas_id', 'gwas_trait'], how='left', indicator=True)
    concat_df = pandas.concat([concat_df, tophits_gwas_df], axis=0)

# Percentage explained all
out_df = concat_df.groupby(['gwas_id', 'gwas_trait', '_merge']).size().reset_index()
out_df = out_df.loc[out_df['_merge'] != 'right_only']
out_df = tophits_df[['gwas_id', 'gwas_trait']].merge(out_df, on=['gwas_id', 'gwas_trait']).drop_duplicates()
out_df = out_df.pivot_table(index=['gwas_id', 'gwas_trait'], columns=['_merge'], values=0)

out_df['n'] = out_df.apply(sum, axis=1)
out_df.rename({'both': 'n_explained', 'left_only': 'n_non_explained'}, axis=1, inplace=True)
out_df['perc_explained_all'] = out_df['n_explained'] / out_df['n'] * 100
out_df['perc_explained_all'] = out_df['perc_explained_all'].apply(int)

# Percentage explained non mch
concat_nomhc_df = concat_df.query('not (chrom==6 & pos38>=25000000 & pos38<=35000000)')
out_nomhc_df = concat_nomhc_df.groupby(['gwas_id', 'gwas_trait', '_merge']).size().reset_index()
out_nomhc_df = out_nomhc_df.loc[out_nomhc_df['_merge'] != 'right_only']
out_nomhc_df = tophits_df[['gwas_id', 'gwas_trait']].merge(out_nomhc_df, on=['gwas_id', 'gwas_trait']).drop_duplicates()
out_nomhc_df = out_nomhc_df.pivot_table(index=['gwas_id', 'gwas_trait'], columns=['_merge'], values=0)

out_nomhc_df['n'] = out_nomhc_df.apply(sum, axis=1)
out_nomhc_df.rename({'both': 'n_explained', 'left_only': 'n_non_explained'}, axis=1, inplace=True)
out_nomhc_df['perc_explained_nomhc'] = out_nomhc_df['n_explained'] / out_nomhc_df['n'] * 100
out_nomhc_df['perc_explained_nomhc'] = out_nomhc_df['perc_explained_nomhc'].apply(int)
out_df = (out_df.merge(out_nomhc_df[['perc_explained_nomhc']], left_index=True, right_index=True))
out_df.to_csv(perc_explained_tsv_path, sep='\t', index=True)
