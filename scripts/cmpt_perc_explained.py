from gwas2eqtl_pleiotropy.PathManager import PathManager

import os
import pandas
import pathlib
import sys


#%%
help_cmd_str = "todo"
try:
    url = sys.argv[1]
    supp_tabl_xlsx_path = sys.argv[2]
    if len(sys.argv) > 3:
        print("""Two many arguments!
        {}""".format(help_cmd_str))
        sys.exit(1)
except IndexError:
    print("""Argument missing!
    {}""".format(help_cmd_str))
    sys.exit(1)

# import pdb; pdb.set_trace()

# outdir_path = os.path.dirname(count_per_rsid_gwas_tsv_path)
# pathlib.Path(outdir_path).mkdir(parents=True, exist_ok=True)

# columns = ['chrom', 'pos38', 'cytoband', 'rsid', 'eqtl_beta', 'eqtl_gene_id', 'gwas_id', 'eqtl_id', 'etissue_class']
# coloc_df = pandas.read_sql(sql, con=url, columns=columns).drop_duplicates()

sql = 'select distinct chrom,pos,ea,gwas_id from tophits'
tophits_df = pandas.read_sql(sql, con=url)
tophits_df.rename({'ea': 'alt', 'pos': 'pos38'}, inplace=True, axis=1)

out_df = pandas.DataFrame()

for gwas_id in sorted(tophits_df['gwas_id'].unique()):
    print(gwas_id)
    sql = "select distinct chrom,pos38,alt,gwas_id from colocpleio where colocpleio.gwas_id='{}'".format(gwas_id)
    coloc_df = pandas.read_sql(sql, con=url)
    tophits_gwas_df = tophits_df.query("gwas_id=='{gwas_id}'".format(gwas_id=gwas_id))
    tophits_gwas_coloc_df = coloc_df.query("gwas_id=='{gwas_id}'".format(gwas_id=gwas_id))
    tophits_gwas_coloc_df = tophits_gwas_df.merge(tophits_gwas_coloc_df, on=['chrom', 'pos38', 'alt', 'gwas_id'])
    tophits_gwas_df = tophits_gwas_df.merge(tophits_gwas_coloc_df, on=['chrom', 'pos38', 'alt', 'gwas_id'], how='left', indicator=True)
    tophits_gwas_n = tophits_gwas_df.shape[0]
    tophits_explained_n = tophits_gwas_df.loc[tophits_gwas_df['_merge'] == 'both'].shape[0]
    out_df = pandas.concat([out_df, tophits_gwas_df], axis=0)

import pdb; pdb.set_trace()
