import os
import pandas
import pathlib
import sys
import seaborn as sns;
import sqlalchemy

from gwas2eqtl_pleiotropy.Logger import Logger

sns.set_theme(color_codes=True)


#%%
help_cmd_str = "todo"
try:
    snp_pp_h4 = float(sys.argv[1])
    sa_url = sys.argv[2]
    disease_corr_tsv_path = sys.argv[3]
    if len(sys.argv) > 4:
        print("""Two many arguments!
        {}""".format(help_cmd_str))
        sys.exit(1)
except IndexError:
    print("""Argument missing!
    {}""".format(help_cmd_str))
    sys.exit(1)

outdir_path = os.path.dirname(disease_corr_tsv_path)
pathlib.Path(outdir_path).mkdir(parents=True, exist_ok=True)

# sql = 'select * from colocpleio where snp_pp_h4>={}'.format(snp_pp_h4)
# columns = ['rsid', 'eqtl_beta', 'eqtl_gene_id', 'gwas_id', 'eqtl_id']
# d_df = pandas.read_sql(sql, con=sa_url, columns=columns).drop_duplicates()

#%%
columns = ['rsid', 'eqtl_beta', 'eqtl_gene_id', 'gwas_id', 'eqtl_id']
sql = 'select * from colocpleio where snp_pp_h4>={}'.format(snp_pp_h4)
engine = sqlalchemy.create_engine(sa_url)
with engine.begin() as conn:
    c_df = pandas.read_sql(sqlalchemy.text(sql), con=conn, columns=columns).drop_duplicates()

#%%
c_df = c_df.pivot_table(values='eqtl_beta', index=['rsid', 'eqtl_gene_id', 'eqtl_id'], columns='gwas_id', fill_value=0)

Logger.info("Spearman correlation")
c_df = c_df.corr(method='spearman')
c_df.to_csv(disease_corr_tsv_path, sep='\t', index=True, header=True)

