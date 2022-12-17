import os
import pandas
import pathlib
import sys
import seaborn as sns;

from gwas2eqtl_pleiotropy.Logger import Logger

sns.set_theme(color_codes=True)


#%%
help_cmd_str = "todo"
try:
    snp_pp_h4 = float(sys.argv[1])
    url = sys.argv[2]
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

sql = 'select * from colocpleio where snp_pp_h4>={}'.format(snp_pp_h4)
columns = ['rsid', 'eqtl_beta', 'eqtl_gene_id', 'gwas_id', 'eqtl_id']
d_df = pandas.read_sql(sql, con=url, columns=columns).drop_duplicates()

#%%
# Logger.info("Reading {}".format(annotated_tsv_path))
# coloc_df = pandas.read_csv(annotated_tsv_path, sep="\t", usecols=['rsid', 'eqtl_beta', 'egene', 'eqtl_id', 'gwas_id'])

# d_df = coloc_df[['rsid', 'eqtl_beta', 'egene', 'gwas_id', 'eqtl_id']].drop_duplicates()
d_df = d_df.pivot_table(values='eqtl_beta', index=['rsid', 'eqtl_gene_id', 'eqtl_id'], columns='gwas_id', fill_value=0)

Logger.info("Spearman correlation")
d_df = d_df.corr(method='spearman')
# import pdb; pdb.set_trace()
d_df.to_csv(disease_corr_tsv_path, sep='\t', index=True, header=True)

