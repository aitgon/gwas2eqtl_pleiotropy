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
    annotated_tsv_path = sys.argv[1]
    disease_corr_tsv_path = sys.argv[2]
    if len(sys.argv) > 3:
        print("""Two many arguments!
        {}""".format(help_cmd_str))
        sys.exit(1)
except IndexError:
    print("""Argument missing!
    {}""".format(help_cmd_str))
    sys.exit(1)

outdir_path = os.path.dirname(disease_corr_tsv_path)
pathlib.Path(outdir_path).mkdir(parents=True, exist_ok=True)

#%%
Logger.info("Reading {}".format(annotated_tsv_path))
coloc_df = pandas.read_csv(annotated_tsv_path, sep="\t", usecols=['rsid', 'eqtl_beta', 'egene', 'eqtl_id', 'gwas_id'])

eqtl_id = 'GTEx_ge_blood'
# for eqtl_id in sorted(coloc_df.eqtl_id.unique()):
# d_df = coloc_df.loc[coloc_df.eqtl_id == eqtl_id, ['rsid', 'eqtl_beta', 'egene', 'gwas_id', 'eqtl_id']].drop_duplicates()
d_df = coloc_df[['rsid', 'eqtl_beta', 'egene', 'gwas_id', 'eqtl_id']].drop_duplicates()
d_df = d_df.pivot_table(values='eqtl_beta', index=['rsid', 'egene', 'eqtl_id'], columns='gwas_id', fill_value=0)

# d_df = d_df.loc[(d_df != 0).sum(axis=1) >= 5]
# import pdb; pdb.set_trace()

Logger.info("Spearman correlation")
d_df = d_df.corr(method='spearman')

d_df.fillna(0, inplace=True)

out_tsv_path = os.path.join(outdir_path, eqtl_id + '.tsv')
d_df.to_csv(disease_corr_tsv_path, sep='\t', index=True, header=True)

