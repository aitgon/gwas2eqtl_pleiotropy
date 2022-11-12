"""In this script, we compute the percentage of GWAS annotate with an eQTL (H4>=0.8)"""
import pdb

import pandas
import sys


#%%
help_cmd_str = "todo"
try:
    coloc_tsv_gz_path = sys.argv[1]
    explained_tsv_path = sys.argv[2]
    explained_perc_tsv_path = sys.argv[3]
    if len(sys.argv) > 4:
        print("""Two many arguments!
        {}""".format(help_cmd_str))
        sys.exit(1)
except IndexError:
    print("""Argument missing!
    {}""".format(help_cmd_str))
    sys.exit(1)

#%%
# coloc_tsv_gz_path = os.path.join(PathManager.get_outdir_path(), "gwas413/genome/5e-08/1000000/annotate_db.py/annotated.tsv.gz")
df = pandas.read_csv(coloc_tsv_gz_path, sep="\t")

#%%
pp_cutoff = 0.8  # min posterior prob for H4 or H3
pval_cutoff = 5e-8  # min posterior prob for H4 or H3

#%%
df = df.loc[df['gwas_pvalue'] < pval_cutoff]
df.sort_values('PP.H4.abf', ascending=False, inplace=True)
df.drop_duplicates(['chrom', 'pos', 'rsid', 'ref', 'alt', 'cytoband', 'gwas_identifier'], keep='first', inplace=True)
df = df[['chrom', 'pos', 'rsid', 'ref', 'alt', 'cytoband', 'gwas_beta', 'gwas_pvalue', 'gwas_identifier', 'gwas_trait', 'gwas_category', 'PP.H4.abf']]

#%%
out_columns = ['gwas_identifier', 'gwas_trait', 'gwas_category', 'loci_n', 'loci_h4_n']
out_df = pandas.DataFrame(columns = out_columns)
for gwas_identifier in sorted(df['gwas_identifier'].unique()):
    df1 = df.loc[df['gwas_identifier'] == gwas_identifier]
    gwas_out_df = df1[['gwas_identifier', 'gwas_trait', 'gwas_category']].drop_duplicates()
    gwas_out_df['loci_n'] = df1.shape[0]
    gwas_out_df['loci_h4_n'] = (df1['PP.H4.abf'] >= pp_cutoff).sum()
    out_df = pandas.concat([out_df, gwas_out_df])

out_df['explained_perc'] = (out_df['loci_h4_n']/out_df['loci_n']*100).apply(lambda x: round(x, 1))

#% Write
df.sort_values(['gwas_identifier', 'chrom', 'pos', 'rsid', 'ref', 'alt', 'cytoband', ], inplace=True)
df.to_csv(explained_tsv_path, sep="\t", index=False)
out_df.to_csv(explained_perc_tsv_path, sep="\t", index=False)

