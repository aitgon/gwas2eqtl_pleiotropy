from eqtl2gwas_pleiotropy.constants import h4_cutoff

import os
import pandas
import pathlib
import sys

#%%
help_cmd_str = "todo"
try:
    annotated_tsv_gz_path = sys.argv[1]
    h4_tsv_gz_path = sys.argv[2]
    if len(sys.argv) > 3:
        print("""Two many arguments!
        {}""".format(help_cmd_str))
        sys.exit(1)
except IndexError:
    print("""Argument missing!
    {}""".format(help_cmd_str))
    sys.exit(1)

outdir_path = os.path.dirname(h4_tsv_gz_path)
pathlib.Path(outdir_path).mkdir(parents=True, exist_ok=True)

#%% input
if not os.path.isfile(annotated_tsv_gz_path):
    print("input file does not exit")
    sys.exit(1)

#%%
coloc_df = pandas.read_csv(annotated_tsv_gz_path, sep="\t")
coloc_df.drop(['id'], inplace=True, axis=1)
coloc_df = coloc_df.loc[coloc_df['PP.H4.abf'] >= h4_cutoff, ]

#%% how many coloc loci?
loci_h4_df = coloc_df.sort_values(by=['PP.H4.abf', 'SNP.PP.H4'], ascending=False)
loci_h4_df = loci_h4_df.drop_duplicates(subset=['PP.H4.abf', 'coloc_window', 'nsnps'], keep='first')
loci_coloc_h4_path = os.path.join(outdir_path, "loci_h4.tsv")
loci_h4_df.to_csv(loci_coloc_h4_path, sep="\t", index=False)

#%% how many causal variants, SNP.PP.H4>=0.5?
variants_h4_df = coloc_df.sort_values(by=['SNP.PP.H4'], ascending=False)
variants_h4_df = variants_h4_df.drop_duplicates(subset=['chrom', 'pos', 'cytoband', 'rsid'], keep='first')
variants_h4_df = variants_h4_df.loc[variants_h4_df['SNP.PP.H4'] >= 0.5, ]
variants_h4_path = os.path.join(outdir_path, "variants_h4.tsv")
variants_h4_df.to_csv(variants_h4_path, sep="\t", index=False)

#%%
coloc_df.to_csv(h4_tsv_gz_path, sep="\t", index=False)
