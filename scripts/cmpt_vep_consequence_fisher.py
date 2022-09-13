import os
import sys

import pandas
import seaborn
from statannot import add_stat_annotation
from statsmodels.stats import multitest as multitest
from matplotlib import pyplot as plt
import numpy

#%%
from scipy.stats import fisher_exact

from gwas2eqtl_pleiotropy.Logger import Logger
from gwas2eqtl_pleiotropy.constants import tick_fontsize, label_fontsize, dpi

upper_var_gwas_cat_count = 5
vep1_path = "out/gwas413/genome/5e-08/1000000/cmpt_vep.py/vep_pleio1.tsv"

#%%
help_cmd_str = "todo"
try:
    vep_input_path = sys.argv[1]
    vep_output_path = sys.argv[2]
    upper_var_gwas_cat_count = int(sys.argv[3])
    consequence_tsv_path = sys.argv[4]
    if len(sys.argv) > 5:
        print("""Two many arguments!
        {}""".format(help_cmd_str))
        sys.exit(1)
except IndexError:
    print("""Argument missing!
    {}""".format(help_cmd_str))
    sys.exit(1)


#%%
vep_input_column_lst = ['chrom', 'start', 'end', 'alleles', 'strand', 'rsid', 'gwas_category_count']
vep_input_df = pandas.read_csv(vep_input_path, sep=" ", header=None, names=vep_input_column_lst)
columns = ['#Uploaded_variation', 'Location', 'Allele', 'Gene', 'Feature', 'Feature_type', 'Consequence', 'cDNA_position', 'CDS_position', 'Protein_position', 'Amino_acids', 'Codons', 'Existing_variation', 'Extra']
vep_output_df = pandas.read_csv(vep_output_path, sep="\t", comment='#', header=None, names=columns)


#%%
df = vep_input_df.merge(vep_output_df, left_on='rsid', right_on='#Uploaded_variation', how='left')

#%%
# import pdb; pdb.set_trace()
# df[['Consequence0', 'Consequence1', 'Consequence2', 'Consequence3']] = df['Consequence'].str.split(',', expand=True)
df['Consequence0'] = df['Consequence'].str.split(',', expand=True)[0]

out_columns = ['consequence', 'gwas_category_count', 'a_pleio_x_with_consequence', 'b_pleio_1_with_consequence', 'c_pleio_x_wout_consequence', 'd_pleio_1_wout_consequence', 'oddsr', 'p']
out_dic = dict(zip(out_columns, [[] for i in range(len(out_columns))]))
out_df = pandas.DataFrame(out_dic)

#%%
# consequence = 'missense_variant'
# consequence_lst = ['3_prime_UTR_variant', '5_prime_UTR_variant', 'downstream_gene_variant', 'frameshift_variant', 'intergenic_variant', 'intron_variant', 'missense_variant', 'splice_region_variant', 'stop_lost', 'upstream_gene_variant']
# consequence_lst = df['Consequence0'].unique().tolist()
for consequence in sorted(df.loc[~df['Consequence0'].isna(), 'Consequence0'].unique()):
    Logger.info('Consenquence: ' + consequence)
    all_df = vep_input_df[['rsid', 'gwas_category_count']].drop_duplicates()
    cons_df = df.loc[df['Consequence0'] == consequence, ['rsid', 'gwas_category_count']].drop_duplicates()
    cons_df = all_df.merge(cons_df, on=['rsid', 'gwas_category_count'], how='left', indicator=True)
    cons_df.loc[cons_df['gwas_category_count'] > upper_var_gwas_cat_count, 'gwas_category_count'] = upper_var_gwas_cat_count


    #%%
    count_df = cons_df.groupby(['gwas_category_count', '_merge']).size().reset_index()
    count_df = count_df.loc[count_df['_merge'] != 'right_only', ]

    #%%
    # pleio 1, consequence
    bb = count_df.loc[(count_df['gwas_category_count'] == 1) & (count_df['_merge'] == 'both'), 0].values[0]
    # pleio 1, non-consequence
    dd = count_df.loc[(count_df['gwas_category_count'] == 1) & (count_df['_merge'] == 'left_only'), 0].values[0] - bb

    #%%
    # gwas_cat_count = 2
    for gwas_cat_count in [*range(2, upper_var_gwas_cat_count + 1)]:
        # pleio x, consequence
        aa = count_df.loc[(count_df['gwas_category_count'] == gwas_cat_count) & (count_df['_merge'] == 'both'), 0].values[0]
        # pleio x, non-consequence
        cc = count_df.loc[(count_df['gwas_category_count'] == gwas_cat_count) & (count_df['_merge'] == 'left_only'), 0].values[0] - aa

        #%% fisher
        table = numpy.array([[aa, bb], [cc, dd]])
        # import pdb; pdb.set_trace()
        oddsr, ppp = fisher_exact(table, alternative='two-sided')

        #%% perc
        if (bb+dd) == 0 or (aa+cc)==0:
            continue

        #%%
        this_row = [consequence, gwas_cat_count, aa, bb, cc, dd, oddsr, ppp]
        this_row_dic = dict(zip(out_columns, this_row))
        out_df = pandas.concat([out_df, pandas.DataFrame(this_row_dic, index=[0])])

rejected, pcorrected = multitest.fdrcorrection(out_df['p'], alpha=0.05)
out_df['pfdr5perc'] = pcorrected

out_df.to_csv(consequence_tsv_path, sep="\t", index=False, header=True)
