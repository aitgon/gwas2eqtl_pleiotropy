import os
import pandas
import pathlib
import sys


#%%
help_cmd_str = "todo"
try:
    gwas_trait_ods = sys.argv[1]
    gwas_category_ods = sys.argv[2]
    etissue_category_ods = sys.argv[3]
    perc_tophits_eqtl_tsv = sys.argv[4]
    count_per_rsid_gwas_egene_etissue_ods = sys.argv[5]
    region_window_100000_ods = sys.argv[6]
    table_s1_xlsx = sys.argv[7]
    table_s2_xlsx = sys.argv[8]
    table_s3_xlsx = sys.argv[9]
    table_s4_xlsx = sys.argv[10]
    table_s5_xlsx = sys.argv[11]
    if len(sys.argv) > 12:
        print("""Two many arguments!
        {}""".format(help_cmd_str))
        sys.exit(1)
except IndexError:
    print("""Argument missing!
    {}""".format(help_cmd_str))
    sys.exit(1)

#%% ST1
sheet_name = 'ST1'
st_df = pandas.read_excel(etissue_category_ods, index_col=None, header=0)
with pandas.ExcelWriter(table_s1_xlsx) as writer:
    st_df.to_excel(writer, index=False, header=True)

#%% ST2
gwas_trait_df = pandas.read_excel(gwas_trait_ods)
gwas_category_df = pandas.read_excel(gwas_category_ods)

gwas_trait_df.drop(['query_ontology', 'query_term', 'ontology_iri', 'category'], axis=1, inplace=True)
gwas_category_df.drop(['query_ontology', 'query_term', 'ontology_iri', 'category'], axis=1, inplace=True)

gwas_trait_df.rename({'id': 'gwas_id', 'trait': 'gwas_trait', 'ontology_id': 'gwas_trait_ontology_id', 'ontology_term': 'gwas_trait_ontology_term'}, axis=1, inplace=True)
gwas_category_df.rename({'id': 'gwas_id', 'trait': 'gwas_trait', 'ontology_id': 'gwas_category_ontology_id', 'ontology_term': 'gwas_category_ontology_term'}, axis=1, inplace=True)

st_df = gwas_trait_df.merge(gwas_category_df, on=['gwas_id', 'gwas_trait'])
with pandas.ExcelWriter(table_s2_xlsx) as writer:
    st_df.to_excel(writer, index=False, header=True)

#%% ST3 ST3_gwas_trait_perc_explained_loci
perc_tophits_eqtl_df = pandas.read_csv(perc_tophits_eqtl_tsv, sep="\t", index_col='gwas_id')
with pandas.ExcelWriter(table_s3_xlsx) as writer:
    perc_tophits_eqtl_df.to_excel(writer, index=True, header=True)

#%% ST4
count_per_rsid_gwas_egene_etissue_df = pandas.read_excel(count_per_rsid_gwas_egene_etissue_ods, index_col=None, header=0)
with pandas.ExcelWriter(table_s4_xlsx) as writer:
    count_per_rsid_gwas_egene_etissue_df.to_excel(writer, index=False, header=True)

#%% ST5
sheet_name = 'ST5'
# # region_window_100000_tsv_path = os.path.join(wdir_path, "cmpt_pleiotropic_regions.py/region_window_100000.tsv")
# count_per_region_df = pandas.read_csv(region_window_100000_ods, sep="\t", header=0)
count_per_region_df = pandas.read_excel(region_window_100000_ods, index_col=None, header=0)
count_per_region_df['eqtl_gene_symbol'] = None
count_per_region_df['etissue_category_term'] = None
count_per_region_df['egene'] = None
#
for rowi, row in count_per_region_df.iterrows():
    chrom = row['chrom']
    start = row['start']
    end = row['end']
    eqtl_gene_symbol_ser = count_per_rsid_gwas_egene_etissue_df.query(
        'chrom=={chrom} & pos38>={start} & pos38<={end}'.format(chrom=chrom, start=start, end=end))['eqtl_gene_symbol_lst']

    # eqtl_gene_symbol_lst = eqtl_gene_symbol_ser.str.split(';').explode().tolist()
    # eqtl_gene_symbol_lst = sorted(eqtl_gene_symbol_ser.str.split(';').explode().unique())
    eqtl_gene_symbol_ser = eqtl_gene_symbol_ser.str.split(';').explode()
    eqtl_gene_symbol_ser = eqtl_gene_symbol_ser[~eqtl_gene_symbol_ser.isna()]
    eqtl_gene_symbol_lst = sorted(eqtl_gene_symbol_ser.unique())
    count_per_region_df.loc[rowi, 'eqtl_gene_symbol'] = ', '.join(eqtl_gene_symbol_lst)
    #
    egene_ser = count_per_rsid_gwas_egene_etissue_df.query('chrom=={chrom} & pos38>={start} & pos38<={end}'.format(chrom=chrom, start=start, end=end))['egene_lst']
    egene_lst = sorted(egene_ser.str.split(';').explode().unique())
    count_per_region_df.loc[rowi, 'egene'] = ', '.join(egene_lst)
    #
    # import pdb; pdb.set_trace()
    etissue_category_ser = count_per_rsid_gwas_egene_etissue_df.query('chrom=={chrom} & pos38>={start} & pos38<={end}'.format(chrom=chrom, start=start, end=end))['etissue_category_term_lst']
    etissue_category_term_lst = sorted(etissue_category_ser.str.split(';').explode().unique())
    count_per_region_df.loc[rowi, 'etissue_category_term'] = ', '.join(etissue_category_term_lst)

count_per_region_df.sort_values(by=['gwas_category_count', 'chrom', 'start'], ascending=[False, True, True], inplace=True)
# import pdb; pdb.set_trace()
count_per_region_df['gwas_category'] = count_per_region_df['gwas_category'].str.replace(',', ', ')
count_per_region_df = count_per_region_df[['chrom', 'cytoband', 'start', 'end', 'gwas_category_count', 'gwas_category', 'eqtl_gene_symbol', 'etissue_category_term', 'egene']]

with pandas.ExcelWriter(table_s5_xlsx) as writer:
    count_per_region_df.to_excel(writer, index=False, header=True)

