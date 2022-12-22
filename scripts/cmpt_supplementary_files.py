from gwas2eqtl_pleiotropy.PathManager import PathManager

import os
import pandas
import pathlib
import sys


#%%
help_cmd_str = "todo"
try:
    gwas_class_ods_path = sys.argv[1]
    etissue_class_ods_path = sys.argv[2]
    count_per_rsid_gwas_egene_etissue_tsv_path = sys.argv[3]
    region_window_100000_tsv_path = sys.argv[4]
    supp_tabl_xlsx_path = sys.argv[5]
    if len(sys.argv) > 6:
        print("""Two many arguments!
        {}""".format(help_cmd_str))
        sys.exit(1)
except IndexError:
    print("""Argument missing!
    {}""".format(help_cmd_str))
    sys.exit(1)

# Create a Pandas Excel writer using XlsxWriter as the engine.
supp_tab_path = os.path.join(supp_tabl_xlsx_path)
pathlib.Path(os.path.dirname(supp_tab_path)).mkdir(exist_ok=True, parents=True)
st_writer = pandas.ExcelWriter(supp_tab_path, engine='xlsxwriter')
sheet_counter = 1

# wdir_path = os.path.join(PathManager.get_project_path(), wdir_path)

#%%
#
sheet_name_lst = []
description_lst = []

sheet_name_lst.append("ST1")
description_lst.append('Classification of eQTL tissues and cell types. The first six were downloaded from the EBI eQTL repository. '
                   'The 7th column "etissue_class" is used here to compute tissue diversity')

sheet_name_lst.append('ST2')
description_lst.append("Metadata and classification of GWAS")

# ST3 ST3_gwas_trait_perc_explained_loci

sheet_name_lst.append("ST4")
description_lst.append("Count and list of GWAS phenotypes, egenes and etissues for each eQTL/GWAS variant")

sheet_name_lst.append("ST5")
description_lst.append("Count and list of GWAS phenotypes, egene symbols and ENSEMBL IDs and etissue classes for each pleiotropic region")

capt_df = pandas.DataFrame({'Supp. Tab.': sheet_name_lst, 'Description': description_lst})
capt_df.to_excel(st_writer, sheet_name='Table descrip.', index=False, header=True)

#%% ST1
sheet_name = 'ST1'
# etissue_class_ods_path = os.path.join(PathManager.get_project_path(), "config", "etissue_class.ods")
st_df = pandas.read_excel(etissue_class_ods_path, index_col=None, header=0)
st_df.drop(['Unnamed: 8', 'etissue_class.1', 'count'], axis=1, inplace=True)
st_df.to_excel(st_writer, sheet_name=sheet_name, index=False, header=True)

#%% ST2
sheet_name = 'ST2'
gwas_class_df = pandas.read_excel(gwas_class_ods_path)
class_pleio_df = gwas_class_df[['icd10_code_level1.1', 'class_pleiotropy']]
class_pleio_df = class_pleio_df.loc[~class_pleio_df.isna().any(axis=1)]
class_pleio_df.rename({'icd10_code_level1.1': 'icd10_code_level1'}, axis=1, inplace=True)
#
st_df = gwas_class_df[['id', 'trait', 'icd10_code_level1']].merge(class_pleio_df, on='icd10_code_level1')
st_df = st_df.drop_duplicates()
st_df = st_df.sort_values(st_df.columns.tolist())
st_df.drop(['icd10_code_level1'], axis=1, inplace=True)
st_df.to_excel(st_writer, sheet_name=sheet_name, index=False, header=True)

#%% ST3 ST3_gwas_trait_perc_explained_loci

#%% ST4
sheet_name = 'ST4'
# region_window_100000_tsv_path = os.path.join(wdir_path, "cmpt_count_per_rsid.py/count_per_rsid_gwas.tsv")
gwas_df = pandas.read_csv(region_window_100000_tsv_path, sep="\t", header=0)

# tsv_path = os.path.join(wdir_path, "cmpt_count_per_rsid.py/count_per_rsid_egene.tsv")
count_per_rsid_gwas_egene_etissue_df = pandas.read_csv(count_per_rsid_gwas_egene_etissue_tsv_path, sep="\t", header=0)
#
# tsv_path = os.path.join(wdir_path, "cmpt_count_per_rsid.py/count_per_rsid_etissue.tsv")
# etissue_df = pandas.read_csv(tsv_path, sep="\t", header=0)
# import pdb; pdb.set_trace()
# count_per_rsid_df = pandas.merge(gwas_df, egene_df, on=['chrom', 'cytoband', 'pos38', 'rsid'])
# count_per_rsid_df = pandas.merge(count_per_rsid_df, etissue_df, on=['chrom', 'cytoband', 'pos38', 'rsid'])
# count_per_rsid_df = count_per_rsid_df[['chrom', 'cytoband', 'pos38', 'rsid', 'gwas_class_count', 'gwas_class_lst', 'egene_count',
#                'eqtl_gene_symbol_lst', 'etissue_label_count', 'etissue_class_lst', 'egene_lst']]
# count_per_rsid_df.sort_values(['gwas_class_count', 'chrom', 'pos38', 'rsid'], ascending=[False, True, True, True], inplace=True)
# count_per_rsid_df['gwas_class_lst'] = count_per_rsid_df['gwas_class_lst'].str.replace(',', ', ')
# count_per_rsid_df['eqtl_gene_symbol_lst'] = count_per_rsid_df['eqtl_gene_symbol_lst'].str.replace(',', ', ')
# count_per_rsid_df['etissue_class_lst'] = count_per_rsid_df['etissue_class_lst'].str.replace(',', ', ')
# count_per_rsid_df['egene_lst'] = count_per_rsid_df['egene_lst'].str.replace(',', ', ')
count_per_rsid_gwas_egene_etissue_df.to_excel(st_writer, sheet_name=sheet_name, index=False, header=True)

#%% ST5
# sheet_name = 'ST5'
# # region_window_100000_tsv_path = os.path.join(wdir_path, "cmpt_pleiotropic_regions.py/region_window_100000.tsv")
count_per_region_df = pandas.read_csv(region_window_100000_tsv_path, sep="\t", header=0)
count_per_region_df['eqtl_gene_symbol'] = None
count_per_region_df['etissue_class'] = None
count_per_region_df['egene'] = None
#
for rowi, row in count_per_region_df.iterrows():
    chrom = row['chrom']
    start = row['start']
    end = row['end']
    eqtl_gene_symbol_ser = count_per_rsid_gwas_egene_etissue_df.query('chrom=={chrom} & pos38>={start} & pos38<={end}'.format(chrom=chrom, start=start, end=end))['eqtl_gene_symbol_lst']
    eqtl_gene_symbol_lst = sorted(eqtl_gene_symbol_ser.str.split(', ').explode().unique())
    count_per_region_df.loc[rowi, 'eqtl_gene_symbol'] = ', '.join(eqtl_gene_symbol_lst)
    #
    egene_ser = count_per_rsid_gwas_egene_etissue_df.query('chrom=={chrom} & pos38>={start} & pos38<={end}'.format(chrom=chrom, start=start, end=end))['egene_lst']
    egene_lst = sorted(egene_ser.str.split(', ').explode().unique())
    count_per_region_df.loc[rowi, 'egene'] = ', '.join(egene_lst)
    #
    etissue_class_ser = count_per_rsid_gwas_egene_etissue_df.query('chrom=={chrom} & pos38>={start} & pos38<={end}'.format(chrom=chrom, start=start, end=end))['etissue_class_lst']
    etissue_class_lst = sorted(etissue_class_ser.str.split(', ').explode().unique())
    count_per_region_df.loc[rowi, 'etissue_class'] = ', '.join(etissue_class_lst)

count_per_region_df.sort_values(by=['gwas_class_count', 'chrom', 'start'], ascending=[False, True, True], inplace=True)
count_per_region_df['gwas_class_lst'] = count_per_region_df['gwas_class_lst'].str.replace(',', ', ')
count_per_region_df = count_per_region_df[['chrom', 'cytoband', 'start', 'end', 'gwas_class_count', 'gwas_class_lst', 'eqtl_gene_symbol', 'etissue_class', 'egene']]
count_per_region_df.to_excel(st_writer, sheet_name=sheet_name, index=False, header=True)

#%%
st_writer.sheets['Table descrip.'].activate()
# st_writer.save()
st_writer.close()
