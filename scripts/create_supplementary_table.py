import pathlib

from eqtl2gwas_pleiotropy.PathManager import PathManager

import os
import pandas
import sys

#%%
help_cmd_str = "todo"
try:
    wdir_path = sys.argv[1]
    supp_tabl_xlsx_path = sys.argv[2]
    if len(sys.argv) > 3:
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
writer = pandas.ExcelWriter(supp_tab_path, engine='xlsxwriter')
sheet_counter = 1

wdir_path = os.path.join(PathManager.get_project_path(), wdir_path)

#%%
#
sheet_name_lst = ['ST1']
description_lst = ['Classification of eQTL tissues and cell types. The first six were downloaded from the EBI eQTL repository. '
                   'The 7th column "etissue_category" is used here to compute tissue diversity']
#
sheet_name_lst.append("ST2")
description_lst.append("Metadata and classification of GWAS")
#
sheet_name_lst.append("ST3")
description_lst.append("Count of GWAS phenotypes of eQTL/GWAS variants")

capt_df = pandas.DataFrame({'Supp. Tab.': sheet_name_lst, 'Description': description_lst})
capt_df.to_excel(writer, sheet_name='Table descrip.', index=False, header=True)

#%%
sheet_name = 'ST{}'.format(sheet_counter)
sheet_counter += 1
etissue_category_ods_path = os.path.join(PathManager.get_project_path(), "config", "etissue_category.ods")
st_df = pandas.read_excel(etissue_category_ods_path, index_col=None, header=0)
st_df.drop(['Unnamed: 7', 'etissue_category.1', 'count'], axis=1, inplace=True)
st_df.to_excel(writer, sheet_name=sheet_name, index=False, header=True)

#%%
sheet_name = 'ST{}'.format(sheet_counter)
sheet_counter += 1
tsv_path = os.path.join(wdir_path, "annot_gwas_metadata.py/gwas413_metadata.tsv")
st_df = pandas.read_csv(tsv_path, sep="\t", header=0)
st_df.sort_values(by=st_df.columns.tolist(), inplace=True)
st_df.to_excel(writer, sheet_name=sheet_name, index=False, header=True)

#%%
sheet_name = 'ST{}'.format(sheet_counter)
sheet_counter += 1
tsv_path = os.path.join(wdir_path, "cmpt_count_per_rsid.py/count_per_rsid_gwas.tsv")
st_df = pandas.read_csv(tsv_path, sep="\t", header=0)
st_df['gwas_category_lst'] = st_df['gwas_category_lst'].str.replace(',', ', ')
# st_df['gwas_category_lst'] = st_df['gwas_category_lst'].str.split(',')
# st_df = st_df.explode('gwas_category_lst')
# st_df['gwas_category'] = 1
# st_df = st_df.pivot_table(index=['chrom', 'cytoband', 'pos', 'rsid', 'gwas_category_count'], columns='gwas_category_lst', values='gwas_category', fill_value=0).reset_index()
# st_df.sort_values(['gwas_category_count', 'chrom', 'pos'], inplace=True, ascending=[False, True, True])
st_df.to_excel(writer, sheet_name=sheet_name, index=False, header=True)

#%%
writer.sheets['Table descrip.'].activate()
writer.save()