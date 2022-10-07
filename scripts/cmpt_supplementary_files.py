import pathlib

from gwas2eqtl_pleiotropy.PathManager import PathManager

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
st_writer = pandas.ExcelWriter(supp_tab_path, engine='xlsxwriter')
sheet_counter = 1

wdir_path = os.path.join(PathManager.get_project_path(), wdir_path)

#%%
#
sheet_name_lst = ['ST1']
description_lst = ['Classification of eQTL tissues and cell types. The first six were downloaded from the EBI eQTL repository. '
                   'The 7th column "etissue_class" is used here to compute tissue diversity']
#
sheet_name_lst.append("ST2")
description_lst.append("Metadata and classification of GWAS")
#
sheet_name_lst.append("ST4")
description_lst.append("Count and list of GWAS phenotypes, egenes and etissues for each eQTL/GWAS variant")
#
sheet_name_lst.append("ST5")
description_lst.append("Count and list of GWAS phenotypes for each pleiotropic region")
#
sheet_name_lst.append("ST6")
description_lst.append("Percentage of explained GWAS loci")

capt_df = pandas.DataFrame({'Supp. Tab.': sheet_name_lst, 'Description': description_lst})
capt_df.to_excel(st_writer, sheet_name='Table descrip.', index=False, header=True)

#%% ST1
sheet_name = 'ST1'
etissue_class_ods_path = os.path.join(PathManager.get_project_path(), "config", "etissue_class.ods")
st_df = pandas.read_excel(etissue_class_ods_path, index_col=None, header=0)
st_df.drop(['Unnamed: 8', 'etissue_class.1', 'count'], axis=1, inplace=True)
st_df.to_excel(st_writer, sheet_name=sheet_name, index=False, header=True)
sheet_counter += 1

#%% ST2
sheet_name = 'ST2'
tsv_path = os.path.join(wdir_path, "annot_gwas_metadata.py/gwas413_metadata.tsv")
st_df = pandas.read_csv(tsv_path, sep="\t", header=0)
st_df.sort_values(by=st_df.columns.tolist(), inplace=True)
st_df.to_excel(st_writer, sheet_name=sheet_name, index=False, header=True)
sheet_counter += 1

#%% ST4
sheet_name = 'ST4'
tsv_path = os.path.join(wdir_path, "cmpt_count_per_rsid.py/count_per_rsid_gwas.tsv")
gwas_df = pandas.read_csv(tsv_path, sep="\t", header=0)

tsv_path = os.path.join(wdir_path, "cmpt_count_per_rsid.py/count_per_rsid_egene.tsv")
egene_df = pandas.read_csv(tsv_path, sep="\t", header=0)

tsv_path = os.path.join(wdir_path, "cmpt_count_per_rsid.py/count_per_rsid_etissue.tsv")
etissue_df = pandas.read_csv(tsv_path, sep="\t", header=0)

st_df = pandas.merge(gwas_df, egene_df, on=['chrom', 'cytoband', 'pos', 'rsid'])
st_df = pandas.merge(st_df, etissue_df, on=['chrom', 'cytoband', 'pos', 'rsid'])
st_df = st_df[['chrom', 'cytoband', 'pos', 'rsid', 'gwas_class_count', 'gwas_class_lst', 'egene_count', 'egene_symbol_lst', 'etissue_label_count', 'etissue_subcategory_lst', 'egene_lst']]
st_df.sort_values(['gwas_class_count', 'chrom', 'pos', 'rsid'], ascending=[False, True, True, True], inplace=True)
st_df['gwas_class_lst'] = st_df['gwas_class_lst'].str.replace(',', ', ')
st_df['egene_symbol_lst'] = st_df['egene_symbol_lst'].str.replace(',', ', ')
st_df['etissue_subcategory_lst'] = st_df['etissue_subcategory_lst'].str.replace(',', ', ')
st_df['egene_lst'] = st_df['egene_lst'].str.replace(',', ', ')
st_df.to_excel(st_writer, sheet_name=sheet_name, index=False, header=True)
sheet_counter += 1

#%% ST5
sheet_name = 'ST5'
tsv_path = os.path.join(wdir_path, "cmpt_pleiotropic_regions.py/region_window_100000.tsv")
st_df = pandas.read_csv(tsv_path, sep="\t", header=0)
st_df['egene'] = None
st_df['egene_symbol'] = None
st_df['etissue'] = None

tsv_path = os.path.join(wdir_path, "filter_h4.py/h4.tsv")
h4annot_df = pandas.read_csv(tsv_path, sep="\t", header=0)

for rowi, row in st_df.iterrows():
    chrom = row['chrom']
    start = row['start']
    end = row['end']
    egene_symbol_lst = h4annot_df.loc[(h4annot_df['chrom'] >= chrom) & (h4annot_df['pos'] >= start) & (h4annot_df['pos'] <= end), 'egene_symbol'].unique()
    egene_lst = h4annot_df.loc[(h4annot_df['chrom'] >= chrom) & (h4annot_df['pos'] >= start) & (h4annot_df['pos'] <= end), 'egene'].unique()
    # import pdb; pdb.set_trace()
    etissue_lst = h4annot_df.loc[(h4annot_df['chrom'] >= chrom) & (h4annot_df['pos'] >= start) & (h4annot_df['pos'] <= end), 'etissue_class'].unique()
    st_df.at[rowi, 'egene'] = ", ".join(egene_lst)
    st_df.at[rowi, 'egene_symbol'] = ", ".join(egene_symbol_lst)
    st_df.at[rowi, 'etissue'] = ", ".join(etissue_lst)

st_df.sort_values(by=['gwas_class_count', 'chrom', 'start'], ascending=[False, True, True], inplace=True)
st_df['gwas_class_lst'] = st_df['gwas_class_lst'].str.replace(',', ', ')
st_df = st_df[['chrom', 'cytoband', 'start', 'end', 'gwas_class_count', 'gwas_class_lst', 'egene_symbol', 'etissue', 'egene']]
st_df.to_excel(st_writer, sheet_name=sheet_name, index=False, header=True)
sheet_counter += 1

#%% ST6
sheet_name = 'ST6'
tsv_path = os.path.join(wdir_path, "cmpt_explained_loci.py/explained_perc.tsv")
st_df = pandas.read_csv(tsv_path, sep="\t", header=0)
st_df.sort_values(by=st_df.columns.tolist(), inplace=True)
st_df.to_excel(st_writer, sheet_name=sheet_name, index=False, header=True)
sheet_counter += 1

#%%
st_writer.sheets['Table descrip.'].activate()
st_writer.save()

#%% ST3
sheet_counter = 3
sheet_name = 'ST{}'.format(sheet_counter)
tsv_path = os.path.join(wdir_path, "filter_h4.py/h4.tsv")
h4_annot_df = pandas.read_csv(tsv_path, sep="\t", header=0)
h4_xlsx_path = os.path.join(os.path.dirname(supp_tabl_xlsx_path), "ST{}_gwas2eqtl_h4.xlsx".format(sheet_counter))
writer_coloc_h4 = pandas.ExcelWriter(h4_xlsx_path, engine='xlsxwriter')
h4_annot_df.sort_values(by=h4_annot_df.columns.tolist(), inplace=True)
with pandas.ExcelWriter(h4_xlsx_path) as fout_h4_annot:
    h4_annot_df.to_excel(fout_h4_annot, sheet_name=sheet_name, index=False)
sheet_counter += 1