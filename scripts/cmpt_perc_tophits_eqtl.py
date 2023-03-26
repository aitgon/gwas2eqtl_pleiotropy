import os
import pandas
import pathlib
import seaborn
import sqlalchemy
import sys

from matplotlib import pyplot as plt

from gwas2eqtl_pleiotropy.constants import tick_fontsize, seaborn_theme_dic

# %%
snp_pp_h4 = "0.5"
db_url = "postgresql://postgres:postgres@0.0.0.0:5435/postgres"
loci_explained_perc_tsv = "out/gwas417/pval_5e-08/r2_0.1/kb_1000/window_1000000/75_50/cmpt_perc_tophits_eqtl.py/perc_tophits_eqtl.tsv"

# #%%
# help_cmd_str = "todo"
# try:
#     snp_pp_h4 = float(sys.argv[1])
#     db_url = sys.argv[2]
#     loci_explained_perc_tsv = sys.argv[3]
#     if len(sys.argv) > 4:
#         print("""Two many arguments!
#         {}""".format(help_cmd_str))
#         sys.exit(1)
# except IndexError:
#     print("""Argument missing!
#     {}""".format(help_cmd_str))
#     sys.exit(1)

#%%
seaborn.set_theme(**seaborn_theme_dic)

#%%
outdir_path = os.path.dirname(loci_explained_perc_tsv)
pathlib.Path(outdir_path).mkdir(exist_ok=True, parents=True)

sql = 'select distinct chrom, pos, nea, ea, tophits.gwas_id, trait as gwas_trait, op.consortium, op.pmid from tophits, open_gwas_info op where op.gwas_id=tophits.gwas_id'
engine = sqlalchemy.create_engine(db_url)
with engine.begin() as conn:
    tophits_df = pandas.read_sql(sqlalchemy.text(sql), con=conn).drop_duplicates()

sql = 'select distinct * from colocpleio where snp_pp_h4>={}'.format(snp_pp_h4)
columns = ['chrom', 'pos38', 'alt', 'eqtl_gene_id', 'gwas_id', 'gwas_trait', 'gwas_trait_ontology_term', 'gwas_category_ontology_term', 'eqtl_id', 'tophits_variant_id']
engine = sqlalchemy.create_engine(db_url)
with engine.begin() as conn:
    coloc_df = pandas.read_sql(sqlalchemy.text(sql), con=conn)
coloc_df = coloc_df[columns].drop_duplicates()
# Remove Locus MHC
coloc_df = coloc_df.query('not (chrom==6 & pos38>=25000000 & pos38<=35000000)')

tophits_df['tophits_variant_id'] = tophits_df['chrom'].astype(str) + '_' + tophits_df['pos'].astype(str) + '_' + tophits_df['nea'] + '_' + tophits_df['ea']
# Remove Locus MHC
tophits_df = tophits_df.query('not (chrom==6 & pos>=25000000 & pos<=35000000)')

m_df = tophits_df[['chrom', 'tophits_variant_id', 'gwas_id']].merge(coloc_df[['chrom', 'tophits_variant_id', 'gwas_id']], on=['chrom', 'tophits_variant_id', 'gwas_id'], how='left', indicator=True).drop_duplicates()
m_df.sort_values(['gwas_id', 'tophits_variant_id'], inplace=True)

# Percentage explained non-MHC
m_df = m_df.groupby(['gwas_id', '_merge']).size().reset_index()
m_df = m_df.loc[m_df['_merge'] != 'right_only']
m_df = m_df.pivot_table(index=['gwas_id'], columns=['_merge'], values=0)

m_df['loci_cnt'] = m_df.apply(sum, axis=1)
m_df.drop(['left_only'], axis=1, inplace=True)
m_df.rename({'both': 'explained_cnt'}, axis=1, inplace=True)
m_df['loci_explained_perc'] = m_df['explained_cnt'] / m_df['loci_cnt'] * 100
m_df['loci_explained_perc'] = m_df['loci_explained_perc'].apply(int)

gwas_annot_df = coloc_df[['gwas_id', 'gwas_trait', 'gwas_trait_ontology_term', 'gwas_category_ontology_term']].drop_duplicates()
m_df.merge(gwas_annot_df, on='gwas_id')
m_df = m_df.merge(gwas_annot_df, on='gwas_id')

xlim = [-5, 105]
xlabel = "Explained loci percentage"
ylabel = "GWAS trait ontology"

# import pdb; pdb.set_trace()
m_df = m_df[['gwas_id', 'gwas_trait', 'gwas_trait_ontology_term', 'gwas_category_ontology_term', 'loci_cnt', 'explained_cnt', 'loci_explained_perc']]
m_df.set_index(['gwas_id'], verify_integrity=True, inplace=True)
m_df.to_csv(loci_explained_perc_tsv, sep='\t', index=True)

#%%
m_df.sort_values(by=['gwas_category_ontology_term', 'gwas_trait_ontology_term', 'gwas_trait'], inplace=True)

#%%
# fig, axes = plt.subplots(4, 1, sharex=True, figsize=(8,8))
# plt.show()
# m2_df = m_df.query('gwas_category_ontology_term=="autoimmune disease" or gwas_category_ontology_term=="cardiovascular disease"')

# #%%
# plt.rcParams["figure.figsize"] = (20, 10)
# order = m2_df['gwas_trait_ontology_term'].tolist()
# seaborn.catplot(data=m2_df, y="gwas_category_ontology_term", x="loci_explained_perc", orient='h', hue="gwas_trait_ontology_term", dodge=True, legend=False)
# plt.legend(loc='upper left', title='Team', bbox_to_anchor=(-1, 1))
# plt.yticks(fontsize=5)
# plt.grid(True)
# plt.subplots_adjust(left=0.4, right=0.9, bottom=0.05, top=0.95)
# plt.show()

#%%
m_df['gwas_trait_ontology_term'] = m_df['gwas_trait_ontology_term'].str.slice(0, 30)

#%%
fig, axes = plt.subplots(9, 4, sharex=False, figsize=(15, 20))
counter = 0
for gwas_category_ontology_term in sorted(m_df['gwas_category_ontology_term'].unique()):
    print(gwas_category_ontology_term)
    m2_df = m_df.query('gwas_category_ontology_term=="{}"'.format(gwas_category_ontology_term))
    order = sorted(m2_df['gwas_trait_ontology_term'].unique())
    # plt.rcParams["figure.figsize"] = (20, 10)
    g = seaborn.stripplot(ax=axes[int(counter/4), counter%4], data=m2_df, y="gwas_trait_ontology_term", x="loci_explained_perc", orient='h', order=order, s=5)
    g.tick_params(axis='y', which='major', labelsize=8)
    g.set(ylabel=None)
    g.set(xlabel="Percent. explain.")
    g.set(title=gwas_category_ontology_term)
    g.set_xticks([0, 25, 50, 75, 100])
    plt.yticks(fontsize=5, rotation=0)
    plt.tight_layout()
    counter = counter + 1
    # if counter == 8:
    #     break

fig.delaxes(axes[8][1])
fig.delaxes(axes[8][2])
fig.delaxes(axes[8][3])

# plt.grid(True)
png_path = os.path.join(outdir_path, "subplots.png".format(gwas_category_ontology_term.replace(' ', '_')))
plt.savefig(png_path)
plt.close()
