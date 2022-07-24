import os

import pandas

#%%
count_per_rsid_gwas_tsv_path = "out/gwas413/genome/5e-08/1000000/cmpt_count_per_rsid.py/count_per_rsid_gwas.tsv"
vep0_path = "out/gwas413/genome/5e-08/1000000/cmpt_vep_consequence.py/vep_pleio{}.tsv"

#%%
df = pandas.read_csv(count_per_rsid_gwas_tsv_path, sep="\t")
cmd_lst = []
for pleio in df['gwas_category_count'].unique():
    vep_path = vep0_path.format(pleio)
    vdf = pandas.read_csv(vep_path, sep="\t", comment='#', header=None)
    vdf = vdf.loc[vdf[0].drop_duplicates(keep='first').index]
    vdf = vdf[[0, 6]]
    vdf[6] = vdf[6].str.split(',', expand=True)[0]
    print(vdf[[0,6]])
    import pdb; pdb.set_trace()

#%%
df1 = pandas.read_csv(pleio1_path, sep="\t")
df2 = pandas.read_csv(pleio2_path, sep="\t")
df3 = pandas.read_csv(pleio3_path, sep="\t")
df4 = pandas.read_csv(pleio4_path, sep="\t")

#%%
df1['Consequence'] = df1['Consequence'].str.split(',', expand=True)[0]
df2['Consequence'] = df2['Consequence'].str.split(',', expand=True)[0]
df3['Consequence'] = df3['Consequence'].str.split(',', expand=True)[0]
df4['Consequence'] = df4['Consequence'].str.split(',', expand=True)[0]

#%%
dfu1 = df1[['#Uploaded_variation', 'Consequence']].drop_duplicates()
dfu2 = df2[['#Uploaded_variation', 'Consequence']].drop_duplicates()
dfu3 = df3[['#Uploaded_variation', 'Consequence']].drop_duplicates()
dfu4 = df4[['#Uploaded_variation', 'Consequence']].drop_duplicates()

#%%
stat_df=dfu1.groupby('Consequence').count()
stat_df.columns = ["Count"]

#%%
out_df = pandas.DataFrame({'Consequence': stat1_df.index, 'Category nb.': [1]*stat1_df.shape[0], 'Perc': round(stat1_df['Count']/stat1_df['Count'].sum()*100)})

#%%

stat_df=dfu2.groupby('Consequence').count()
stat_df.columns = ["Count"]
stat_long_df = pandas.DataFrame({'Consequence': stat1_df.index, 'Category nb.': [2]*stat1_df.shape[0], 'Perc': round(stat1_df['Count']/stat1_df['Count'].sum()*100)})

#%%
stat3_df=dfu3.groupby('Consequence').count()
stat3_df.columns = ["Count"]
stat3_df['Perc3'] = round(stat3_df['Count']/stat3_df['Count'].sum()*100)

#%%
stat4_df=dfu4.groupby('Consequence').count()
stat4_df.columns = ["Count"]
stat4_df['Perc4'] = round(stat4_df['Count']/stat4_df['Count'].sum()*100)

#%%
m_df = pandas.merge(stat1_df['Perc1'], stat2_df['Perc2'], left_index=True, right_index=True, how='outer')
m_df = pandas.merge(m_df, stat3_df['Perc3'], left_index=True, right_index=True, how='outer')
m_df = pandas.merge(m_df, stat4_df['Perc4'], left_index=True, right_index=True, how='outer')

#%%
m_df.fillna(0, inplace=True)
m_df = m_df.loc[m_df.sum(axis=1)>10]
m_df = m_df.sort_index()

#%%
