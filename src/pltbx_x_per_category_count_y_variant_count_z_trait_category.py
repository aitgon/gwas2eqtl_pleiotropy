import sqlalchemy
from constants import label_fontsize, tick_fontsize, boxenplot_kws, palette_r, palette
from matplotlib.ticker import MaxNLocator

import os
import pandas
import pathlib
import seaborn
import sys
import matplotlib.pyplot as plt


#%%
plt.rcParams["figure.figsize"] = (8, 6)
from constants import seaborn_theme_dic
seaborn.set_theme(**seaborn_theme_dic)

#%%
help_cmd_str = "todo"
try:
    count_per_rsid_gwas_ods_path = sys.argv[1]
    png_path = sys.argv[2]
    if len(sys.argv) > 3:
        print("""Two many arguments!
        {}""".format(help_cmd_str))
        sys.exit(1)
except IndexError:
    print("""Argument missing!
    {}""".format(help_cmd_str))
    sys.exit(1)


#%%
if not os.path.isfile(count_per_rsid_gwas_ods_path):
    print("input file does not exit")
    sys.exit(1)

outdir_path = os.path.dirname(png_path)
pathlib.Path(outdir_path).mkdir(parents=True, exist_ok=True)

#%%
# count_per_rsid_gwas_df = pandas.read_csv(count_per_rsid_gwas_ods_path, sep="\t")
df = pandas.read_excel(count_per_rsid_gwas_ods_path, engine='odf')
df = df[['rsid', 'gwas_category_lst', 'gwas_category_count']].drop_duplicates()
df['gwas_category_lst'] = df['gwas_category_lst'].str.split(';')
df = df.explode('gwas_category_lst', ignore_index=False)
df.sort_values(by=["gwas_category_lst", "rsid", "gwas_category_count"], inplace=True)

#%%
# order = [str(x) for x in range(1, max(m3_df['gwas_category_count'].unique()) + 1)]
# xticklabels = order.copy()
y = "proportion"
x = "gwas_category_count"

count_data = df['gwas_category_count'].value_counts().reset_index()
count_data.columns = ['gwas_category_count', 'count']
count_data['proportion'] = count_data['count'] / count_data['count'].sum()

# Create the bar plot
plt.figure(figsize=(8, 6))
ax = seaborn.barplot(x=x, y=y, data=count_data, palette=palette)

# Annotate each bar with the actual count of variants
for row_i, row in count_data.iterrows():
    ax.text(row_i, row['proportion'] + 0.01, f"n={int(row['count'])}", ha='center', color='black', fontsize=tick_fontsize)

label_fontsize = 26
plt.title("eQTLs by trait category count", fontsize=label_fontsize)
plt.xlabel("Trait category count", fontsize=label_fontsize)
plt.ylabel("Proportion of eQTLs", fontsize=label_fontsize)
plt.yticks(fontsize=tick_fontsize)
plt.xticks(fontsize=tick_fontsize)
plt.ylim([0, 1])
# plt.yscale("log")

plt.tight_layout()
# png_path = os.path.join(outdir_path, "0proportion_eqtls_by_category_count.png")
plt.savefig(png_path)
plt.close()

for gwas_category_i, gwas_category in enumerate(sorted(df["gwas_category_lst"].unique().tolist()), start=1):
    df2 = df[df["gwas_category_lst"] == gwas_category]

    # Calculate the count and proportion of each GWAS category count
    count_data = df2['gwas_category_count'].value_counts().reset_index()
    count_data.columns = ['gwas_category_count', 'count']
    count_data['proportion'] = count_data['count'] / count_data['count'].sum()

    # Create the bar plot
    plt.figure(figsize=(8, 6))
    ax = seaborn.barplot(x=x, y=y, data=count_data, palette=palette)

    # Annotate each bar with the actual count of variants
    for row_i, row in count_data.iterrows():
        ax.text(row_i, row['proportion'] + 0.01, f"n={int(row['count'])}", ha='center', color='black', fontsize=tick_fontsize)

    label_fontsize = 26
    plt.title("{}".format(gwas_category), fontsize=label_fontsize)
    plt.xlabel("Trait category count", fontsize=label_fontsize)
    plt.ylabel("Proportion of eQTLs", fontsize=label_fontsize)
    plt.yticks(fontsize=tick_fontsize)
    plt.xticks(fontsize=tick_fontsize)
    plt.ylim([0, 1])

    plt.tight_layout()
    png_path = os.path.join(outdir_path, "{}.png".format(gwas_category_i))
    plt.savefig(png_path)
    plt.close()
