import matplotlib.pyplot as plt
import os
import pandas
import pathlib
import seaborn
import sys

from constants import tick_fontsize, dpi
from constants import seaborn_theme_dic
from matplotlib.ticker import MaxNLocator

plt.rcParams["figure.figsize"] = (8, 6)
seaborn.set_theme(**seaborn_theme_dic)

#%%
help_cmd_str = "todo"
try:
    count_per_rsid_gwas_egene_etissue_ods = sys.argv[1]
    hist_rsid_gwas_path = sys.argv[2]
    hist_rsid_egene_path = sys.argv[3]
    hist_rsid_etissue_path = sys.argv[4]
    if len(sys.argv) > 5:
        print("""Two many arguments!
        {}""".format(help_cmd_str))
        sys.exit(1)
except IndexError:
    print("""Argument missing!
    {}""".format(help_cmd_str))
    sys.exit(1)

pathlib.Path(os.path.dirname(hist_rsid_gwas_path)).mkdir(parents=True, exist_ok=True)
pathlib.Path(os.path.dirname(hist_rsid_egene_path)).mkdir(parents=True, exist_ok=True)
pathlib.Path(os.path.dirname(hist_rsid_etissue_path)).mkdir(parents=True, exist_ok=True)

count_df = pandas.read_excel(count_per_rsid_gwas_egene_etissue_ods, engine='odf')

#%%
<<<<<<< HEAD
ylabel = "Proportion"
ylim = [0.0001, 1]
=======
ylabel = "Proportion of eQTLs"
ylim = [0.0001, 0.6]
>>>>>>> 7780b854f4c03d61cfedb2434aa1fd98189fe736
edgecolor = 'k'
linewidth = 2
grid_axis = 'y'
hist_kwargs = {'density': 1, 'edgecolor': edgecolor, 'linewidth': linewidth}
stat = 'proportion'
label_fontsize = 28

#%% gwas
data_ser = count_df['gwas_category_count']
ax = seaborn.histplot(data_ser, stat=stat, discrete=True)

ax.xaxis.set_major_locator(MaxNLocator(integer=True))
plt.grid(visible=True, axis='y')
title = "Trait specificity"
plt.title(title, fontsize=label_fontsize)
plt.xlabel("Trait category count", fontsize=label_fontsize)
plt.xticks(fontsize=tick_fontsize)
plt.ylabel(ylabel, fontsize=label_fontsize)
plt.yticks(fontsize=tick_fontsize)
<<<<<<< HEAD
plt.yscale('log')
plt.ylim(ylim)

=======
# plt.yscale('log')
plt.ylim(ylim)

# Access the container with the bars
bars = ax.containers[0]
# Add counts (n) on top of each bar
for bari, bar in enumerate(bars):
    height = bar.get_height()
    n = int(data_ser.value_counts()[bari+1])
    if height > 0:  # Only label bars with height > 0
        ax.text(
            bar.get_x() + bar.get_width() / 2,  # X position
            height,                             # Y position
            f'{n}',                   # The count text
            ha='center',                        # Horizontal alignment
            va='bottom'                         # Vertical alignment
        )

>>>>>>> 7780b854f4c03d61cfedb2434aa1fd98189fe736
plt.tight_layout()
plt.savefig(hist_rsid_gwas_path, dpi=dpi)
plt.close()

with open(hist_rsid_gwas_path + ".txt", 'w') as fout:
    fout.write('\n'.join([str(h.get_height()) for h in ax.patches]))

#%% eQTL gene
data_ser = count_df['egene_count']
ax = seaborn.histplot(data_ser, stat=stat, discrete=True)

ax.xaxis.set_major_locator(MaxNLocator(integer=True))
plt.grid(visible=True, axis='y')
title = "eQTL gene specificity"
plt.title(title, fontsize=label_fontsize)
plt.xlabel("eQTL gene count", fontsize=label_fontsize)
plt.xticks(fontsize=tick_fontsize)
plt.ylabel(ylabel, fontsize=label_fontsize)
plt.yticks(fontsize=tick_fontsize)
<<<<<<< HEAD
plt.yscale('log')
=======
# plt.yscale('log')
>>>>>>> 7780b854f4c03d61cfedb2434aa1fd98189fe736
plt.ylim(ylim)

plt.tight_layout()
plt.savefig(hist_rsid_egene_path, dpi=dpi)
plt.close()

with open(hist_rsid_egene_path + ".txt", 'w') as fout:
    fout.write('\n'.join([str(h.get_height()) for h in ax.patches]))

#%% eQTL biosample
data_ser = count_df['etissue_category_term_count']
ax = seaborn.histplot(data_ser, stat=stat, discrete=True)

ax.xaxis.set_major_locator(MaxNLocator(integer=True))
plt.grid(visible=True, axis='y')
title = "eQTL gene specificity"
plt.title(title, fontsize=label_fontsize)
plt.xlabel("eQTL biosample count", fontsize=label_fontsize)
plt.xticks(fontsize=tick_fontsize)
plt.ylabel(ylabel, fontsize=label_fontsize)
plt.ylim(ylim)
plt.yticks(fontsize=tick_fontsize)
<<<<<<< HEAD
plt.yscale('log')
=======
# plt.yscale('log')
>>>>>>> 7780b854f4c03d61cfedb2434aa1fd98189fe736

plt.tight_layout()
plt.savefig(hist_rsid_etissue_path, dpi=dpi)
plt.close()

with open(hist_rsid_etissue_path + ".txt", 'w') as fout:
    fout.write('\n'.join([str(h.get_height()) for h in ax.patches]))
