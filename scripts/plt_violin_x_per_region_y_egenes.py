import shlex
import subprocess
import sys

import seaborn
from statannot import add_stat_annotation

from gwas2eqtl_pleiotropy.Logger import Logger
from gwas2eqtl_pleiotropy.PathManager import PathManager
from matplotlib import pyplot as plt
from gwas2eqtl_pleiotropy.constants import tick_fontsize, label_fontsize, scatter_dot_size, dpi

import os
import pandas
import pathlib

plt.rcParams["figure.figsize"] = (8, 6)
from gwas2eqtl_pleiotropy.constants import seaborn_theme_dic
seaborn.set_theme(**seaborn_theme_dic)

#%%
help_cmd_str = "todo"
try:
    region_window_100000_bed_path = sys.argv[1]
    annotate_h4_bed_path = sys.argv[2]
    egenes_per_window_png_path = sys.argv[3]
    if len(sys.argv) > 4:
        print("""Two many arguments!
        {}""".format(help_cmd_str))
        sys.exit(1)
except IndexError:
    print("""Argument missing!
    {}""".format(help_cmd_str))
    sys.exit(1)

outdir_path = os.path.join(os.path.dirname(egenes_per_window_png_path))
pathlib.Path(outdir_path).mkdir(parents=True, exist_ok=True)

#%%
intersect_bed_path = os.path.join(outdir_path, "egenes_per_window.bed")
cmd_stf = "bedtools intersect -a {a_path} -b {b_path} -wa -wb"
cmd = cmd_stf.format(**{'a_path': region_window_100000_bed_path, 'b_path': annotate_h4_bed_path, 'output_bed': intersect_bed_path})
Logger.info(cmd)
with open(intersect_bed_path, 'w') as fout:
    result = subprocess.run(shlex.split(cmd), stdout=fout)

#%%
df = pandas.read_csv(intersect_bed_path, sep='\t', header=None)
df = df[[0, 1, 2, 3, 11, 12]].drop_duplicates()
df.columns = ['chrom', 'start', 'end', 'gwas_class_count', 'egene', 'egene_symbol']
df['start'] = df['start'] + 1
df['length'] = df['end'] - df['start']
df.loc[df['length']==0, 'length']=1

#%%
count_df = df[['chrom', 'start', 'end', 'length', 'gwas_class_count', 'egene']].groupby(['chrom', 'start', 'end', 'length', 'gwas_class_count']).count()
# count_df['egene'] = count_df['egene']/count_df.index.get_level_values('length')
count_df.reset_index(inplace=True)
count_df.to_csv("t.tsv", sep='\t')
#%%
m_df = count_df[['gwas_class_count', 'egene']].drop_duplicates()

# import pdb; pdb.set_trace()
#%%
order = sorted(m_df['gwas_class_count'].unique())
xticklabels = order.copy()
# xticklabels[-1] = '≥{}'.format(order[-1])
title = "eGenes per region"
xlabel = "GWAS category count"
ylabel = "eGenes count"
y = "egene"

#%%
box_pairs = [(1, i) for i in order if not i==1 ]
# import pdb; pdb.set_trace()
ax = seaborn.violinplot(x="gwas_class_count", y=y, data=m_df, order=order, palette="rocket_r")
test_results = add_stat_annotation(ax, data=m_df, x="gwas_class_count", y=y, order=order,
                                   box_pairs=box_pairs,
                                   test='Mann-Whitney', text_format='star',
                                   loc='inside', verbose=2)

plt.title(title, fontsize=label_fontsize)
plt.xlabel(xlabel, fontsize=label_fontsize)
plt.xticks(fontsize=tick_fontsize, rotation=0)
plt.ylabel(ylabel, fontsize=label_fontsize)
plt.yticks(fontsize=tick_fontsize)
ax.set_xticklabels(xticklabels)

plt.tight_layout()
plt.savefig(egenes_per_window_png_path)
plt.close()
