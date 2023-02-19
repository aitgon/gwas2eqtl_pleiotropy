import os
import pandas
import pathlib
import seaborn
import sys

from matplotlib import pyplot as plt


#%%
help_cmd_str = "todo"
try:
    snp_pp_h4 = float(sys.argv[1])
    sa_url = sys.argv[2]
    count_per_rsid_gwas_egene_etissue_ods = sys.argv[3]
    track_tsv = sys.argv[4]
    if len(sys.argv) > 5:
        print("""Two many arguments!
        {}""".format(help_cmd_str))
        sys.exit(1)
except IndexError:
    print("""Argument missing!
    {}""".format(help_cmd_str))
    sys.exit(1)

outdir_path = os.path.dirname(track_tsv)
pathlib.Path(outdir_path).mkdir(parents=True, exist_ok=True)

df = pandas.read_excel(count_per_rsid_gwas_egene_etissue_ods, engine='odf')
df = df[['chrom', 'pos38', 'rsid', 'gwas_category_count']].copy()
df['chromStart'] = df['pos38'] - 1
df['chrom'] = 'chr' + df['chrom'].astype(str)
df['rsid'] = 'rs' + df['rsid'].astype(str)
df.rename({'chrom': '#chrom', 'pos38': 'chromEnd', 'rsid': 'name', 'gwas_category_count': 'score'}, axis=1, inplace=True)
df['score'] = df['score'] * 1000 / 4
# import pdb; pdb.set_trace()

track_config_str = """browser position chr5:132239646-132497907
browser pack Pleiotropic_eQTLs
track name=Pleiotropic_eQTLs description="Pleiotropic eQTLs" visibility=2 color=0,128,0 useScore=1
"""
with open(track_tsv, 'w') as fout:
    fout.write(track_config_str)

df = df[['#chrom', 'chromStart', 'chromEnd', 'name', 'score']]
df.to_csv(track_tsv, sep='\t', header=True, mode='a', index=False)
