Copy input data to the config folder; eg out/config/genome/5e-08/1000000

# Snakemake

~~~
snakemake --cores all -p -d ${PWD} -s tools/Snakefile2.yml --config outdir=out/gwas413/genome/5e-08/1000000 raw_coloc_tsv=config/coloc413.tsv gwas_cat_ods=config/gwas413.ods etissue_cat_ods=config/etissue_category.ods upper_var_gwas_cat_count=5 --rerun-incomplete -n
~~~

# Command by commande 

Create h08 file

~~~
python scripts/filter_h4.py
~~~

Compute rsid, egene, etissue and gwas counts

per rsid

~~~
python scripts/cmpt_count_per_rsid.py  
~~~

~~~
python scripts/cmpt_count_per_bin.py
~~~

Histograms of GWAS subcategory, egenes and etissues

~~~
python scripts/plt_hist_gwas_etissue_egene.py
~~~

Scatter plot of disease, etissue and egene per RSID

~~~
python scripts/plt_scttr_count_per_rsid_gwas.py
python scripts/plt_scttr_count_per_rsid_etissue.py
python scripts/plt_scttr_count_per_rsid_egene.py
~~~
sudo apt install libreoffice
Boxplots of GWAS subcategory counts with egenes and etissues

~~~
python scripts/plt_bxplt_gwas_egene_etissue.py
~~~

Compute tissue enrichment

input: coloc p_h4>0.8 ( "out/cmpt_count_per_rsid.py/coloc08.tsv" ) 
output: "out/cmpt_tissue_enrich.py"

~~~
python scripts/cmpt_tissue_enrich.py
~~~

List pleiotropic regions

~~~
python scripts/cmpt_pleiotropic_regions.py
~~~

Compute pleiotropic regions

~~~
python scripts/cmpt_pleiotropic_regions.py 
~~~
