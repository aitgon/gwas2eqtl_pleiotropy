# Snakemake

~~~
snakemake --cores all -p -d ${PWD} -s tools/Snakefile.yml
~~~

# Command by commande 

Create h08 file

~~~
python scripts/filterh4.py
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
