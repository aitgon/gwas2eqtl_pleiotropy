Copy input data to the config folder; eg out/config/genome/5e-08/1000000

# Snakemake

~~~
snakemake --cores 1 -p -d ${PWD} -s tools/Snakefile2.yml --config outdir=out/gwas413/genome/5e-08/1000000 raw_coloc_tsv=config/coloc413.tsv gwas_cat_ods=config/gwas413.ods etissue_cat_ods=config/etissue_category.ods upper_var_gwas_cat_count=5 public_data_dir=/home/gonzalez/Software/public
~~~
