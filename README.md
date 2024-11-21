Updated Nov 21, 2024

The GWAS/eQTL colocalization data is the result of the gwas2eqtl repository and can be downloaded from ...

A postgresql database is necessary to annotate the colocalization data.

~~~
cd containers
docker compose --project-name gwas2eqtl_pleiotropy_dev -f docker-compose.yml up --build --force-recreate --remove-orphans -d
cd ..
~~~

This micromamba environment is only used to run apptainer, but any other apptainer executable is also ok.

~~~
micromamba create -n gwas2eqtl_pleiotropy apptainer==1.3.2 -c conda-forge python=3.11
micromamba activate gwas2eqtl_pleiotropy
~~~

This container "containers/gwas2eqtl_pleiotropy.def" is very important to make sure that everything will work.

~~~
mkdir -p results/containers
sudo apptainer build results/containers/gwas2eqtl_pleiotropy.sif  containers/gwas2eqtl_pleiotropy.def
~~~

The workflow "workflows/11snkfl_insrt_postgres.yml" is used to populate the database with the colocalization and annotation data.

The annotation data in tar.gz format is given as argument

The result of this workflow is a "colocpleio" materialized view that is use to run all downstream analysis.

~~~
apptainer exec --env PUBLIC_DATA_DIR=$HOME/Software/public results/containers/gwas2eqtl_pleiotropy.sif snakemake -p --cores all -s workflows/11snkfl_insrt_postgres.yml  --configfile config/11snkfl_insrt_postgres.yml
~~~

Based on the previous "colocpleio" materialized view, this workflow will run all the analysis and create figures.

The main figures in the paper are computed with SNP_PP_H4=0.5, but some comparison in the supplementary are run with parameters "SNP_PP_H4=0.25" and "SNP_PP_H4=0.75".

~~~
snakemake -p --cores all -s workflows/20snkfl_all.yml --configfile config/20snkfl_all_30snkfl_vep_fig3f.yml
~~~

This extra step is required as dedicated container is used for VEP (Variant effect predictor).

~~~
apptainer pull --name results/containers/vep.sif docker://ensemblorg/ensembl-vep

snakemake -p --cores all -s workflows/20snkfl_all.yml --configfile config/20snkfl_all_30snkfl_vep_fig3f.yml --use-singularity --singularity-args '-B results/data/vep-cache:/opt/vep/vep-cache'
~~~

First we need to modify this file: results/20241003/pval_5e-08/r2_0.1/kb_1000/window_1000000/75_50/plthtmp_disease_comorbidity_matrix.py/corr.svg
by hand in inkscape and modify the figure for the MS.

Copy it to ms/corr_inkscape.png.

Then run latex with another container.

~~~
cd ms
sudo apptainer build ../results/containers/latex.sif ../containers/latex.def
rm -f *.aux *.dvi *.log *.out; apptainer run ../results/containers/latex.sif texi2pdf ms00_fig_table_suppl.tex; rm -f *.aux *.dvi *.log *.out
~~~
