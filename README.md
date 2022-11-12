# gwas2eqtl_pleiotropy

## Colocalization using gwas2eqtl

The input data was generated using https://github.com/aitgon/gwas2eqtl, tag, 0.1.1), run with pval=5e-8, window 1000000

~~~
git clone git@github.com:aitgon/gwas2eqtl.git
cd gwas2eqtl
PYTHONPATH=.:$PYTHONPATH snakemake -j all -s workflow/Snakefile_eqtl.yml -p --rerun-incomplete  --config  public_data_dir=$HOME/Software/public process_data_dir=$HOME/Software/process region=genome eqtl_fdr=0.05 window=500000
PYTHONPATH=.:$PYTHONPATH snakemake -j all -s workflow/Snakefile_gwas.yml -p --config gwas_ods=gwas413.ods gwas_pval=5e-8 public_data_dir=/home/gonzalez/Software/public process_data_dir=/home/gonzalez/Software/process image_sif=out/eqt2gwas.sif --rerun-incomplete
snakemake -c all -s workflow/Snakefile.yml -p --config gwas_ods=gwas413.ods gwas_pval=5e-8 public_data_dir=/home/gonzalez/Software/public process_data_dir=/home/gonzalez/Software/process region=genome window=1000000 eqtl_fdr=0.05 image_sif=out/gwas2eqtl.sif --rerun-incomplete
~~~

The previous commands result in this file coloc413.tsv.gz that can be also foud in the OSF site

These are the columns:

- chrom
- pos
- rsid
- ref
- alt
- egene
- egene_symbol
- eqtl_beta
- eqtl_pvalue
- eqtl_identifier
- gwas_beta
- gwas_pvalue
- gwas_identifierpp_h4
- PP.H4.abf
- coloc_window
- nsnps
- PP.H3.abf
- PP.H2.abf
- PP.H1.abf
- PP.H0.abf

## Analyse pleiotropy using this repository.

Then snamake is run with:

~~~
snakemake --cores all -p -d ${PWD} -s tools/Snakefile.yml --config outdir=out/gwas413/genome/5e-08/1000000 coloc_tsv_gz=out/gwas413/genome/5e-08/1000000/coloc413.tsv.gz gwas_cat_ods=config/gwas413.ods etissue_cat_ods=config/etissue_class.ods upper_var_gwas_cat_count=5 public_data_dir=/home/gonzalez/Software/public
~~~

# MS

The MS template is taken from here: <https://github.com/quantixed/manuscript-templates>

Compilation works like this:

~~~
rm -f ms_00_merge.tex.pdf; singularity exec -u ../singularity/out/latex.sif texi2pdf ms/ms00_merge.tex; rm -f *.aux *.dvi *.log *.out
~~~

Copy input data to the config folder; eg out/config/genome/5e-08/1000000

# Run analyses and generate MS figures

~~~
export david_email=aitor.gonzalez@inserm.fr; snakemake --cores all -p -d ${PWD} -s tools/00snkfl_all.yml --config coloc_tsv_gz=../gwas2eqtl/out/gwas420/coloc_gwas418.tsv.gz outdir=out/gwas418/pval_5e-08/r2_0.1/kb_1000/window_1000000 etissue_class_ods=config/etissue_class.ods max_gwas_class_count=5 gwas_class_ods=config/gwas418.ods public_data_dir=${HOME}/Software/public david_email=${david_email}
~~~

~~~
export david_email=aitor.gonzalez@inserm.fr; snakemake --cores all -p -d ${PWD} -s tools/snkfl_vep.yml --config coloc_tsv_gz=../gwas2eqtl/out/gwas420/coloc_gwas418.tsv.gz outdir=out/gwas418/pval_5e-08/r2_0.1/kb_1000/window_1000000 etissue_class_ods=config/etissue_class.ods max_gwas_class_count=5 gwas_class_ods=config/gwas418.ods public_data_dir=${HOME}/Software/public david_email=${david_email}
~~~

# MS

The MS template is taken from here: <https://github.com/quantixed/manuscript-templates>

Compilation works like this:

~~~
rm -f ms_00_merge.tex.pdf; texi2pdf ms/ms00_merge.tex; rm -f *.aux *.dvi *.log *.out
~~~

Download dependencies

~~~
wget -nc http://mirrors.ctan.org/fonts/ifsym/ifsym.sty -P ms
wget -nc http://mirrors.ctan.org/macros/latex/contrib/chemgreek/chemgreek.sty -P ms
wget -nc http://mirrors.ctan.org/macros/latex/contrib/mhchem/mhchem.sty -P ms
wget -nc http://mirrors.ctan.org/macros/latex/contrib/siunitx/siunitx-abbreviations.cfg -P ms
wget -nc http://mirrors.ctan.org/macros/latex/contrib/siunitx/siunitx-v2.sty -O ms/siunitx.sty
wget -nc https://raw.githubusercontent.com/ldbc/ldbc_graphalytics_docs/master/bbding.sty -P ms
wget -nc https://raw.githubusercontent.com/quantixed/manuscript-templates/master/bioRxiv.cls -P ms
wget -nc https://raw.githubusercontent.com/quantixed/manuscript-templates/master/orcidlink.sty -P ms
wget -nc https://raw.githubusercontent.com/quantixed/manuscript-templates/master/bioRxiv_logo.png -P ms
convert ms/bioRxiv_logo.png ms/bioRxiv_logo.eps

wget -nc http://mirrors.ctan.org/macros/latex/contrib/textgreek/textgreek.ins -P ms
wget -nc http://mirrors.ctan.org/macros/latex/contrib/textgreek/textgreek.dtx -P ms
cd ms
latex textgreek.ins
cd $OLDPWD
~~~

# Poster

Go to poster_eccb22_barcelona and follow the readme.md instructions

python scripts/plthtmp_disease_comorbidity_matrix.py out/gwas420/pval_5e-08/r2_0.1/kb_1000/window_1000000/cmpt_disease_comorbidity_matrix2.py/corr3.tsv config/gwas418.ods out/gwas420/pval_5e-08/r2_0.1/kb_1000/window_1000000/plthtmp_disease_comorbidity_matrix.py/corr3.png

# Presentation GOLD meeting 2022

~~~
texi2pdf presentation_gold2022_paris/presentation_gold2022_paris.tex
~~~


