# eqtl2gwas_pleiotropy

## Colocalization using eqtl2gwas

The input data was generated using https://github.com/aitgon/eqtl2gwas, tag, 0.1.1), run with pval=5e-8, window 1000000

~~~
git clone git@github.com:aitgon/eqtl2gwas.git
cd eqtl2gwas
PYTHONPATH=.:$PYTHONPATH snakemake -j all -s workflow/Snakefile_eqtl.yml -p --rerun-incomplete  --config  public_data_dir=$HOME/Software/public process_data_dir=$HOME/Software/process region=genome eqtl_fdr=0.05 window=500000
PYTHONPATH=.:$PYTHONPATH snakemake -j all -s workflow/Snakefile_gwas.yml -p --config gwas_ods=gwas413.ods gwas_pval=5e-8 public_data_dir=/home/gonzalez/Software/public process_data_dir=/home/gonzalez/Software/process image_sif=out/eqt2gwas.sif --rerun-incomplete
snakemake -c all -s workflow/Snakefile.yml -p --config gwas_ods=gwas413.ods gwas_pval=5e-8 public_data_dir=/home/gonzalez/Software/public process_data_dir=/home/gonzalez/Software/process region=genome window=1000000 eqtl_fdr=0.05 image_sif=out/eqtl2gwas.sif --rerun-incomplete
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
snakemake --cores all -p -d ${PWD} -s tools/Snakefile.yml --config outdir=out/gwas413/genome/5e-08/1000000 coloc_tsv_gz=out/gwas413/genome/5e-08/1000000/coloc413.tsv.gz gwas_cat_ods=config/gwas413.ods etissue_cat_ods=config/etissue_category.ods upper_var_gwas_cat_count=5 public_data_dir=/home/gonzalez/Software/public
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
snakemake --cores all -p -d ${PWD} -s tools/Snakefile.yml --config outdir=out/gwas413/genome/5e-08/1000000 raw_coloc_tsv=config/coloc413.tsv gwas_cat_ods=config/gwas413.ods etissue_cat_ods=config/etissue_category.ods upper_var_gwas_cat_count=5 public_data_dir=/home/gonzalez/Software/public
~~~

# MS

The MS template is taken from here: <https://github.com/quantixed/manuscript-templates>

Compilation works like this:

~~~
rm -f ms_00_merge.tex.pdf; singularity exec -u ../singularity/out/latex.sif texi2pdf ms/ms00_merge.tex; rm -f *.aux *.dvi *.log *.out
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

# MS Supplementary table

~~~
python scripts/create_supplementary_table.py out/gwas413/genome/5e-08/1000000 out/gwas413/genome/5e-08/1000000/ms_supp_tab/supplementary_table.xlsx
~~~
