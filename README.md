# Run analyses and generate MS figures

The input path is given as "raw_coloc_tsv" and it must be copied somewhere in this folder

~~~
rsync -avzP gonzalez@warsaw:/home/gonzalez/Repositories/eqtl2gwas/out/merged/genome/5e-08/1000000/coloc413.tsv /home/gonzalez/Repositories/eqtl2gwas_pleiotropy/out/gwas413/genome/5e-08/1000000/.
~~~

The columns are these:

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

Then snamake is run with:

~~~
snakemake --cores all -p -d ${PWD} -s tools/Snakefile.yml --config outdir=out/gwas413/genome/5e-08/1000000 raw_coloc_tsv=out/gwas413/genome/5e-08/1000000/coloc413.tsv gwas_cat_ods=config/gwas413.ods etissue_cat_ods=config/etissue_category.ods upper_var_gwas_cat_count=5 public_data_dir=/home/gonzalez/Software/public
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
