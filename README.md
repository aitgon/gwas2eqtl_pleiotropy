Copy input data to the config folder; eg out/config/genome/5e-08/1000000

# Run analyses and generate MS figures

~~~
snakemake --cores 1 -p -d ${PWD} -s tools/Snakefile2.yml --config outdir=out/gwas413/genome/5e-08/1000000 raw_coloc_tsv=config/coloc413.tsv gwas_cat_ods=config/gwas413.ods etissue_cat_ods=config/etissue_category.ods upper_var_gwas_cat_count=5 public_data_dir=/home/gonzalez/Software/public
~~~

# MS

The MS template is taken from here: <https://github.com/quantixed/manuscript-templates>

Compilation works like this:

~~~
rm -f ms_00_merge.tex.pdf; singularity exec -u ../singularity/out/latex.sif texi2pdf ms_00_merge.tex; rm -f *.aux *.dvi *.log
~~~

Download dependencies

~~~
wget -nc http://mirrors.ctan.org/macros/latex/contrib/siunitx/siunitx-v2.sty -O siunitx.sty
wget -nc http://mirrors.ctan.org/macros/latex/contrib/mhchem/mhchem.sty
wget -nc http://mirrors.ctan.org/macros/latex/contrib/chemgreek/chemgreek.sty
wget -nc http://mirrors.ctan.org/fonts/ifsym/ifsym.sty
wget -nc https://raw.githubusercontent.com/ldbc/ldbc_graphalytics_docs/master/bbding.sty
wget -nc http://mirrors.ctan.org/macros/latex/contrib/siunitx/siunitx-abbreviations.cfg
wget -nc https://raw.githubusercontent.com/quantixed/manuscript-templates/master/bioRxiv.cls
wget -nc http://mirrors.ctan.org/macros/latex/contrib/textgreek/textgreek.ins
wget -nc http://mirrors.ctan.org/macros/latex/contrib/textgreek/textgreek.dtx
latex textgreek.ins
wget -nc https://raw.githubusercontent.com/quantixed/manuscript-templates/master/orcidlink.sty
~~~

