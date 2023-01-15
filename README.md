Follow the instructions in "doc/prod_montreal.md"

## Analyse pleiotropy using this repository.

First the database must be constructed

# docker postgres

Dev

~~~
cd container
docker compose --project-name gwas2eqtl_pleiotropy_dev -f docker-compose.yml up --build --force-recreate --remove-orphans -d
~~~

# Run analyses and generate MS figures

# MS

Texlive dependencies

~~~
tlmgr install orcidlink
tlmgr install lipsum
tlmgr install preprint
tlmgr install csvsimple
~~~

Compilation works like this:

~~~
rm -f ms_00_merge.tex.pdf; texi2pdf ms00_merge.tex; rm -f *.aux *.dvi *.log *.out
~~~

# Poster

Go to poster_eccb22_barcelona and follow the readme.md instructions

python scripts/plthtmp_disease_comorbidity_matrix.py out/gwas420/pval_5e-08/r2_0.1/kb_1000/window_1000000/cmpt_disease_comorbidity_matrix2.py/corr3.tsv config/gwas418.ods out/gwas420/pval_5e-08/r2_0.1/kb_1000/window_1000000/plthtmp_disease_comorbidity_matrix.py/corr3.png

# Presentation GOLD meeting 2022

~~~
cd presentation_230120_gold2022_paris
texi2pdf presentation_gold2022_paris.tex
~~~

# Presentation Internal Seminar 20 jan 2023

~~~
cd presentation_230120_internal_seminar
texi2pdf intern_sem.tex
~~~
