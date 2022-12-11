## Analyse pleiotropy using this repository.

First the database must be constructed

# docker postgres

Dev

~~~
cd container
docker compose --project-name gwas2eqtl_pleiotropy_dev -f docker-compose.yml up --build --force-recreate --remove-orphans -d
~~~

Prod

~~~
cd container
docker compose --project-name gwas2eqtl_pleiotropy_prod --env-file env_prod -f docker-compose.yml up --build --force-recreate --remove-orphans -d
~~~

Dev - 5433 - One chrom

~~~
snakemake -p -j all -s snkfl_all.yml --config postgres_port=5433 outdir=out/5433 public_data_dir=/home/gonzalez/Software/public process_data_dir=/home/gonzalez/Software/process CHR=22 --resources db=1
~~~

Dev - 5433 - All chrom

~~~
snakemake -p -j all -s snkfl_all.yml --config postgres_port=5433 outdir=out/5433 public_data_dir=/home/gonzalez/Software/public process_data_dir=/home/gonzalez/Software/process CHR=22 --resources db=1
~~~

Prod - 5434 - One chrom

~~~
snakemake -p -j all -s snkfl_all.yml --config postgres_port=5434 outdir=out/5434 public_data_dir=/home/gonzalez/Software/public process_data_dir=/home/gonzalez/Software/process CHR=22 --resources db=1
~~~

Prod - 5434 - All chrom

~~~
snakemake -p -j all -s snkfl_all.yml --config postgres_port=5434 outdir=out/5434 public_data_dir=/home/gonzalez/Software/public process_data_dir=/home/gonzalez/Software/process --resources db=1
~~~

Insert to db

~~~
python scripts/insrt_cytoband.py postgresql://chiliconcarne101usr:chiliconcarne101pwd@0.0.0.0:5435/gwas2eqtl
python scripts/insrt_geneid2symbol.py postgresql://chiliconcarne101usr:chiliconcarne101pwd@0.0.0.0:5435/gwas2eqtl
python scripts/insrt_gwas_annot.py postgresql://chiliconcarne101usr:chiliconcarne101pwd@0.0.0.0:5435/gwas2eqtl config/gwas418.ods
python scripts/insrt_pos19.py postgresql://chiliconcarne101usr:chiliconcarne101pwd@0.0.0.0:5435/gwas2eqtl
~~~

From the "gwas2eqtl" project

~~~
python scripts/insrt_coloc.py 0.7 0.7  postgresql://chiliconcarne101usr:chiliconcarne101pwd@0.0.0.0:5435/gwas2eqtl config/gwas418.ods /home/gonzalez/Software/public/raw.githubusercontent.com/eQTL-Catalogue/eQTL-Catalogue-resources/master/tabix/tabix_ftp_paths.tsv /home/gonzalez/Repositories/gwas2eqtl/out/gwas418/coloc/{gwas_id}/pval_5e-08/r2_0.1/kb_1000/window_1000000/{eqtl_id}.tsv
python workflow/scripts/insrt_tophits.py postgresql://chiliconcarne101usr:chiliconcarne101pwd@0.0.0.0:5435/gwas2eqtl config/gwas418.ods  /home/gonzalez/Repositories/gwas2eqtl/out/gwas418/tophits/{gwas_id}/pval_5e-08/r2_0.1/kb_1000/hg38.tsv
~~~

Create "colocweb" view

~~~
SELECT DISTINCT co.chrom,
    pos19.pos19,
    co.pos38,
    co.cytoband,
    co.rsid,
    co.ref,
    co.alt,
    gw.gwas_trait,
    gw.gwas_class,
    co.gwas_beta,
    en.symbol AS eqtl_gene_symbol,
    co.eqtl_beta,
    co.eqtl_id,
    co.eqtl_gene_id,
    co.gwas_id,
    co.gwas_pval,
    co.eqtl_pval,
    co.pp_h4_abf,
    co.snp_pp_h4,
    co.coloc_variant_id AS tophits_variant_id,
    co.nsnps
   FROM (((( SELECT DISTINCT co0.chrom,
            co0.pos AS pos38,
            concat_ws(''::text, co0.chrom, cy.cytoband) AS cytoband,
            concat('rs', co0.rsid) AS rsid,
            co0.ref,
            co0.alt,
            co0.gwas_beta,
            co0.eqtl_beta,
            co0.eqtl_id,
            co0.eqtl_gene_id,
            co0.gwas_id,
            co0.gwas_pval,
            co0.eqtl_pval,
            co0.pp_h4_abf,
            co0.snp_pp_h4,
            co0.coloc_variant_id,
            co0.nsnps
           FROM coloc co0,
            cytoband cy
          WHERE ((co0.chrom = cy.chrom) AND (co0.pos <@ cy.start_end38))) co
     LEFT JOIN ensg2symbol en ON (((en.gene_id)::text = (co.eqtl_gene_id)::text)))
     LEFT JOIN gwas_annot gw ON (((gw.gwas_id)::text = (co.gwas_id)::text)))
     LEFT JOIN pos19 ON ((co.pos38 = pos19.pos)))
  ORDER BY co.chrom, pos19.pos19, co.pos38, co.alt, gw.gwas_trait, en.symbol, co.eqtl_id;
~~~

Create "colocpleio" view

~~~
SELECT DISTINCT co.chrom,
    pos19.pos19,
    co.pos38,
    co.cytoband,
    co.rsid,
    co.ref,
    co.alt,
    gw.gwas_trait,
    gw.gwas_class,
    co.gwas_beta,
    en.symbol AS eqtl_gene_symbol,
    co.eqtl_beta,
    co.eqtl_id,
    co.eqtl_gene_id,
    co.gwas_id,
    co.gwas_pval,
    co.eqtl_pval,
    co.pp_h4_abf,
    co.snp_pp_h4,
    co.coloc_variant_id AS tophits_variant_id,
    co.nsnps,
    eqtl_annot.etissue_class
   FROM (((( SELECT DISTINCT co0.chrom,
            co0.pos AS pos38,
            concat_ws(''::text, co0.chrom, cy.cytoband) AS cytoband,
            concat('rs', co0.rsid) AS rsid,
            co0.ref,
            co0.alt,
            co0.gwas_beta,
            co0.eqtl_beta,
            co0.eqtl_id,
            co0.eqtl_gene_id,
            co0.gwas_id,
            co0.gwas_pval,
            co0.eqtl_pval,
            co0.pp_h4_abf,
            co0.snp_pp_h4,
            co0.coloc_variant_id,
            co0.nsnps
           FROM coloc co0,
            cytoband cy
          WHERE ((co0.chrom = cy.chrom) AND (co0.pos <@ cy.start_end38))) co
     LEFT JOIN ensg2symbol en ON (((en.gene_id)::text = (co.eqtl_gene_id)::text)))
     LEFT JOIN gwas_annot gw ON (((gw.gwas_id)::text = (co.gwas_id)::text)))
     LEFT JOIN pos19 ON ((co.pos38 = pos19.pos))
     LEFT JOIN eqtl_annot ON ((co.eqtl_id = eqtl_annot.eqtl_id))
)
  ORDER BY co.chrom, pos19.pos19, co.pos38, co.alt, gw.gwas_trait, en.symbol, co.eqtl_id
~~~

Then snakemake is run with:

~~~
snakemake -j all -s tools/00snkfl_all.yml --config david_email=${DAVID_EMAIL} db_url=postgresql://chiliconcarne101usr:chiliconcarne101pwd@0.0.0.0:5435/gwas2eqtl etissue_class_ods=config/etissue_class.ods gwas_class_ods=config/gwas418.ods max_gwas_class_count=4 outdir=out/gwas418/pval_5e-08/r2_0.1/kb_1000/window_1000000/7_7 public_data_dir=/home/gonzalez/Software/publicsnakemake -j all -s tools/00snkfl_all.yml --config david_email=${DAVID_EMAIL} db_url=postgresql://chiliconcarne101usr:chiliconcarne101pwd@0.0.0.0:5435/gwas2eqtl etissue_class_ods=config/etissue_class.ods gwas_class_ods=config/gwas418.ods max_gwas_class_count=4 outdir=out/gwas418/pval_5e-08/r2_0.1/kb_1000/window_1000000/7_7 public_data_dir=/home/gonzalez/Software/public
~~~

~~~
export david_email=aitor.gonzalez@inserm.fr; snakemake --cores all -p -d ${PWD} -s tools/snkfl_vep.yml --config coloc_tsv_gz=../gwas2eqtl/out/gwas420/coloc_gwas418.tsv.gz outdir=out/gwas418/pval_5e-08/r2_0.1/kb_1000/window_1000000 etissue_class_ods=config/etissue_class.ods max_gwas_class_count=5 gwas_class_ods=config/gwas418.ods public_data_dir=${HOME}/Software/public david_email=${david_email}
~~~

# MS

The MS template is taken from here: <https://github.com/quantixed/manuscript-templates>

Compilation works like this:

~~~
rm -f ms_00_merge.tex.pdf; singularity exec -u ../singularity/out/latex.sif texi2pdf ms/ms00_merge.tex; rm -f *.aux *.dvi *.log *.out
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
rm -f ms_00_merge.tex.pdf; texi2pdf ms/ms00_merge.tex; rm -f *.aux *.dvi *.log *.out
~~~

# Poster

Go to poster_eccb22_barcelona and follow the readme.md instructions

python scripts/plthtmp_disease_comorbidity_matrix.py out/gwas420/pval_5e-08/r2_0.1/kb_1000/window_1000000/cmpt_disease_comorbidity_matrix2.py/corr3.tsv config/gwas418.ods out/gwas420/pval_5e-08/r2_0.1/kb_1000/window_1000000/plthtmp_disease_comorbidity_matrix.py/corr3.png

# Presentation GOLD meeting 2022

~~~
texi2pdf presentation_gold2022_paris/presentation_gold2022_paris.tex
~~~