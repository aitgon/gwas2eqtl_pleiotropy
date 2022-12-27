Updated Dec 27, 2022

Annotate GWAS with disease ontologies

~~~
python scripts/query_ontology_ols.py config/gwas417.ods out/gwas417/pval_5e-08/r2_0.1/kb_1000/window_1000000/75_50/query_ontology.py/gwas_ontology.ods
~~~

Docker up

~~~
cd container
docker compose --project-name gwas2eqtl_pleiotropy_dev -f docker-compose.yml up --build --force-recreate --remove-orphans -d
cd ..
~~~

From the "gwas2eqtl" project

~~~
python workflow/scripts/insrt_coloc.py 0.75 0  postgresql://postgres:postgres@0.0.0.0:5435/postgres ../gwas2eqtl_pleiotropy/config/gwas417.ods /home/gonzalez/Software/public/raw.githubusercontent.com/eQTL-Catalogue/eQTL-Catalogue-resources/master/tabix/tabix_ftp_paths.tsv /home/gonzalez/Repositories/gwas2eqtl/out/gwas418/coloc/{gwas_id}/pval_5e-08/r2_0.1/kb_1000/window_1000000/{eqtl_id}.tsv
python workflow/scripts/insrt_tophits.py postgresql://postgres:postgres@0.0.0.0:5435/postgres ../gwas2eqtl_pleiotropy/config/gwas417.ods  /home/gonzalez/Repositories/gwas2eqtl/out/gwas418/tophits/{gwas_id}/pval_5e-08/r2_0.1/kb_1000/hg38.tsv
~~~

python scripts/insrt_vep_consequence.py postgresql://postgres:postgres@0.0.0.0:5435/postgres

Insert to db

~~~
python scripts/insrt_cytoband.py postgresql://postgres:postgres@0.0.0.0:5435/postgres
python scripts/insrt_geneid2symbol.py postgresql://postgres:postgres@0.0.0.0:5435/postgres
python scripts/insrt_gwas_annot.py postgresql://postgres:postgres@0.0.0.0:5435/postgres config/gwas417.ods
python scripts/insrt_pos19.py postgresql://postgres:postgres@0.0.0.0:5435/postgres
python scripts/insrt_etissue_category.py postgresql://postgres:postgres@0.0.0.0:5435/postgres config/etissue_category.ods
python scripts/insrt_open_gwas.py postgresql://postgres:postgres@0.0.0.0:5435/postgres
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
    gw.gwas_category,
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
    gw.gwas_category,
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
    eqtl_annot.etissue_category_term,
    op.pmid
   FROM (((((( SELECT DISTINCT co0.chrom,
            co0.pos AS pos38,
            concat_ws(''::text, co0.chrom, cy.cytoband) AS cytoband,
            co0.rsid,
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
     LEFT JOIN open_gwas_info op ON (((op.gwas_id)::text = (co.gwas_id)::text)))
     LEFT JOIN pos19 ON ((co.pos38 = pos19.pos)))
     LEFT JOIN eqtl_annot ON (((co.eqtl_id)::text = (eqtl_annot.eqtl_id)::text)))
  ORDER BY co.chrom, pos19.pos19, co.pos38, co.alt, gw.gwas_trait, en.symbol, co.eqtl_id;
~~~

Optional for the annotation

~~~
python scripts/query_ontology_ols.py config/gwas417.ods out/gwas418/pval_5e-08/r2_0.1/kb_1000/window_1000000/75_50/query_ontology_ols.py/gwas_ontology.ods
~~~

Then snakemake is run with:

~~~
snakemake -j all -s tools/00snkfl_all.yml --config david_email=${DAVID_EMAIL} db_url=postgresql://postgres:postgres@0.0.0.0:5435/postgres etissue_class_ods=config/etissue_class.ods gwas_class_ods=config/gwas417.ods max_gwas_class_count=4 outdir=out/gwas418/pval_5e-08/r2_0.1/kb_1000/window_1000000/75_50 public_data_dir=/home/gonzalez/Software/public snp_pp_h4=0.5 -p 
~~~

~~~
snakemake --cores all -p -d ${PWD} -s tools/snkfl_vep.yml --config outdir=out/gwas418/pval_5e-08/r2_0.1/kb_1000/window_1000000/75_50 max_gwas_class_count=4 db_url=postgresql://postgres:postgres@0.0.0.0:5435/postgres snp_pp_h4=0.5 -p
~~~

# MS

~~~
rm -f ms_00_merge.tex.pdf; texi2pdf ms/ms00_merge.tex; rm -f *.aux *.dvi *.log *.out
~~~

python scripts/query_ontology_ols.py config/gwas417.ods out/gwas418/pval_5e-08/r2_0.1/kb_1000/window_1000000/75_50/annot_gwas_metadata.py/gwas_ontology.ods
