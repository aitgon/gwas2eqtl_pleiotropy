Updated Feb 5, 2023

Install

~~~
~~~

Annotate GWAS trait with disease trait and category ontologies

~~~
python scripts/query_ontology_ols.py config/gwas417_trait_query.ods out/gwas417/query_ontology.py/gwas_trait_ontology.ods
python scripts/query_ontology_ols.py config/gwas417_category_query.ods out/gwas417/query_ontology.py/gwas_category_ontology.ods
~~~

Docker up

~~~
cd container
docker compose --project-name gwas2eqtl_pleiotropy_dev -f docker-compose.yml up --build --force-recreate --remove-orphans -d
cd ..
~~~

From the "gwas2eqtl" project

~~~
python workflow/scripts/insrt_coloc.py 0.75 0  postgresql://postgres:postgres@0.0.0.0:5435/postgres ../gwas2eqtl_pleiotropy/config/gwas417_query_category.ods /home/gonzalez/Software/public/raw.githubusercontent.com/eQTL-Catalogue/eQTL-Catalogue-resources/master/tabix/tabix_ftp_paths.tsv /home/gonzalez/Repositories/gwas2eqtl/out/gwas417/coloc/{gwas_id}/pval_5e-08/r2_0.1/kb_1000/window_1000000/{eqtl_id}.tsv
python workflow/scripts/insrt_tophits.py postgresql://postgres:postgres@0.0.0.0:5435/postgres ../gwas2eqtl_pleiotropy/config/gwas_trait_ontology.ods  /home/gonzalez/Repositories/gwas2eqtl/out/gwas418/tophits/{gwas_id}/pval_5e-08/r2_0.1/kb_1000/hg38.tsv
~~~

Insert GWAS trait and category ontologies

~~~
python scripts/insrt_opengwas2trait_ontology.py postgresql://postgres:postgres@0.0.0.0:5435/postgres out/gwas417/query_ontology.py/gwas_trait_ontology.ods
python scripts/insrt_opengwas2category_ontology.py postgresql://postgres:postgres@0.0.0.0:5435/postgres out/gwas417/query_ontology.py/gwas_category_ontology.ods
~~~

python scripts/insrt_vep_consequence.py postgresql://postgres:postgres@0.0.0.0:5435/postgres

Insert to db

~~~
python scripts/insrt_cytoband.py postgresql://postgres:postgres@0.0.0.0:5435/postgres
python scripts/insrt_geneid2symbol.py postgresql://postgres:postgres@0.0.0.0:5435/postgres
#python scripts/insrt_gwas_annot.py postgresql://postgres:postgres@0.0.0.0:5435/postgres config/gwas417_query_precise.ods
python scripts/insrt_pos19.py postgresql://postgres:postgres@0.0.0.0:5435/postgres
python scripts/insrt_etissue_category.py postgresql://postgres:postgres@0.0.0.0:5435/postgres config/etissue_category.ods
python scripts/insrt_open_gwas.py postgresql://postgres:postgres@0.0.0.0:5435/postgres
python scripts/insrt_entrezgene2ensg2symbol.py postgresql://postgres:postgres@0.0.0.0:5435/postgres
python scripts/insrt_entrezgene2pubmed_count.py postgresql://postgres:postgres@0.0.0.0:5435/postgres
insrt_ncbirefseq -u postgresql://postgres:postgres@0.0.0.0:5435/postgres
~~~

Pleiotropy from Watanabe et al. 2019

~~~
wget -r -c https://static-content.springer.com/esm/art%3A10.1038%2Fs41588-019-0481-0/MediaObjects/41588_2019_481_MOESM3_ESM.xlsx -P ~/Software/public
insrt_watanabe_posthuma2019 -u postgresql://postgres:postgres@0.0.0.0:5435/postgres -i ~/Software/public/static-content.springer.com/esm/art:10.1038%2Fs41588-019-0481-0/MediaObjects/41588_2019_481_MOESM3_ESM.xlsx
~~~

Create "ensg2pubmed_count" view

~~~
SELECT entrezgene2ensg2symbol.gene_id,
    entrezgene2pubmed_count.pubmed_count
   FROM entrezgene2ensg2symbol,
    entrezgene2pubmed_count
  WHERE (entrezgene2ensg2symbol.entrezgene = entrezgene2pubmed_count.entrezgene);
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
    gwtron.gwas_trait,
    gwtron.gwas_ontology_term AS gwas_trait_ontology_term,
    gwtron.gwas_ontology_id AS gwas_trait_ontology_id,
    gwcaon.gwas_ontology_term AS gwas_category_ontology_term,
    gwcaon.gwas_ontology_id AS gwas_category_ontology_id,
    op.batch,
    op.pmid,
    co.gwas_beta,
    refseq.symbol AS eqtl_gene_symbol,
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
    en2.pubmed_count,
    af.eas_af,
    af.amr_af,
    af.afr_af,
    af.eur_af,
    af.sas_af,
    refseq.refseq_transcript_id AS eqtl_refseq_transcript_id,
    refseq.refseq_transcript_start38 AS eqtl_refseq_transcript_start38,
    refseq.refseq_transcript_end38 AS eqtl_refseq_transcript_end38,
    refseq.refseq_transcript_strand AS eqtl_refseq_transcript_strand,
    wp.domains AS domains_watanabe2019
   FROM (((((((((( SELECT DISTINCT co0.chrom,
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
     LEFT JOIN "genome-mysql.soe.ucsc.edu/hg38/ncbirefseq" refseq ON (((refseq.gene_id)::text = (co.eqtl_gene_id)::text)))
     LEFT JOIN opengwas2trait_ontology gwtron ON (((gwtron.gwas_id)::text = (co.gwas_id)::text)))
     LEFT JOIN opengwas2category_ontology gwcaon ON (((gwcaon.gwas_id)::text = (co.gwas_id)::text)))
     LEFT JOIN open_gwas_info op ON (((op.gwas_id)::text = (co.gwas_id)::text)))
     LEFT JOIN pos19 ON ((co.pos38 = pos19.pos)))
     LEFT JOIN eqtl_annot ON (((co.eqtl_id)::text = (eqtl_annot.eqtl_id)::text)))
     LEFT JOIN ensg2pubmed_count en2 ON (((co.eqtl_gene_id)::text = (en2.gene_id)::text)))
     LEFT JOIN af_1000genomes af ON ((((co.rsid)::text = (af.rsid)::text) AND ((co.ref)::text = (af.ref)::text) AND ((co.alt)::text = (af.alt)::text))))
     LEFT JOIN watanabe_posthuma2019 wp ON ((((co.rsid)::text = (wp.rsid)::text) AND ((co.ref)::text = (wp.ref)::text) AND ((co.alt)::text = (wp.alt)::text))))
  ORDER BY co.chrom, pos19.pos19, co.pos38, co.alt, gwtron.gwas_trait, refseq.symbol, co.eqtl_id;
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
    gwtron.gwas_trait,
    gwcaon.gwas_ontology_term as gwas_class,
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
   FROM (((((((( SELECT DISTINCT co0.chrom,
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
     LEFT JOIN opengwas2trait_ontology gwtron ON (((gwtron.gwas_id)::text = (co.gwas_id)::text)))
     LEFT JOIN opengwas2category_ontology gwcaon ON (((gwcaon.gwas_id)::text = (co.gwas_id)::text)))
     LEFT JOIN open_gwas_info op ON (((op.gwas_id)::text = (co.gwas_id)::text)))
     LEFT JOIN pos19 ON ((co.pos38 = pos19.pos)))
     LEFT JOIN eqtl_annot ON (((co.eqtl_id)::text = (eqtl_annot.eqtl_id)::text)))
     LEFT JOIN ensg2pubmed_count en2 ON (((co.eqtl_gene_id)::text = (en2.gene_id)::text)))
  ORDER BY co.chrom, pos19.pos19, co.pos38, co.alt, gwtron.gwas_trait, en.symbol, co.eqtl_id;
~~~

Then snakemake is run with:

snp_pp_h4 0.25 0.5 0.75

~~~
export DAVID_EMAIL=myemaildavid@example.com
export SNP_PP_H4=0.5
export OUTDIR=out/gwas417/pval_5e-08/r2_0.1/kb_1000/window_1000000/75_50
~~~

~~~
snakemake -p --cores all -s tools/00snkfl_all.yml --config david_email=${DAVID_EMAIL} db_url=postgresql://postgres:postgres@0.0.0.0:5435/postgres etissue_category_ods=config/etissue_category.ods gwas_trait_ods=out/gwas417/query_ontology.py/gwas_trait_ontology.ods  gwas_category_ods=out/gwas417/query_ontology.py/gwas_category_ontology.ods max_gwas_category_count=4 outdir="${OUTDIR}" public_data_dir="${HOME}"/Software/public snp_pp_h4="${SNP_PP_H4}"

snakemake -p --cores all -s tools/snkfl_vep.yml --config db_url=postgresql://postgres:postgres@0.0.0.0:5435/postgres max_gwas_category_count=4 outdir="${OUTDIR}" public_data_dir=/home/gonzalez/Software/public process_data_dir=/home/gonzalez/Software/process snp_pp_h4="${SNP_PP_H4}"

snakemake -p --cores all -s tools/06_ucsc_hubs.yml --config public_data_dir=/home/gonzalez/Software/public db_url=postgresql://postgres:postgres@0.0.0.0:5435/postgres outdir="${OUTDIR}" snp_pp_h4="${SNP_PP_H4}"
~~~

UCSC track hub

~~~
snakemake -p --cores 1 -s tools/06_ucsc_hubs.yml --config outdir=out/gwas417/pval_5e-08/r2_0.1/kb_1000/window_1000000/75_50 public_data_dir=/home/gonzalez/Software/public snp_pp_h4=0.50 db_url=postgresql://postgres:postgres@0.0.0.0:5435/postgres
python scripts/cmpt_track_hub.py out/gwas417/pval_5e-08/r2_0.1/kb_1000/window_1000000/75_50/cmpt_ucsc_tracks_eqtl_inter_bed.py/bigbed out/gwas417/pval_5e-08/r2_0.1/kb_1000/window_1000000/75_50/cmpt_count_per_rsid.py/count_per_rsid_gwas_egene_etissue.ods  ../gwas2eqtl_trackhub
~~~

# MS

~~~
cd ms
rm -f *.aux *.dvi *.log *.out; texi2pdf ms00_bmc_article.tex; rm -f *.aux *.dvi *.log *.out
rm -f *.aux *.dvi *.log *.out; texi2pdf ms00_fig_suppl.tex; rm -f *.aux *.dvi *.log *.out
rm -f *.aux *.dvi *.log *.out; texi2pdf ms00_fig_table_suppl.tex; rm -f *.aux *.dvi *.log *.out
~~~
