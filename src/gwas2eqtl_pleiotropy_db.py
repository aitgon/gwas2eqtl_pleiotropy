from sqlalchemy import Column, Integer, String, SmallInteger, Float, UniqueConstraint
from sqlalchemy.orm import declarative_base
from sqlalchemy.dialects.postgresql import INT4RANGE

# declarative base class
Base = declarative_base()

class tophits(Base):
   """scripts/tophits2db2.py"""
   __tablename__ = "tophits"
   __table_args__ = (UniqueConstraint('chrom', 'pos', 'ea', 'gwas_id', name='_chrom_pos_uc'),)

   id = Column('id', SmallInteger, primary_key=True)
   chrom = Column('chrom', SmallInteger, nullable=False)
   pos = Column('pos', Integer, nullable=False)
   rsid= Column('rsid', Integer, nullable=False)
   nea = Column('nea', String(255), nullable=False)
   ea = Column('ea', String(255), nullable=False)
   pval = Column('pval', Float, nullable=False)  # comma sep list of pvals
   beta = Column('beta', Float, nullable=False)  # comma sep list of betas
   se = Column('se', Float, nullable=False)
   eaf = Column('eaf', Float, nullable=True)
   n = Column('n', Integer, nullable=True)
   gwas_id = Column('gwas_id', String(63), primary_key=True)
   pos19 = Column('pos19', Integer, nullable=False)
   Column('coloc_lead_pos', Integer, nullable=False)
   Column('coloc_variant_id', String(50), nullable=False)
   Column('coloc_region', String(50), nullable=False)

class coloc(Base):
   """scripts/insrt_coloc.py"""
   __tablename__ = "coloc"
   __table_args__ = (UniqueConstraint('chrom', 'pos', 'alt', 'eqtl_gene_id', 'gwas_id', 'eqtl_id', name='_coloc_uc'),)

   id = Column('id', Integer, primary_key=True)
   chrom = Column('chrom', SmallInteger, nullable=False)
   pos = Column('pos', Integer, nullable=False)
   rsid= Column('rsid', Integer, nullable=False)
   ref = Column('ref', String(127), nullable=False)
   alt = Column('alt', String(127), nullable=False)
   pval = Column('gwas_pval', Float, nullable=False)
   beta = Column('gwas_beta', Float, nullable=False)
   eqtl_gene_id = Column('eqtl_gene_id', String(15), nullable=False)
   gwas_id = Column('gwas_id', String(127), nullable=False)
   eqtl_pval = Column('eqtl_pval', Float, nullable=False)
   eqtl_beta = Column('eqtl_beta', Float, nullable=False)
   eqtl_id = Column('eqtl_id', String(127), nullable=False)
   pp_h4_abf = Column('pp_h4_abf', Float, nullable=False)
   snp_pp_h4 = Column('snp_pp_h4', Float, nullable=False)
   nsnps = Column('nsnps', SmallInteger, nullable=False)
   coloc_variant_id = Column('coloc_variant_id', String(63), nullable=False)
   coloc_region = Column('coloc_region', String(63), nullable=False)
   pp_h3_abf = Column('pp_h3_abf', Float, nullable=False)
   pp_h2_abf = Column('pp_h2_abf', Float, nullable=False)
   pp_h1_abf = Column('pp_h1_abf', Float, nullable=False)
   pp_h0_abf = Column('pp_h0_abf', Float, nullable=False)

class pos19(Base):
   """src/insrt_pos19.py"""
   __tablename__ = "pos19"
   __table_args__ = (UniqueConstraint('chrom', 'pos', name='_pos19_uc'),)

   id = Column('id', Integer, primary_key=True)
   chrom = Column('chrom', SmallInteger, nullable=False)
   pos = Column('pos', Integer, nullable=False)
   pos19 = Column('pos19', Integer, nullable=True)


class opengwas2trait_ontology(Base):
   """src/insrt_opengwas2trait_ontology.py"""
   __tablename__ = "opengwas2trait_ontology"

   gwas_id = Column('gwas_id', String(63), primary_key=True)
   gwas_trait = Column('gwas_trait', String(255), nullable=False)
   gwas_ontology_term = Column('gwas_ontology_term', String(127), nullable=False)
   gwas_ontology_id = Column('gwas_ontology_id', String(31), nullable=False)
   gwas_ontology_iri = Column('gwas_ontology_iri', String(127), nullable=False)


class opengwas2category_ontology(Base):
   """src/insrt_opengwas2category_ontology.py"""
   __tablename__ = "opengwas2category_ontology"

   gwas_id = Column('gwas_id', String(63), primary_key=True)
   gwas_trait = Column('gwas_trait', String(255), nullable=False)
   gwas_ontology_term = Column('gwas_ontology_term', String(127), nullable=False)
   gwas_ontology_id = Column('gwas_ontology_id', String(15), nullable=False)
   gwas_ontology_iri = Column('gwas_ontology_iri', String(127), nullable=False)


<<<<<<< HEAD
class gwas_annot(Base):
   """src/insrt_gwas_annot.py"""
   __tablename__ = "gwas_annot"

   gwas_id = Column('gwas_id', String(63), primary_key=True)
   gwas_trait = Column('gwas_trait', String(255), nullable=False)
   gwas_ontology_term = Column('gwas_ontology_term', String(127), nullable=False)
   gwas_ontology_id = Column('gwas_ontology_id', String(15), nullable=False)
   gwas_category = Column('gwas_category', String(127), nullable=False)


class gwascatalog(Base):
   """src/insrt_pos19.py"""
   __tablename__ = "gwascatalog"
   __table_args__ = (UniqueConstraint('pmid', 'study', 'trait', 'mapped_trait', 'accession', name='_gwascatalog_uc'),)

   id = Column('id', Integer, primary_key=True)
   pmid = Column('pmid', Integer, nullable=False)
   study = Column('study', String(511), nullable=False)
   trait = Column('trait', String(511), nullable=False)
   mapped_trait = Column('mapped_trait', String(511), nullable=True)
   mapped_trait_uri = Column('mapped_trait_uri', String(1024), nullable=True)
   accession = Column('accession', String(15), nullable=False)
=======
# class gwas_annot(Base):
#    """src/insrt_gwas_annot.py"""
#    __tablename__ = "gwas_annot"
#
#    gwas_id = Column('gwas_id', String(63), primary_key=True)
#    gwas_trait = Column('gwas_trait', String(255), nullable=False)
#    gwas_ontology_term = Column('gwas_ontology_term', String(127), nullable=False)
#    gwas_ontology_id = Column('gwas_ontology_id', String(15), nullable=False)
#    gwas_category = Column('gwas_category', String(127), nullable=False)


# class gwascatalog(Base):
#    """src/insrt_pos19.py"""
#    __tablename__ = "gwascatalog"
#    __table_args__ = (UniqueConstraint('pmid', 'study', 'trait', 'mapped_trait', 'accession', name='_gwascatalog_uc'),)
#
#    id = Column('id', Integer, primary_key=True)
#    pmid = Column('pmid', Integer, nullable=False)
#    study = Column('study', String(511), nullable=False)
#    trait = Column('trait', String(511), nullable=False)
#    mapped_trait = Column('mapped_trait', String(511), nullable=True)
#    mapped_trait_uri = Column('mapped_trait_uri', String(1024), nullable=True)
#    accession = Column('accession', String(15), nullable=False)
>>>>>>> 7780b854f4c03d61cfedb2434aa1fd98189fe736



class ensg2symbol(Base):
   """src/insrt_gwas_annot.py"""
   __tablename__ = "ensg2symbol"
   gene_id = Column('gene_id', String(15), primary_key=True)
   gene_symbol = Column('symbol', String(63), nullable=False, unique=False)


class cytoband(Base):
   """src/insrt_cytoband.py"""
   __tablename__ = "cytoband"
   __table_args__ = (UniqueConstraint('chrom', 'start_end38', name='_cytobad2_uc'),)

   id = Column('id', String(15), primary_key=True)
   chrom = Column('chrom', SmallInteger, nullable=False)
   start_end38 = Column('start_end38', INT4RANGE, nullable=False)
   cytoband = Column('cytoband', String(7), nullable=False)


class entrezgene2ensg2symbol(Base):
   """src/insrt_entrezgene2ensg2symbol.py"""
   __tablename__ = "entrezgene2ensg2symbol"

   entrezgene = Column('entrezgene', Integer, primary_key=True)
   gene_id = Column('gene_id', String(15), nullable=False)
   gene_symbol = Column('gene_symbol', String(63), nullable=False, unique=True)


class entrezgene2pubmed_count(Base):
   """src/insrt_entrezgene2pubmed_count.py"""
   __tablename__ = "entrezgene2pubmed_count"

   entrezgene = Column('entrezgene', Integer, primary_key=True)
   pubmed_count = Column('pubmed_count', Integer, nullable=False)


class eqtl_annot(Base):
   """src/insrt_eqtl_annot.py"""
   __tablename__ = "eqtl_annot"

   eqtl_id = Column('eqtl_id', String(63), primary_key=True)
   study = Column('study', String(63), nullable=False)
   qtl_group = Column('qtl_group', String(63), nullable=False)
   tissue_ontology_id = Column('tissue_ontology_id', String(63), nullable=False)
   tissue_ontology_term = Column('tissue_ontology_term', String(63), nullable=False)
   tissue_label = Column('tissue_label', String(63), nullable=False)
   condition_label = Column('condition_label', String(63), nullable=False)
   quant_method = Column('quant_method', String(63), nullable=False)
   sample_size = Column('sample_size', String(63), nullable=False)
   ftp_path = Column('ftp_path', String(255), nullable=False)
   etissue_category_id = Column('etissue_category_id', String(15), nullable=False)
   etissue_category_term = Column('etissue_category_term', String(63), nullable=False)
   ref = Column('ref', String(63), nullable=False)


class open_gwas_info(Base):
   """src/insrt_open_gwas.py"""
   __tablename__ = "open_gwas_info"

   gwas_id = Column('gwas_id', String(127), primary_key=True)
   batch = Column('batch', String(7), nullable=False)
   note = Column('note', String(511), nullable=True)
   group_name = Column('group_name', String(127), nullable=True)
   mr = Column('mr', SmallInteger, nullable=True)
   year = Column('year', SmallInteger, nullable=True)
   author = Column('author', String(63), nullable=True)
   sex = Column('sex', String(63), nullable=True)
   pmid = Column('pmid', Integer, nullable=True)
   population = Column('population', String(127), nullable=True)
   unit = Column('unit', String(63), nullable=True)
   sample_size = Column('sample_size', Integer, nullable=True)
   nsnp = Column('nsnp', Integer, nullable=True)
   build = Column('build', String(63), nullable=True)
   trait = Column('trait', String(511), nullable=False)
   category = Column('category', String(63), nullable=True)
   subcategory = Column('subcategory', String(63), nullable=True)
   ontology = Column('ontology', String(127), nullable=True)
   consortium = Column('consortium', String(511), nullable=True)
   ncontrol = Column('ncontrol', Integer, nullable=True)
   ncase = Column('ncase', Integer, nullable=True)
   priority = Column('priority', SmallInteger, nullable=True)
   sd = Column('sd', Float, nullable=True)

class mysql_ucsc_hg38_ncbirefseq(Base):

   __tablename__ = ("genome-mysql.soe.ucsc.edu/hg38/ncbirefseq")[0:63]

   id = Column('id', Integer, primary_key=True)
   chrom = Column('chrom', SmallInteger, nullable=False)
   start38 = Column('refseq_transcript_start38', Integer, nullable=False)
   end38 = Column('refseq_transcript_end38', Integer, nullable=False)
   strand = Column('refseq_transcript_strand', String(1), nullable=False)
   refseq_transcript_id = Column('refseq_transcript_id', String(15), nullable=False)
   symbol = Column('symbol', String(63), nullable=False)
   transcript_id = Column('transcript_id', String(15), nullable=False)
   gene_id = Column('gene_id', String(15), nullable=False)
   start_end38 = Column('refseq_transcript_start_end38', INT4RANGE, nullable=False)

class af_1000genomes(Base):

   __tablename__ = "ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502"

   id = Column('id', Integer, primary_key=True)
   variant_id = Column(String(255), nullable=False)
   chrom = Column(SmallInteger, nullable=False)
   pos19 = Column(Integer, nullable=False)
   ref = Column(String(255), nullable=False)
   alt = Column(String(255), nullable=False)
   rsid = Column(Integer, nullable=False)
   eas_af = Column(Float, nullable=False)
   amr_af = Column(Float, nullable=False)
   afr_af = Column(Float, nullable=False)
   eur_af = Column(Float, nullable=False)
   sas_af = Column(Float, nullable=False)

class watanabe_posthuma2019(Base):

   __tablename__ = 'watanabe_posthuma2019'

   hgvs_id = Column(String(31), primary_key=True)
   chrom = Column(SmallInteger, nullable=False)
   pos19 = Column(Integer, nullable=False)
   ref = Column(String(255), nullable=False)
   alt = Column(String(255), nullable=False)
   rsid = Column(Integer, nullable=True)
   traits = Column(SmallInteger, nullable=False)
   domains = Column(SmallInteger, nullable=False)
   Type = Column(String(31), nullable=False)
   Domain = Column(String(31), nullable=False)
   Activities = Column(SmallInteger, nullable=False)
   Body_Structures = Column(SmallInteger, nullable=False)
   Cardiovascular = Column(SmallInteger, nullable=False)
   Cognitive = Column(SmallInteger, nullable=False)
   Connective_tissue = Column(SmallInteger, nullable=False)
   Dermatological = Column(SmallInteger, nullable=False)
   Ear_Nose_Throat = Column(SmallInteger, nullable=False)
   Endocrine = Column(SmallInteger, nullable=False)
   Environment = Column(SmallInteger, nullable=False)
   Gastrointestinal = Column(SmallInteger, nullable=False)
   Immunological = Column(SmallInteger, nullable=False)
   Infection = Column(SmallInteger, nullable=False)
   Metabolic = Column(SmallInteger, nullable=False)
   Mortality = Column(SmallInteger, nullable=False)
   Muscular = Column(SmallInteger, nullable=False)
   Neoplasms = Column(SmallInteger, nullable=False)
   Neurological = Column(SmallInteger, nullable=False)
   Nutritional = Column(SmallInteger, nullable=False)
   Ophthalmological = Column(SmallInteger, nullable=False)
   Psychiatric = Column(SmallInteger, nullable=False)
   Reproduction = Column(SmallInteger, nullable=False)
   Respiratory = Column(SmallInteger, nullable=False)
   Skeletal = Column(SmallInteger, nullable=False)
   Social_Interactions = Column(SmallInteger, nullable=False)
