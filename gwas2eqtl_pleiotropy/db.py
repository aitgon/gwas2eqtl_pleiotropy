from sqlalchemy import Column, Integer, String, SmallInteger, Float, UniqueConstraint
from sqlalchemy.orm import declarative_base
from sqlalchemy.dialects.postgresql import INT4RANGE

# declarative base class
Base = declarative_base()


class pos19(Base):
   """scripts/insrt_pos19.py"""
   __tablename__ = "pos19"
   __table_args__ = (UniqueConstraint('chrom', 'pos', name='_pos19_uc'),)

   id = Column('id', Integer, primary_key=True)
   chrom = Column('chrom', SmallInteger, nullable=False)
   pos = Column('pos', Integer, nullable=False)
   pos19 = Column('pos19', Integer, nullable=True)


class ensg2symbol(Base):
   """scripts/insrt_gwas_annot.py"""
   __tablename__ = "ensg2symbol"
   gene_id = Column('gene_id', String(15), primary_key=True)
   gene_symbol = Column('symbol', String(63), nullable=False, unique=False)


class cytoband(Base):
   """scripts/insrt_cytoband.py"""
   __tablename__ = "cytoband"
   __table_args__ = (UniqueConstraint('chrom', 'start_end38', name='_cytobad2_uc'),)

   id = Column('id', String(15), primary_key=True)
   chrom = Column('chrom', SmallInteger, nullable=False)
   start_end38 = Column('start_end38', INT4RANGE, nullable=False)
   cytoband = Column('cytoband', String(7), nullable=False)


class gwas_annot(Base):
   """scripts/insrt_gwas_annot.py"""
   __tablename__ = "gwas_annot"

   gwas_id = Column('gwas_id', String(63), primary_key=True)
   gwas_trait = Column('gwas_trait', String(255), nullable=False)
   gwas_class = Column('gwas_class', String(127), nullable=False)


class eqtl_annot(Base):
   """scripts/insrt_eqtl_annot.py"""
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
   etissue_class = Column('etissue_class', String(63), nullable=False)
   ref = Column('ref', String(63), nullable=False)


class open_gwas_info(Base):
   """scripts/insrt_open_gwas_info.py"""
   __tablename__ = "open_gwas_info"

   gwas_id = Column('gwas_id', String(127), primary_key=True)
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
   ontology = Column('ontology', String(63), nullable=True)
   consortium = Column('consortium', String(511), nullable=True)
   ncontrol = Column('ncontrol', Integer, nullable=True)
   ncase = Column('ncase', Integer, nullable=True)
   priority = Column('priority', SmallInteger, nullable=True)
   sd = Column('sd', Float, nullable=True)
