from sqlalchemy import Column, Integer, String, SmallInteger, Float, UniqueConstraint
from sqlalchemy.orm import declarative_base
from sqlalchemy.dialects.postgresql import INT4RANGE

# declarative base class
Base = declarative_base()


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


class tophits(Base):
   """scripts/tophits2db2.py"""
   __tablename__ = "tophits"
   __table_args__ = (UniqueConstraint('chrom', 'pos', 'gwas_id', name='_chrom_pos_uc'),)

   id = Column('id', String(63), primary_key=True)
   chrom = Column('chrom', SmallInteger, nullable=False)
   pos = Column('pos', Integer, nullable=False)
   rsid= Column('rsid', Integer, nullable=False)
   nea = Column('nea', String(255), nullable=False)
   ea = Column('ea', String(255), nullable=False)
   pval = Column('pval', String(32), nullable=False)  # comma sep list of pvals
   beta = Column('beta', String(32), nullable=False)  # comma sep list of betas
   se = Column('se', Float, nullable=False)
   eaf = Column('eaf', Float, nullable=True)
   n = Column('n', Integer, nullable=False)
   gwas_id = Column('gwas_id', String(63), primary_key=True)
   pos19 = Column('pos19', Integer, nullable=False)
