from sqlalchemy import Column, Integer, String, SmallInteger, Float, UniqueConstraint
from sqlalchemy.orm import declarative_base

# declarative base class
Base = declarative_base()


class coloc(Base):
   """scripts/ins_coloc.py"""
   __tablename__ = "coloc"
   __table_args__ = (UniqueConstraint('chrom', 'pos', 'egene_id', 'gwas_id', 'eqtl_id', 'coloc_lead_pos', name='_coloc2_uc'),)

   id = Column('id', String(255), primary_key=True)
   chrom = Column('chrom', SmallInteger, nullable=False)
   pos = Column('pos', Integer, nullable=False)
   rsid= Column('rsid', Integer, nullable=False)
   ref = Column('ref', String(255), nullable=False)
   alt = Column('alt', String(255), nullable=False)
   egene_id = Column('egene_id', String(15), nullable=False)
   gwas_beta = Column('gwas_beta', Float, nullable=False)
   gwas_pval = Column('gwas_pval', Float, nullable=False)
   gwas_id = Column('gwas_id', String(255), nullable=False)
   eqtl_id = Column('eqtl_id', String(255), nullable=False)
   eqtl_beta = Column('eqtl_beta', Float, nullable=False)
   eqtl_pval = Column('eqtl_pval', Float, nullable=False)
   pp_h4_abf = Column('PP.H4.abf', Float, nullable=False)
   snp_pp_h4 = Column('SNP.PP.H4', Float, nullable=False)
   pp_h3_abf = Column('PP.H3.abf', Float, nullable=False)
   pp_h2_abf = Column('PP.H2.abf', Float, nullable=False)
   pp_h1_abf = Column('PP.H1.abf', Float, nullable=False)
   pp_h0_abf = Column('PP.H0.abf', Float, nullable=False)
   nsnps = Column('nsnps', SmallInteger, nullable=False)
   coloc_lead_pos = Column('coloc_lead_pos', Integer, nullable=False)
   rscoloc_lead_rsidid= Column('coloc_lead_rsid', Integer, nullable=False)
   coloc_region= Column('coloc_region', String(255), nullable=False)


class ensg2symbol(Base):
   """scripts/annotate_db2.py"""
   __tablename__ = "ensg2symbol"
   gene_id = Column('gene_id', String(15), primary_key=True)
   gene_symbol = Column('symbol', String(63), nullable=False, unique=True)


class cytoband(Base):
   """scripts/annotate_db2.py"""
   __tablename__ = "cytoband"
   __table_args__ = (
   UniqueConstraint('chrom', 'start', name='_cytobad_uc'),)

   id = Column('id', String(15), primary_key=True)
   chrom = Column('chrom', SmallInteger, nullable=False)
   start = Column('start', Integer, nullable=False)
   end = Column('end', Integer, nullable=False)
   cytoband = Column('cytoband', String(7), nullable=False)


class gwas_annot(Base):
   """scripts/annotate_db2.py"""
   __tablename__ = "gwas_annot"

   gwas_id = Column('gwas_id', String(63), primary_key=True)
   gwas_trait = Column('gwas_trait', String(255), nullable=False, unique=True)
   gwas_class = Column('gwas_class', String(127), nullable=False)

