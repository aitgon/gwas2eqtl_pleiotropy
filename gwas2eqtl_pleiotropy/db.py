from sqlalchemy import Table, Column, Integer, String, MetaData, Float, UniqueConstraint, ForeignKey

meta = MetaData()

rsid2cytoband_tbl = Table(
   'rsid2cytoband', meta,
   Column('rsid', String(50), primary_key=True),
   Column('cytoband', String(50), nullable=False),
)

ensg2symbol_tbl = Table(
   'ensg2symbol', meta,
   Column('gene_id', String(50), primary_key=True),
   Column('symbol', String(50), nullable=False),
)

gwas_annot_tbl = Table(
   'gwas_annot', meta,
   Column('gwas_id', String(50), primary_key=True),
   Column('gwas_class', String(50), nullable=False),
   Column('gwas_trait', String(50), nullable=False),
)

coloc = Table(
   'coloc', meta,
   Column('id', Integer, primary_key = True),
   Column('chrom', Integer, nullable=False),
   Column('pos', Integer, nullable=False),
   Column('rsid', String(50), ForeignKey("rsid2cytoband.rsid"), nullable=False),
   Column('ref', String(50), nullable=False),
   Column('alt', String(50), nullable=False),
   Column('egene', String(50), ForeignKey("ensg2symbol.gene_id"), nullable=False),
   Column('gwas_beta', Float, nullable=False),
   Column('gwas_pval', Float, nullable=False),
   Column('gwas_id', String(50), ForeignKey("gwas_annot.gwas_id"), nullable=False),
   Column('eqtl_id', String(50), nullable=False),
   Column('eqtl_beta', Float, nullable=False),
   Column('eqtl_pval', Float, nullable=False),
   Column('PP.H4.abf', Float, nullable=False),
   Column('SNP.PP.H4', Float, nullable=False),
   Column('PP.H3.abf', Float, nullable=False),
   Column('PP.H2.abf', Float, nullable=False),
   Column('PP.H1.abf', Float, nullable=False),
   Column('PP.H0.abf', Float, nullable=False),
   Column('nsnps', Integer, nullable=False),
   Column('coloc_lead_pos', Integer, nullable=False),
   Column('coloc_lead_rsid', String(50), nullable=False),
   Column('coloc_region', String(50), nullable=False),
)
