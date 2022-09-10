import shlex
import subprocess
import sys

from sqlalchemy import create_engine, text, table
from sqlalchemy.testing.schema import Table

from gwas2eqtl_pleiotropy.db import meta, coloc

help_cmd_str = "todo"
try:
    sqlalchemy_url = sys.argv[1]
    coloc_tsv_gz = sys.argv[2]
    if len(sys.argv) > 3:
        print("""Two many arguments!
        {}""".format(help_cmd_str))
        sys.exit(1)
except IndexError:
    print("""Argument missing!
    {}""".format(help_cmd_str))
    sys.exit(1)


engine = create_engine(sqlalchemy_url, echo = True)
meta.create_all(engine)
cols = ["chrom", "pos", "rsid", "ref", "alt", "egene", "gwas_beta", "gwas_pval", "gwas_id", "eqtl_beta", "eqtl_pval", "eqtl_id", "PP.H4.abf", "SNP.PP.H4", "nsnps", "PP.H3.abf", "PP.H2.abf", "PP.H1.abf", "PP.H0.abf", "coloc_lead_pos", "coloc_lead_rsid", "coloc_region"]

colocimport = Table('colocimport', meta, autoload_with=engine)
ins = coloc.insert().from_select(cols, colocimport.select())
with engine.connect() as con:
    con.execute(ins)
import pdb; pdb.set_trace()
