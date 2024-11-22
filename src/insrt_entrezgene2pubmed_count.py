import os
import pathlib
import pooch
import pandas
# import requests
import sqlalchemy
import sys

# from gwas2eqtl_pleiotropy.constants import public_data_dir
from gwas2eqtl_pleiotropy_db import Base, entrezgene2pubmed_count
from sqlalchemy import create_engine


#%%
help_cmd_str = "todo"
try:
    url = sys.argv[1]
    if len(sys.argv) > 2:
        print("""Two many arguments!
        {}""".format(help_cmd_str))
        sys.exit(1)
except IndexError:
    print("""Argument missing!
    {}""".format(help_cmd_str))
    sys.exit(1)

#%% Create all tables
engine = create_engine(url)
if sqlalchemy.inspect(engine).has_table("entrezgene2pubmed_count"):
    entrezgene2pubmed_count.__table__.drop(engine)
# Base.metadata.create_all(engine)
Base.metadata.tables["entrezgene2pubmed_count"].create(bind=engine)

# #%% Download
# url = "https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2pubmed.gz"
# tsv_path = os.path.join(public_data_dir, url.split('https://')[1])
# if not os.path.isfile(tsv_path):
#     pathlib.Path(os.path.dirname(tsv_path)).mkdir(exist_ok=True, parents=True)
#     resp = requests.get(url)  # making requests to server
#     with open(tsv_path, "wb") as f:  # opening a file handler to create new file
#         f.write(resp.content)  # writing content to file

#%%
# Download a file and save it locally, returning the path to it.
# Running this again will not cause a download. Pooch will check the hash
# (checksum) of the downloaded file against the given value to make sure
# it's the right file (not corrupted or outdated).
"""
tsv_path = pooch.retrieve(
    url="https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2pubmed.gz",
    known_hash="md5:39d68c1d5a7cf79e3c1785ad41bc1c69",
)
"""

poochobj = pooch.create(
    # Download dir
    path=os.path.join(os.environ['PUBLIC_DATA_DIR'], "ftp.ncbi.nlm.nih.gov/gene/DATA"),
    base_url="https://ftp.ncbi.nlm.nih.gov/gene/DATA/",
    version_dev="main",
    registry={
        "gene2pubmed.gz": "md5:2eab466e616d6b2d2d280dd5ee7480cf",
    },
)
gene2pubmed_gz_path = poochobj.fetch("gene2pubmed.gz")

#%%
df = pandas.read_csv(gene2pubmed_gz_path, sep="\t", header=0)
df = df.loc[df['#tax_id'] == 9606]
df.drop(['#tax_id'], inplace=True, axis=1)

df = df.groupby('GeneID').size().reset_index()
df.rename({'GeneID': 'entrezgene', 0: 'pubmed_count'}, inplace=True, axis=1)

df.to_sql('entrezgene2pubmed_count', con=engine, if_exists='append', index=False)
