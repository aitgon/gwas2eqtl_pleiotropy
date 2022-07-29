import os
import pathlib
import pandas

from mysql.connector import connect
from eqtl2gwas_pleiotropy.Logger import Logger
from eqtl2gwas_pleiotropy.PathManager import PathManager


class UCSC:

    def __init__(self, database='hg38'):

        self.host = 'genome-euro-mysql.soe.ucsc.edu'
        self.user = 'genome'
        self.database = database

    def get_ref_gene_table(self, force=False):

        """Get full refgene tss table

        Force True will download data"""

        wdir = PathManager.get_outdir_path()
        ucsc_refgene_tsv_path = os.path.join(wdir, "refseq", "ucsc_refseq_genes.bed")

        if not os.path.isfile(ucsc_refgene_tsv_path) or force:
            Logger.debug("Download gene informations from Refseq")

            pathlib.Path(os.path.dirname(ucsc_refgene_tsv_path)).mkdir(exist_ok=True, parents=True)

            output_list = list()

            sql = "select chrom,txStart,txEnd,name,strand,name2 from ncbiRefSeq where (name like 'NM_%')"

            connection = connect(host=self.host, user=self.user, database=self.database)
            cursor = connection.cursor()
            Logger.info("SQL select UCSC...")
            cursor.execute(sql)
            result = cursor.fetchall()
            for row in result:
                output_list.append(row)
            cursor.close()
            connection.close()

            ucsc_refgene_df = pandas.DataFrame.from_records(output_list, columns=['chrom', 'txStart', 'txEnd', 'name', 'strand', 'name2'])
            allowed_chroms = ['chr' + str(i) for i in [*range(1, 23)] + ['Y', 'X']]
            ucsc_refgene_df = ucsc_refgene_df.loc[ucsc_refgene_df['chrom'].isin(allowed_chroms)]
            ucsc_refgene_df.to_csv(ucsc_refgene_tsv_path, index=False, header=True, sep="\t")

        ucsc_refgene_df = pandas.read_csv(ucsc_refgene_tsv_path, sep='\t', header=0, index_col='name')

        return ucsc_refgene_df



    def get_ensg2enst2ensp(self, force=False):

        """Get full refgene tss table

        Force True will download data"""

        # wdir = PathManager.get_outdir_path()
        ucsc_ensg2enst2ensp_tsv_path = os.path.join(PathManager.get_outdir_path(), "ensg2enst2ensp.tsv")

        if not os.path.isfile(ucsc_ensg2enst2ensp_tsv_path) or force:
            Logger.debug("Download gene informations from Refseq")

            pathlib.Path(os.path.dirname(ucsc_ensg2enst2ensp_tsv_path)).mkdir(exist_ok=True, parents=True)

            output_list = list()

            # sql = "select chrom,txStart,txEnd,name,strand,name2 from ncbiRefSeq where (name like 'NM_%' or name like 'XM_%')"
            sql = "select * from ensGtp"

            connection = connect(host=self.host, user=self.user, database=self.database)
            cursor = connection.cursor()
            Logger.info("SQL select UCSC...")
            cursor.execute(sql)
            result = cursor.fetchall()
            for row in result:
                output_list.append(row)
            cursor.close()
            connection.close()

            df = pandas.DataFrame.from_records(output_list, columns=['gene', 'transcript', 'protein'])
            df.to_csv(ucsc_ensg2enst2ensp_tsv_path, index=False, header=True, sep="\t")

        df = pandas.read_csv(ucsc_ensg2enst2ensp_tsv_path, sep='\t', header=0)

        return df
