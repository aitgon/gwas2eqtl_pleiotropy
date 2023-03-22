# import ssl
# ssl._create_default_https_context = ssl._create_unverified_context
import os
import pathlib

import pandas
import sys

import logging
import traceback as tb

import pandas
import suds.metrics as metrics
# from tests import *
from suds import *
from suds.client import Client
from datetime import datetime

#%%
from gwas2eqtl_pleiotropy.Logger import Logger

help_cmd_str = "todo"
try:
    count_per_rsid_gwas_egene_etissue_ods_path = sys.argv[1]
    david_email = sys.argv[2]
    davidgo_tsv_path = sys.argv[3]
    if len(sys.argv) > 4:
        print("""Two many arguments!
        {}""".format(help_cmd_str))
        sys.exit(1)
except IndexError:
    print("""Argument missing!
    {}""".format(help_cmd_str))
    sys.exit(1)

outdir_path = os.path.dirname(davidgo_tsv_path)
pathlib.Path(outdir_path).mkdir(parents=True, exist_ok=True)

#%%
# fin_df = pandas.read_csv(count_per_rsid_gwas_egene_etissue_ods_path, sep="\t", header=0)
fin_df = pandas.read_excel(count_per_rsid_gwas_egene_etissue_ods_path, engine='odf')

#%%
p_back_str = ",".join(fin_df.loc[fin_df['gwas_category_count'] == 1, "egene_lst"].str.split(';').explode().unique())

max_gwas_category_count = fin_df['gwas_category_count'].max()
for pleio_i in range(2, max_gwas_category_count+1):
    Logger.info(pleio_i)
    p_input_str = ",".join(fin_df.loc[fin_df['gwas_category_count'] == pleio_i, "egene_lst"].str.split(';').explode().unique())

    #%%

    #%%
    url = 'https://david.ncifcrf.gov/webservice/services/DAVIDWebService?wsdl'
    #
    # create a service client using the wsdl.
    client = Client(url)
    client.wsdl.services[0].setlocation(
        'https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap11Endpoint/')

    # authenticate user email
    client.service.authenticate(david_email)

    #%%
    # add a list
    idType = 'ENSEMBL_GENE_ID'
    listName = 'pleio{}'.format(pleio_i)
    listType = 0
    print(client.service.addList(p_input_str, idType, listName, listType))

    #%%
    # add a list
    idType = 'ENSEMBL_GENE_ID'
    listName = 'pleio1'
    listType = 1
    print(client.service.addList(p_back_str, idType, listName, listType))

    # %%
    # print client.service.getDefaultCategoryNames()
    categorySting = str(client.service.setCategories('GOTERM_BP_DIRECT'))
    categories = categorySting.split(',')

    #getChartReport
    thd = 0.1
    ct = 2
    chartReport = client.service.getChartReport(thd, ct)
    chartRow = len(chartReport)
    print('Total chart records:',chartRow)

    #%%
    #parse and print chartReport
    # resF = 'list1.chartReport.txt'
    davidgo_pleio_tsv_path = davidgo_tsv_path + "_{}.tsv".format(pleio_i)
    with open(davidgo_pleio_tsv_path, 'w') as fout:
        fout.write('Category\tTerm\tCount\t%\tPvalue\tGenes\tList Total\tPop Hits\tPop Total\tFold Enrichment\tBonferroni\tBenjamini\tFDR\n')
        for simpleChartRecord in chartReport:
            if not simpleChartRecord is None:
                categoryName = simpleChartRecord.categoryName
                termName = simpleChartRecord.termName
                listHits = simpleChartRecord.listHits
                percent = simpleChartRecord.percent
                ease = simpleChartRecord.ease
                Genes = simpleChartRecord.geneIds
                listTotals = simpleChartRecord.listTotals
                popHits = simpleChartRecord.popHits
                popTotals = simpleChartRecord.popTotals
                foldEnrichment = simpleChartRecord.foldEnrichment
                bonferroni = simpleChartRecord.bonferroni
                benjamini = simpleChartRecord.benjamini
                FDR = simpleChartRecord.afdr
                rowList = [categoryName,termName,str(listHits),str(percent),str(ease),Genes,str(listTotals),str(popHits),str(popTotals),str(foldEnrichment),str(bonferroni),str(benjamini),str(FDR)]
                fout.write('\t'.join(rowList) + '\n')
    print('write file:', davidgo_pleio_tsv_path, 'finished!')
