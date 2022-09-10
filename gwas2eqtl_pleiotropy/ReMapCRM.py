from gwas2eqtl_pleiotropy.Logger import Logger

import os
import pandas
import pathlib
import requests

from gwas2eqtl_pleiotropy.constants import public_data_dir


class ReMapCRM:

    # api-endpoint
    url = "https://remap.univ-amu.fr/storage/remap2022/hg38/MACS2/remap2022_crm_macs2_hg38_v1_0.bed.gz"
    remap_crm_path = os.path.join(public_data_dir, url.replace('http://', ''))

    def __init__(self):

        self._df = None

    @property
    def df(self):

        """Loads file as DF with metadata compatible header"""

        if self._df is None:
            Logger.info("Opening ReMap CRMs: " + self.remap_crm_path)

            if not os.path.isfile(self.remap_crm_path):
                self.download()
            self._df = pandas.read_excel(self.remap_crm_path, engine="odf")
            self._df = self._df.loc[self._df["population"] == "European"]

        return self._df

    @classmethod
    def download(cls):

        if not os.path.isfile(cls.remap_crm_path):
            # Logger.info("Downloading ReMap CRMs: " + cls.remap_crm_path)
            pathlib.Path(os.path.dirname(cls.remap_crm_path)).mkdir(parents=True, exist_ok=True)
            # # sending get request and saving the response as response object
            # r = requests.get(url=ReMapCRM.url)
            with requests.get(ReMapCRM.url, stream=True) as r:
                r.raise_for_status()
                with open(cls.remap_crm_path, 'wb') as f:
                    for chunk in r.iter_content(chunk_size=1024):
                        f.write(chunk)

        return cls.remap_crm_path
