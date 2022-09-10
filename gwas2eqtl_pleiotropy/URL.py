import gzip
import os
import pathlib
import sys
import urllib.request
import requests

from gwas2eqtl_pleiotropy.Logger import Logger
from gwas2eqtl_pleiotropy.constants import public_data_dir


class URL:

    def __init__(self, url, download_dir=None):

        self.url = url
        if download_dir is None:
            download_dir = public_data_dir

        url_path = url.split('//')[1]
        self.local_file_path = os.path.join(download_dir, url_path)
        self.local_dir_path = os.path.dirname(self.local_file_path)
        pathlib.Path(self.local_dir_path).mkdir(parents=True, exist_ok=True)

    def download(self):

        if not os.path.isfile(self.local_file_path):

            Logger.info(
                msg="Downloading: {}".format(self.url))

            try:  # method 1 to download
                urllib.request.urlretrieve(self.url, self.local_file_path)
            except urllib.error.HTTPError:  # method 2 to download
                r = requests.get(self.url)
                with open(self.local_file_path, 'wb') as outfile:
                    outfile.write(r.content)

        if not os.path.isfile(self.local_file_path):
            Logger.error(msg='File could not be downloaded. Script will exit')
            sys.exit(1)

        return self.local_file_path

    def gunzip(self, gunzip_path=None):

        if gunzip_path is None:
            gunzip_path = (self.local_file_path).split('.gz')[-2]  # gunzipped:

        pathlib.Path(os.path.dirname(gunzip_path)).mkdir(parents=True, exist_ok=True)

        if not os.path.isfile(gunzip_path):
            Logger.debug(msg="Gunzipping: {}".format(gunzip_path))
            with gzip.open(self.local_file_path, 'rb') as f:
                with open(gunzip_path, 'wb') as fout:
                    fout.write(f.read())

        return gunzip_path
