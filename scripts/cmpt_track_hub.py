import glob, os
import pathlib
import sys

import pandas
import trackhub

"""
https://genome.ucsc.edu/goldenPath/help/trackDb/trackDbHub.html
https://s3.mpi-bn.mpg.de/data-tobias-ucsc/hg38/trackDb.txt
"""

#%%
help_cmd_str = "todo"
try:
    hub_name = sys.argv[1]
    bigbed_dir_path = sys.argv[2]
    hub_txt_path = sys.argv[3]
    if len(sys.argv) > 4:
        print("""Two many arguments!
        {}""".format(help_cmd_str))
        sys.exit(1)
except IndexError:
    print("""Argument missing!
    {}""".format(help_cmd_str))
    sys.exit(1)

eqtl_tsv_path = 'https://raw.githubusercontent.com/eQTL-Catalogue/eQTL-Catalogue-resources/master/tabix/tabix_ftp_paths.tsv'
eqtl_df = pandas.read_csv(eqtl_tsv_path, sep="\t", usecols=[0, 3, 4, 6, 8])
eqtl_df.sort_values(by=['tissue_ontology_term', 'tissue_label', 'study'], inplace=True)
eqtl_df = eqtl_df.loc[eqtl_df['ftp_path'].str.contains('/ge/|/microarray/',regex=True,na=False), ]

# keep urls containing "ge" or "microarray"
eqtl_df['index'] = (eqtl_df['ftp_path'].str.replace('.all.tsv.gz', '', regex=True)).str.split('/', expand=True)[10].tolist()
eqtl_df.set_index('index', drop=True, verify_integrity=True, inplace=True)
eqtl_id_lst = (eqtl_df.index).tolist()

short_label = "{} colocalization".format(hub_name)
hub, genomes_file, genome, trackdb = trackhub.default_hub(
     hub_name=hub_name,
     short_label="{} colocalization".format(hub_name),
     long_label="{}: colocalization analysis of IEU OpenGWAS and EBI eQTL Catalogue".format(hub_name),
     genome="hg38",
     email="aitor.gonzalez@univ-amu.fr")

composite = trackhub.CompositeTrack(
     name=hub_name,
     short_label=short_label,
     tracktype='bigInteract',
     visibility='pack')

trackdb.add_tracks(composite)

for eqtl_id in eqtl_id_lst:
    # eqtl2_id = trackhub.helpers.sanitize(os.path.basename(eqtl_id))
    bigbed = os.path.join(bigbed_dir_path, "{}.inter.bb".format(eqtl_id))
    if os.path.getsize(bigbed) > 2278:

        tissue_ontology_term = eqtl_df.loc[eqtl_id, 'tissue_ontology_term']
        track_name = trackhub.helpers.sanitize(tissue_ontology_term) + "_" + trackhub.helpers.sanitize(eqtl_id)

        # We're keeping this relatively simple, but arguments can be
        # programmatically determined (color tracks based on sample; change scale
        # based on criteria, etc).
        short_label = "{}_{}".format(tissue_ontology_term, eqtl_id)
        long_label = "{} {}".format(tissue_ontology_term, eqtl_id)
        track = trackhub.Track(
            name=track_name,          # track names can't have any spaces or special chars.
            short_label=short_label,
            long_label="{} {}".format(tissue_ontology_term, eqtl_id),
            source=bigbed,      # filename to build this track from
            visibility='full',  # shows the full signal
            color='128,0,5',    # brick red
            # autoScale='on',     # allow the track to autoscale
            tracktype='bigInteract',  # required when making a track
        )

        # Each track is added to the trackdb

        composite.add_tracks(track)

# In this example we "upload" the hub locally. Files are created in the
# "example_hub" directory, along with symlinks to the tracks' data files.
# This directory can then be pushed to GitHub or rsynced to a server.

remote_dir = hub_txt_path
pathlib.Path(remote_dir).mkdir(parents=True, exist_ok=True)
trackhub.upload.upload_hub(hub=hub, host='localhost', remote_dir=remote_dir)