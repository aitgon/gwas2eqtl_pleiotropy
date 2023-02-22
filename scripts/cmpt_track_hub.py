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
    bigbed_dir_path = sys.argv[1]
    hub_txt_path = sys.argv[2]
    if len(sys.argv) > 3:
        print("""Two many arguments!
        {}""".format(help_cmd_str))
        sys.exit(1)
except IndexError:
    print("""Argument missing!
    {}""".format(help_cmd_str))
    sys.exit(1)

# First we initialize the components of a track hub

# hub, genomes_file, genome, trackdb = trackhub.default_hub(
#     hub_name="gwas2eqtl",
#     short_label='GWAS/eQTL colocalization',
#     long_label='Colocalization analysis of IEU OpenGWAS and EBI eQTL Catalogue',
#     genome="hg38",
#     email="aitor.gonzalez@univ-amu.fr")

eqtl_tsv_path = 'https://raw.githubusercontent.com/eQTL-Catalogue/eQTL-Catalogue-resources/master/tabix/tabix_ftp_paths.tsv'
eqtl_df = pandas.read_csv(eqtl_tsv_path, sep="\t", usecols=[0, 3, 6, 8])
eqtl_df = eqtl_df.loc[eqtl_df['ftp_path'].str.contains('/ge/|/microarray/',regex=True,na=False), ]
# if "eqtl_id" in config:
#     eqtl_id = config["eqtl_id"] # path to the eqtl< ods list
#     eqtl_df = eqtl_df.loc[eqtl_df['ftp_path'].str.contains(eqtl_id)]
# keep urls containing "ge" or "microarray"
eqtl_df['index'] = (eqtl_df['ftp_path'].str.replace('.all.tsv.gz', '', regex=True)).str.split('/', expand=True)[10].tolist()
eqtl_df.set_index('index', drop=True, verify_integrity=True, inplace=True)
eqtl_id_lst = (eqtl_df.index).tolist()

hub, genomes_file, genome, trackdb = trackhub.default_hub(
     hub_name="gwas2eqtl",
     short_label="gwas2eqtl colocalization",
     long_label="gwas2eqtl: colocalization analysis of IEU OpenGWAS and EBI eQTL Catalogue",
     genome="hg38",
     email="aitor.gonzalez@univ-amu.fr")

composite = trackhub.CompositeTrack(
     name='gwas2eqtl',
     short_label='gwas2eqtl colocalization',
     tracktype='bigInteract',
     visibility='full')

trackdb.add_tracks(composite)

# signal_view = trackhub.ViewTrack(
#      name='group',
#      view='group_view',
#      visibility='full',
#      tracktype='bigInteract',
#      short_label='Signal')
# composite.add_tracks(signal_view)

# Next we add tracks for some bigWigs. These can be anywhere on the
# filesystem; symlinks will be made to them. Here we use some example data
# included with the trackhub package; in practice you'd point to your own
# data.

# track footprint_scores
# compositeTrack on
# shortLabel footprint scores
# longLabel TOBIAS footprint scores
# visibility full
# priority 3
# html docs/footprints
# type bigWig

# track = trackhub.Track(
#     name='gwas2eqtl',  # track names can't have any spaces or special chars.
#     compositeTrack=True,
#     visibility='full',  # shows the full signal
#     priority=True,  # brick red
#     tracktype='bigInteract',  # required when making a track
# )
#
# trackdb.add_tracks(track)


for eqtl_id in eqtl_id_lst:
    bigbed = os.path.join(bigbed_dir_path, "{}.inter.bb".format(eqtl_id))
    track_name = "gwas2eqtl_{}".format(eqtl_id)

    # for bigbed in glob.glob(os.path.join(bigbed_dir_path, '*.bb')):

    # track names can't have any spaces or special characters. Since we'll
    # be using filenames as names, and filenames have non-alphanumeric
    # characters, we use the sanitize() function to remove them.

    track_name = trackhub.helpers.sanitize(os.path.basename(track_name))

    # We're keeping this relatively simple, but arguments can be
    # programmatically determined (color tracks based on sample; change scale
    # based on criteria, etc).
    # import pdb; pdb.set_trace()
    track = trackhub.Track(
        name=track_name,          # track names can't have any spaces or special chars.
        short_label=eqtl_id,
        long_label="gwas2eqtl {}: Colocalization with 417 IEU OpenGWAS".format(eqtl_id, eqtl_id),
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
