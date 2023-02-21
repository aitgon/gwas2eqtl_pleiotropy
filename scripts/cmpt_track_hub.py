import glob, os
import pathlib
import sys

import trackhub

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

hub, genomes_file, genome, trackdb = trackhub.default_hub(
    hub_name="gwas2eqtl",
    short_label='GWAS/eQTL colocalization',
    long_label='Colocalization analysis of IEU OpenGWAS and EBI eQTL Catalogue',
    genome="hg38",
    email="aitor.gonzalez@univ-amu.fr")

# Next we add tracks for some bigWigs. These can be anywhere on the
# filesystem; symlinks will be made to them. Here we use some example data
# included with the trackhub package; in practice you'd point to your own
# data.

for bigbed in glob.glob(os.path.join(bigbed_dir_path, '*.bb')):

    # track names can't have any spaces or special characters. Since we'll
    # be using filenames as names, and filenames have non-alphanumeric
    # characters, we use the sanitize() function to remove them.

    name = trackhub.helpers.sanitize(os.path.basename(bigbed))

    # We're keeping this relatively simple, but arguments can be
    # programmatically determined (color tracks based on sample; change scale
    # based on criteria, etc).

    track = trackhub.Track(
        name=name,          # track names can't have any spaces or special chars.
        source=bigbed,      # filename to build this track from
        visibility='full',  # shows the full signal
        color='128,0,5',    # brick red
        # autoScale='on',     # allow the track to autoscale
        tracktype='bigInteract',  # required when making a track
    )

    # Each track is added to the trackdb

    trackdb.add_tracks(track)

# In this example we "upload" the hub locally. Files are created in the
# "example_hub" directory, along with symlinks to the tracks' data files.
# This directory can then be pushed to GitHub or rsynced to a server.

remote_dir = hub_txt_path
pathlib.Path(remote_dir).mkdir(parents=True, exist_ok=True)
trackhub.upload.upload_hub(hub=hub, host='localhost', remote_dir=remote_dir)
