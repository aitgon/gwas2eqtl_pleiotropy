import os
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
    count_per_rsid_gwas_egene_etissue_ods = sys.argv[2]
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

#########################
# default hub
#########################

hub_name = 'gwas2eqtl'
short_label = "{} colocalization".format(hub_name)
hub, genomes_file, genome, trackdb = trackhub.default_hub(
     hub_name=hub_name,
     short_label="{} colocalization".format(hub_name),
     long_label="{}: colocalization analysis of IEU OpenGWAS and EBI eQTL Catalogue".format(hub_name),
     genome="hg38",
     email="aitor.gonzalez@univ-amu.fr")

#########################
# Pleiotropic eQTLs
#########################

df = pandas.read_excel(count_per_rsid_gwas_egene_etissue_ods, engine='odf')
df = df[['chrom', 'pos38', 'rsid', 'gwas_category_count']].copy()
df['chromStart'] = df['pos38'] - 1
df['chrom'] = 'chr' + df['chrom'].astype(str)
df['rsid'] = 'rs' + df['rsid'].astype(str)

eqtl_pleio_composite_name = hub_name + "_eqtl_pleio"
eqtl_pleio_composite_short_label = "Pleiotropic eQTLs"
eqtl_pleio_composite = trackhub.CompositeTrack(
     name=eqtl_pleio_composite_name,
     short_label=eqtl_pleio_composite_short_label,
     tracktype='bigBed',
     visibility='pack')

trackdb.add_tracks(eqtl_pleio_composite)

level_colors = {1: '0,0,255', 2: '0,128,128', 3: '128,128,0', 4: '255,0,0'}
for gwas_category_count in sorted(df['gwas_category_count'].unique()):
    eqtl_pleio_path = os.path.join(os.path.dirname(count_per_rsid_gwas_egene_etissue_ods), "bigbed", "eqtl_pleio{}.bb".format(gwas_category_count))
    track_pleio = trackhub.Track(
        name="eqtl_pleitropy_{}".format(gwas_category_count),  # track names can't have any spaces or special chars.
        source=eqtl_pleio_path,  # filename to build this track from
        visibility='pack',  # shows the full signal
        color=level_colors[gwas_category_count],  # black
        tracktype='bigBed 5',  # required when making a track
    )
    eqtl_pleio_composite.add_tracks(track_pleio)

#########################
# beta equal
#########################

composite_name = hub_name + "_beta_equal"
composite_short_label = "gwas2eqtl: beta equal"
composite = trackhub.CompositeTrack(
     name=composite_name,
     short_label=composite_short_label,
     tracktype='bigInteract',
     visibility='pack')

trackdb.add_tracks(composite)

for eqtl_id in eqtl_id_lst:
    bigbed_beta_equal = os.path.join(bigbed_dir_path, "beta_equal", "{}.inter.bb".format(eqtl_id))
    if os.path.getsize(bigbed_beta_equal) > 2278:

        tissue_ontology_term = eqtl_df.loc[eqtl_id, 'tissue_ontology_term']
        track_name = trackhub.helpers.sanitize(tissue_ontology_term) + "_" + trackhub.helpers.sanitize(eqtl_id) + "_beta_equal"

        # We're keeping this relatively simple, but arguments can be
        # programmatically determined (color tracks based on sample; change scale
        # based on criteria, etc).
        short_label = "{}_{}_beta_equal".format(tissue_ontology_term, eqtl_id)
        long_label = "{} {}, beta equal".format(tissue_ontology_term, eqtl_id)
        track = trackhub.Track(
            name=track_name,          # track names can't have any spaces or special chars.
            short_label=short_label,
            long_label="{} {} beta_equal".format(tissue_ontology_term, eqtl_id),
            source=bigbed_beta_equal,      # filename to build this track from
            visibility='pack',  # shows the full signal
            color='0,0,0',    # black
            # autoScale='on',     # allow the track to autoscale
            tracktype='bigInteract',  # required when making a track
        )

        # Each track is added to the trackdb
        composite.add_tracks(track)

#########################
# beta unequal
#########################

composite_name = hub_name + "_beta_unequal"
composite_short_label = "gwas2eqtl: beta unequal"
composite = trackhub.CompositeTrack(
     name=composite_name,
     short_label=composite_short_label,
     tracktype='bigInteract',
     visibility='pack')

trackdb.add_tracks(composite)

for eqtl_id in eqtl_id_lst:
    bigbed_beta_unequal = os.path.join(bigbed_dir_path, "beta_unequal", "{}.inter.bb".format(eqtl_id))
    if os.path.getsize(bigbed_beta_unequal) > 2278:

        tissue_ontology_term = eqtl_df.loc[eqtl_id, 'tissue_ontology_term']
        track_name = trackhub.helpers.sanitize(tissue_ontology_term) + "_" + trackhub.helpers.sanitize(eqtl_id) + "_beta_unequal"

        # We're keeping this relatively simple, but arguments can be
        # programmatically determined (color tracks based on sample; change scale
        # based on criteria, etc).
        short_label = "{}_{}_beta_unequal".format(tissue_ontology_term, eqtl_id)
        long_label = "{} {}, beta unequal".format(tissue_ontology_term, eqtl_id)
        track = trackhub.Track(
            name=track_name,          # track names can't have any spaces or special chars.
            short_label=short_label,
            long_label="{} {} beta_unequal".format(tissue_ontology_term, eqtl_id),
            source=bigbed_beta_unequal,      # filename to build this track from
            visibility='pack',  # shows the full signal
            color='0,0,0',    # black
            # autoScale='on',     # allow the track to autoscale
            tracktype='bigInteract',  # required when making a track
        )

        # Each track is added to the trackdb
        composite.add_tracks(track)

#########################
# upload
#########################

# In this example we "upload" the hub locally. Files are created in the
# "example_hub" directory, along with symlinks to the tracks' data files.
# This directory can then be pushed to GitHub or rsynced to a server.

remote_dir = hub_txt_path
pathlib.Path(remote_dir).mkdir(parents=True, exist_ok=True)
trackhub.upload.upload_hub(hub=hub, host='localhost', remote_dir=remote_dir)
