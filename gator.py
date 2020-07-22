import sys
import os
import argparse

import pandas as pd
import metadata
from jakomics import utilities, kegg, colors

version = "v0.1.0"

print(f'{colors.bcolors.GREEN}Gator {version}{colors.bcolors.END}')

parser = argparse.ArgumentParser(description="", formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('--in_dir',
                    help="Directory with faa genomes",
                    required=False,
                    default="")

parser.add_argument('-f', '--files',
                    help="Paths to individual faa files",
                    nargs='*',
                    required=False,
                    default=[])

args = parser.parse_args()


# prep metadata and databases
gator_path = (os.path.dirname(os.path.abspath(sys.argv[0])))
metadata = metadata.Metadata(os.path.join(gator_path, "gator.xlsx"))
metadata.create_hal_files()

# get genomes
file_list = utilities.get_files(args.files, args.in_dir, ["faa"])
for genome in file_list:
    print(f'working on {genome.short_name}')
    genome.annotations = {}

    for id, db in metadata.db_info.iterrows():
        genome.annotations[db['DB_NAME']] = {}

        if db['METHOD'] == 'kofam':
            hits = kegg.run_kofam(genome.file_path, db['hal_path'])
            genome.annotations[db['DB_NAME']] = kegg.parse_kofam_hits(hits)

metadata.remove_temp_files()

# temp results
for genome in file_list:
    print("---")
    for id, gene_info in metadata.gene_info.iterrows():
        for db in genome.annotations:
            for annotation in genome.annotations[db]:
                for hit in genome.annotations[db][annotation]:
                    if hit.result()['annotation'] in gene_info[db]:
                        print(genome.short_name, gene_info['GENE_NAME'],
                              db, annotation, hit.view(), sep="\t")
