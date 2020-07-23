import sys
import os
import argparse
import pandas as pd
from multiprocessing import Pool, Manager

import metadata
import annotation
from jakomics import utilities, kegg, colors

version = "v0.2.0"

print(f'{colors.bcolors.GREEN}Genome annotATOR (GATOR) {version} (Under active development!!){colors.bcolors.END}')

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

manager = Manager()
annotated_genomes = manager.list()

# prep metadata and databases
gator_path = (os.path.dirname(os.path.abspath(sys.argv[0])))
metadata = metadata.Metadata(os.path.join(gator_path, "gator.xlsx"))
metadata.create_hal_files()

# get genomes
unannotated_genomes = utilities.get_files(args.files, args.in_dir, ["faa"])

# annotation genomes


def main(genome):
    global annotated_genomes

    # run genome against databases. each method type will need its own logic
    genome.raw_results = {}
    for id, db in metadata.db_info.iterrows():
        genome.raw_results[db['DB_NAME']] = {}

        # kofam method
        if db['METHOD'] == 'kofam':
            hits = kegg.run_kofam(genome.file_path, db['hal_path'])
            genome.raw_results[db['DB_NAME']] = kegg.parse_kofam_hits(hits)

    # add raw results to gene annotations - should be generic to each database
    genome.annotations = {}
    for id, reference_gene in metadata.gene_info.iterrows():
        gene = reference_gene['GENE_NAME']
        genome.annotations[gene] = annotation.Annotation(reference_gene)

    # #
    for gene_index, gene in metadata.gene_info.iterrows():
        for db_index, db in metadata.db_info.iterrows():
            term_list = [x.strip() for x in gene[db['DB_NAME']].split(',')]
            for term in term_list:
                if term in genome.raw_results[db['DB_NAME']]:
                    for hit in genome.raw_results[db['DB_NAME']][term]:
                        r = hit.result()
                        print(genome.short_name,
                              gene['GENE_NAME'],
                              r['gene'],
                              r['annotation'],
                              gene['GENE_PRODUCT'],
                              gene['GENE_NOTE'],
                              gene['COMPLEX'],
                              gene['REACTION'],
                              db['DB_NAME'],
                              r['score'],
                              r['evalue'],
                              sep="\t")

    annotated_genomes.append(genome)


#
pool = Pool()
pool.map(main, unannotated_genomes)
pool.close()

# cleanup
metadata.remove_temp_files()
