import sys
import os
import argparse
import pandas as pd
from multiprocessing import Pool, Manager

import metadata
import annotation
from jakomics import utilities, kegg, colors

version = "v0.3.0"

print(f'{colors.bcolors.GREEN}Genome annotATOR (GATOR) {version} (Under active development!!){colors.bcolors.END}')

# OPTIONS #####################################################################

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


def annotate(genome):
    global annotated_genomes

    # run genome against databases. each method type will need its own logic
    genome.raw_results = {}
    for id, db in metadata.db_info.iterrows():
        genome.raw_results[db['DB_NAME']] = {}

        # kofam method
        if db['METHOD'] == 'kofam':
            hits = kegg.run_kofam(genome.file_path, db['hal_path'])
            genome.raw_results[db['DB_NAME']] = kegg.parse_kofam_hits(hits)

    details = pd.DataFrame(columns=['GENOME', 'GENE', 'PRODUCT', 'TYPE', 'ID',
                                    'LOCUS_TAG', 'SCORE', 'EVAL', 'NOTE', 'COMPLEX', "REACTION"])

    for gene_index, gene in metadata.gene_info.iterrows():
        for db_index, db in metadata.db_info.iterrows():
            term_list = [x.strip() for x in gene[db['DB_NAME']].split(',')]
            for term in term_list:
                if term in genome.raw_results[db['DB_NAME']]:
                    for hit in genome.raw_results[db['DB_NAME']][term]:
                        r = hit.result()

                        details = details.append(
                            pd.Series(data={'GENOME': genome.short_name,
                                            'GENE': gene['GENE_NAME'],
                                            'PRODUCT': gene['GENE_PRODUCT'],
                                            'TYPE': db['DB_NAME'],
                                            'ID': r['annotation'],
                                            'LOCUS_TAG': r['gene'],
                                            'SCORE': r['score'],
                                            'EVAL': r['evalue'],
                                            'NOTE': gene['GENE_NOTE'],
                                            'COMPLEX': gene['COMPLEX'],
                                            'REACTION': gene['REACTION']
                                            }
                                      ),
                            ignore_index=True)
    genome.annotations = details
    annotated_genomes.append(genome)


# MAIN ########################################################################
pool = Pool()
pool.map(annotate, unannotated_genomes)
pool.close()

# cleanup
metadata.remove_temp_files()

pathway_results = pd.DataFrame(columns=[
    'GENOME',
    'PATHWAY',
    'PRESENT',
    'PATHWAY_STEPS',
    'STEPS_PRESENT',
    'REACTION',
    'GENES',
    'PATHWAY_DEFINITION'])

for genome in annotated_genomes:

    genes = list(set(genome.annotations['GENE']))

    for p in metadata.pathways:
        results = p.score_pathway(genes, genome.short_name)
        pathway_results = pathway_results.append(results, ignore_index=True)

pathway_results.to_csv("pathway_results.txt", sep="\t", index=False)
