import sys
import os
import argparse
import pandas as pd
from multiprocessing import Pool
from tqdm import tqdm
import ast
import datetime

import metadata as m
from jakomics import utilities, kegg, colors, blast, hmm
from jakomics.genome import GENOME
from jakomics.patric import Patric_Gene
from jakomics.file import FILE

version = "v1.3.1"

print(f'{colors.bcolors.GREEN}Genome annotATOR (GATOR) {version}{colors.bcolors.END}')

# OPTIONS #####################################################################

parser = argparse.ArgumentParser(description="",
                                 formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('--in_dir',
                    help="Directory with faa genomes",
                    required=False,
                    default="")

parser.add_argument('-db', '--gator_db',
                    help="Excel file with custom gator db",
                    required=False,
                    default="default")

parser.add_argument('-f', '--files',
                    help="Paths to individual faa files",
                    nargs='*',
                    required=False,
                    default=[])

parser.add_argument('--verify_db',
                    action='store_true',
                    help='Just check the database')

parser.add_argument('-p', '--patric',
                    action='store_true',
                    help='Genbank files are from patric. Adds patric and EC db support')

parser.add_argument('--save_raw',
                    action='store_true',
                    help='Save the raw search data to files')

parser.add_argument('-n', '--name',
                    help="Name to prepend the output files. Defaults to a timestamp",
                    required=False,
                    default=datetime.datetime.now().strftime("%Y%m%d_%H%M%S"))

args = parser.parse_args()

# Functions

def prep_metadata():
    # prep metadata and databases
    
    if args.gator_db == "default":
        gator_path = (os.path.dirname(os.path.abspath(sys.argv[0])))
        metadata = m.Metadata(os.path.join(gator_path, "gator_db.xlsx"))
    else:
        metadata = m.Metadata(args.gator_db)

    metadata.create_hal_files()
    metadata.make_blast_dbs()
    metadata.verify_metadata()

    metadata.summary()

    return(metadata)

metadata = prep_metadata()

def annotate(genome):

    if genome.suffix in ['.gb', '.gbk', '.gbff']:
        gbk = GENOME(genome)

        # write genes to genomes and gene class dictionary
        gbk.faa_path = genome.short_name + "_" + genome.id + ".faa"
        genome.patric = gbk.genbank_to_fasta(
            write_faa=gbk.faa_path, return_gene_dict=True)
        genome.file_path = gbk.faa_path
        genome.temp_files['faa'] = gbk.faa_path

    # run genome against databases. each method type will need its own logic
    genome.raw_results = {}
    for id, db in metadata.db_info.iterrows():
        genome.raw_results[db['DB_NAME']] = {}

        # raw results are dicts with term as key, and list of objects as values
        if db['METHOD'] == 'kofam':
            hits = kegg.run_kofam(genome.file_path, db['hal_path'], os.path.join(
                db['DB_PATH'], 'ko_list'), echo=False)
            genome.raw_results[db['DB_NAME']] = kegg.parse_kofam_hits(hits)

        elif db['METHOD'] == 'blastp':
            genome.raw_results[db['DB_NAME']] = blast.run_blast(type="prot",
                                                                q=genome.file_path,
                                                                db=db['DB_PATH'],
                                                                e=1e-15,
                                                                make=False,
                                                                return_query_results=False)
        elif db['METHOD'] == 'hmm':
            genome.temp_log = genome.id + '.hmm.log'
            genome.temp_output = genome.id + '.hmm.temp.txt'

            r = hmm.run_hmmsearch(genome.file_path,
                              genome.temp_log,
                              genome.temp_output,
                              db['DB_PATH'],
                              cut_tc=False,  # hot fix
                              echo=False,
                              run=True)

            if r != ['']:
                print(f"{colors.bcolors.YELLOW}ERROR: hmmsearch threw the following error while searching {genome.name} against {db['DB_PATH']}{colors.bcolors.END}")
                for line in r:
                    if len(line) > 0:
                        print(f'\t-{colors.bcolors.YELLOW}{line.rstrip()}{colors.bcolors.END}')

            genome.raw_results[db['DB_NAME']] = hmm.parse_hmm_hits(
                genome.temp_output)
            os.system('rm ' + genome.temp_log)
            os.system('rm ' + genome.temp_output)

        elif db['METHOD'] == 'EC':
            if args.patric:
                for gene in genome.patric:
                    if hasattr(genome.patric[gene], 'EC_number'):
                        patric_annotation = Patric_Gene(genome.patric[gene].id)
                        for EC in genome.patric[gene].EC_number:
                            patric_annotation.annotation = EC

                            if EC in genome.raw_results[db['DB_NAME']]:
                                genome.raw_results[db['DB_NAME']][EC].append(
                                    patric_annotation)
                            else:
                                genome.raw_results[db['DB_NAME']][EC] = [
                                    patric_annotation]

        elif db['METHOD'] == 'PATRIC':
            if args.patric:
                for gene in genome.patric:
                    if hasattr(genome.patric[gene], 'product'):
                        patric_annotation = Patric_Gene(genome.patric[gene].id)
                        for product in genome.patric[gene].product:
                            patric_annotation.annotation = product

                            if product in genome.raw_results[db['DB_NAME']]:
                                genome.raw_results[db['DB_NAME']][product].append(
                                    patric_annotation)
                            else:
                                genome.raw_results[db['DB_NAME']][product] = [
                                    patric_annotation]

    # Get Results
    details = pd.DataFrame(columns=['GENOME', 'GENE', 'PRODUCT', 'TYPE', 'ID',
                                    'LOCUS_TAG', 'SCORE', 'EVAL', 'NOTE', 'COMPLEX', "REACTION"])

    # parse details output file
    ## the result method needs to return the same data regardless of class
    for gene_index, gene in metadata.gene_info.iterrows():
        for db_index, db in metadata.db_info.iterrows():
            if pd.notnull(gene[db['DB_NAME']]):
                term_list = [x.strip() for x in gene[db['DB_NAME']].split(',')]

                if db['DB_NAME'] == "PATRIC":
                    '''
                    Patric db is formatted as a list, rather than comma-delimited
                    '''
                    term_list = [x.strip()
                                 for x in ast.literal_eval(gene[db['DB_NAME']])]

                for term in term_list:
                    if term in genome.raw_results[db['DB_NAME']]:
                        for hit in genome.raw_results[db['DB_NAME']][term]:
                            r = hit.result()

                            results = pd.Series(data={'GENOME': genome.short_name,
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
                                          )

                            details = pd.concat([details, results.to_frame().T], ignore_index=True)

    details.to_csv(f"{genome.name}_gator_detail.txt", sep="\t", index=False)


    # Score potatoes
    genes = list(set(details['GENE']))
    pathway_results = pd.DataFrame(columns=['GENOME',
                                            'PATHWAY',
                                            'PRESENT',
                                            'PATHWAY_STEPS',
                                            'STEPS_PRESENT',
                                            'REACTION',
                                            'GENES',
                                            'PATHWAY_DEFINITION'])
    
    for p in metadata.pathways:
        results = p.score_pathway(genes, genome.short_name)
        pathway_results = pd.concat([pathway_results, results.to_frame().T], ignore_index=True)

    pathway_results.to_csv(f"{genome.name}_gator_pathway.txt", sep="\t", index=False)

    # Clean Up
    genome.remove_temp()



# MAIN ########################################################################

if __name__ == "__main__":

    if args.verify_db:
        print("Verifying db only")
        metadata.remove_temp_files()

    else:
        
        # # save raw data
        # file_out_paths = {}
        # if args.save_raw:
        #     for id, db in metadata.db_info.iterrows():
        #         file_out_paths[db['DB_NAME']] = open(db['DB_NAME'] + ".txt", "w")

        if args.patric:
            print(f'{colors.bcolors.GREEN}Patric mode enabled{colors.bcolors.END}')

        # get genomes, extract .faa file from genbank files
        unannotated_genomes = utilities.get_files(
            args.files, args.in_dir, ["faa", "gb", "gbk", "gbff"])

        print(f'{colors.bcolors.GREEN}Starting GATOR{colors.bcolors.END}')
        gator_pool = Pool()

        for _ in tqdm(gator_pool.imap_unordered(annotate, unannotated_genomes), total=len(unannotated_genomes), desc="Annotated", unit=" genomes"):
            pass
        gator_pool.close()

        # cleanup
        metadata.remove_temp_files()
        print(f'{colors.bcolors.GREEN}Thanks for using the Genome annotATOR!{colors.bcolors.END}')
