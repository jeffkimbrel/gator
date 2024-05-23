import sys
import os
import argparse
import pandas as pd
from multiprocessing import Pool
from tqdm import tqdm
import ast
import datetime
import pickle

import metadata as m
from jakomics import utilities, kegg, colors, blast, hmm
from jakomics.genome import GENOME
from jakomics.patric import Patric_Gene

# OPTIONS #####################################################################

parser = argparse.ArgumentParser(description="",
                                 formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('--in_dir',
                    help="Directory with faa genomes",
                    required=False,
                    default="")

parser.add_argument('--out_dir',
                    help="Directory to write results to",
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

parser.add_argument('--file_list',
                    help="File with full paths",
                    required=False,
                    default = "default")

parser.add_argument('--verify_db',
                    action='store_true',
                    help='Just check the database')

parser.add_argument('--docker',
                    action='store_true',
                    help='Run gator in a docker container')

parser.add_argument('-p', '--patric',
                    action='store_true',
                    help='Genbank files are from patric. Adds patric and EC db support')

parser.add_argument('--save_raw',
                    action='store_true',
                    help='Save the raw search data to files')

parser.add_argument('--score_as_ratio',
                    action='store_true',
                    help='Return KOFAM scores as ratios of the treshold')

parser.add_argument('-n', '--name',
                    help="Name to prepend the output files. Defaults to a timestamp",
                    required=False,
                    default=datetime.datetime.now().strftime("%Y%m%d_%H%M%S"))

parser.add_argument('-t', '--threads',
                    help="Threads to use",
                    type = int,
                    required=False,
                    default=8)

parser.add_argument('--pickle',
                    help="Path to save gator metadata pickle file",
                    required=False,
                    default=f"gator_metadata.pickle")

args = parser.parse_args()

if args.docker:
    print("gator is running in docker mode...")
    args.out_dir = "/gator-app"

# Functions

def prep_metadata(version, genome_paths):
    # prep metadata and databases
    
    if args.gator_db == "default":
        gator_path = (os.path.dirname(os.path.abspath(sys.argv[0])))
        metadata = m.Metadata(os.path.join(gator_path, "gator_db.xlsx"))
    else:
        metadata = m.Metadata(args.gator_db)

    metadata.create_hal_files(args.out_dir)
    metadata.make_blast_dbs()
    metadata.verify_metadata()

    metadata.summary()

    if args.verify_db == False:
        metadata.version = version
        metadata.genome_paths = genome_paths
        file = open(os.path.join(args.out_dir, args.pickle), 'wb')
        pickle.dump(metadata, file)
        file.close()

    return(metadata)

def annotate(genome):

    file = open(os.path.join(args.out_dir, args.pickle), 'rb')
    metadata = pickle.load(file)
    file.close()

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

            # if utilities.check_executable("exec_annotation"):
            #     print("exec_annotation found!")

            hits = kegg.run_kofam(faa_path = genome.file_path, 
                                  hal_path = db['hal_path'], 
                                  temp_dir = os.path.join(args.out_dir, genome.id + "_KO"),
                                  ko_list = os.path.join(db['DB_PATH'], 'ko_list'), 
                                  score_as_ratio = args.score_as_ratio,
                                  echo = False)
            
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

        # elif db['METHOD'] == 'EC':
        #     if args.patric:
        #         for gene in genome.patric:
        #             if hasattr(genome.patric[gene], 'EC_number'):
        #                 patric_annotation = Patric_Gene(genome.patric[gene].id)
        #                 for EC in genome.patric[gene].EC_number:
        #                     patric_annotation.annotation = EC

        #                     if EC in genome.raw_results[db['DB_NAME']]:
        #                         genome.raw_results[db['DB_NAME']][EC].append(
        #                             patric_annotation)
        #                     else:
        #                         genome.raw_results[db['DB_NAME']][EC] = [
        #                             patric_annotation]

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

    details_file = os.path.join(args.out_dir, f"{genome.short_name}_gator_{args.name}_detail.txt")
    details.to_csv(details_file, sep="\t", index=False)


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

    pathway_file = os.path.join(args.out_dir, f"{genome.short_name}_gator_{args.name}_pathway.txt")
    pathway_results.to_csv(pathway_file, sep="\t", index=False)

    # Clean Up
    genome.remove_temp()

# MAIN ########################################################################

if __name__ == "__main__":
    
    version = "v1.6.2"

    print(f'{colors.bcolors.GREEN}Genome annotATOR (GATOR) {version}{colors.bcolors.END}')

    if args.verify_db:
        print(f"{colors.bcolors.PURPLE}Verifying db only{colors.bcolors.END}")
        metadata = prep_metadata(version, "")
        metadata.remove_temp_files()

    else:
        if args.score_as_ratio:
            print(f"{colors.bcolors.YELLOW}WARNING: Gator is running in `score_as_ratio` mode, so all kofam results will be returned. Only kofam scores with a value >= 1 would be considered `passed` with this mode disabled{colors.bcolors.END}")

        if args.file_list != "default":
            with open(args.file_list) as file:
                file_list = [line.rstrip() for line in file]
        else:
            file_list = args.files

        # get genomes, extract .faa file from genbank files
        unannotated_genomes = utilities.get_files(
            file_list, args.in_dir, ["faa", "gb", "gbk", "gbff"])


        metadata = prep_metadata(version, [[o.short_name, o.file_path, o.id] for o in unannotated_genomes]) # gator version and list of genomes for pickle file

        # # save raw data
        # file_out_paths = {}
        # if args.save_raw:
        #     for id, db in metadata.db_info.iterrows():
        #         file_out_paths[db['DB_NAME']] = open(db['DB_NAME'] + ".txt", "w")

        if args.patric:
            print(f'{colors.bcolors.GREEN}Patric mode enabled{colors.bcolors.END}')

        # decide number of threads
        threads = args.threads
        if len(unannotated_genomes) < threads:
            threads = len(unannotated_genomes)
        
        if threads > os.cpu_count():
            threads = os.cpu_count()

        print(f'{colors.bcolors.GREEN}Starting GATOR with {threads} thread(s){colors.bcolors.END}')
        gator_pool = Pool(threads)

        for _ in tqdm(gator_pool.imap_unordered(annotate, unannotated_genomes), total=len(unannotated_genomes), desc="Annotated", unit=" genomes"):
            pass
        gator_pool.close()

        # cleanup
        metadata.remove_temp_files()
        print(f'{colors.bcolors.GREEN}Thanks for using the Genome annotATOR!{colors.bcolors.END}')
