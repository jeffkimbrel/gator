import os
import sys
import uuid

import pandas as pd
import numpy as np
from jakomics import colors, blast
import pathway


class Metadata:
    def __init__(self, metadata_file_path):
        self.db_info = pd.read_excel(metadata_file_path,
                                     sheet_name="db",
                                     engine='openpyxl')
        self.gene_info = pd.read_excel(metadata_file_path,
                                       sheet_name="gene",
                                       engine='openpyxl')
        self.pathway_info = pd.read_excel(metadata_file_path,
                                          sheet_name="pathway",
                                          engine='openpyxl')

        # remove empty rows
        self.db_info.dropna(subset=["DB_NAME"], inplace=True)
        self.gene_info.dropna(subset=["GENE_NAME"], inplace=True)
        self.pathway_info.dropna(subset=["PATHWAY_NAME"], inplace=True)

        self.parse_paths()

    def __str__(self):
        return "<GATOR Metadata Class>"

    def summary(self):
        print(
            f"{colors.bcolors.GREEN}GATOR DB GENES: {len(self.gene_info)}\nGATOR DB PATHWAYS: {len(self.pathway_info)}{colors.bcolors.END}")

    def verify_metadata(self):

        # are all genes unique?
        genes = self.gene_info['GENE_NAME']
        if self.gene_info['GENE_NAME'].is_unique:
            print("All GENE_NAME values are unique.")
            self.genes = genes
        else:
            print(
                f"{colors.bcolors.RED}ERROR: Duplicate GENE_NAMEs found: {set(genes[genes.duplicated(keep=False)])}{colors.bcolors.END}")
            sys.exit()

        # are all pathway genes actually being searched for?
        all_genes_searched_for = True
        for p in self.pathways:
            if set(p.genes).issubset(self.genes):
                pass
            else:
                all_genes_searched_for = False
                print(
                    f"{colors.bcolors.RED}ERROR: {p.name} GENE_NAMEs ({list(set(p.genes) - set(self.genes))}) are not found on the genes sheet.{colors.bcolors.END}")

        if all_genes_searched_for:
            print("All GENE_NAMEs in the pathway search are found on the genes sheet.")
        else:
            print(
                f"{colors.bcolors.RED}ERROR: Problems found... please check the gene names on the Excel pathway sheet{colors.bcolors.END}")
            sys.exit()

    def remove_temp_files(self):
        for id, db in self.db_info.iterrows():
            if pd.notnull(db['hal_path']):
                os.remove(db['hal_path'])

    def parse_paths(self):
        self.pathways = []
        for id, path in self.pathway_info.iterrows():
            self.pathways.append(pathway.Pathway(path))

    def make_blast_dbs(self):
        for id, db in self.db_info.iterrows():
            if db['METHOD'] == 'blastp':
                print(f"Making blast database at {db['DB_PATH']}")
                blast.make_blast_db('prot', db['DB_PATH'])

    def create_hal_files(self, out_dir):
        '''
        finds dbs of method kofam and makes a .hal file of all KOs under that db name on the gene sheet.
        Adds hal path to hal_path line in db_info
        '''

        self.db_info['hal_path'] = None

        for id, db in self.db_info.iterrows():

            if db['METHOD'] == 'kofam':
                temp_file = os.path.join(out_dir, uuid.uuid4().hex + ".hal")
                print(f'Making .hal file at {temp_file}')
                hal_target = open(temp_file, 'w')
                self.db_info.at[id, 'hal_path'] = temp_file

                found_problem = False

                # the ko_raw_list to ko_list deals with comma-separated list of KOs in one cell
                for ko_raw_list in self.gene_info[db['DB_NAME']]:
                    if pd.notnull(ko_raw_list):
                        ko_list = [x.strip() for x in ko_raw_list.split(',')]
                        for ko in ko_list:
                            hmm_path = os.path.join(db['DB_PATH'], ko + ".hmm")

                            if os.path.exists(hmm_path):
                                hal_target.write(hmm_path + "\n")
                            else:
                                found_problem = True
                                print(
                                    f'WARNING: {colors.bcolors.RED}{ko}.hmm{colors.bcolors.END} is not found at {db["DB_PATH"]}')

                hal_target.close()
                
                if found_problem:
                    print("There were some issues here...")
                    os.remove(temp_file)
                    sys.exit()
                else:
                    print(f'All .hmm files are found in {db["DB_PATH"]}')
