import os
import sys
import uuid

import pandas as pd
import numpy as np
from jakomics import colors
import pathway


class Metadata:
    def __init__(self, metadata_file_path):
        self.db_info = pd.read_excel(metadata_file_path, sheet_name="db")
        self.gene_info = pd.read_excel(metadata_file_path, sheet_name="gene")
        self.pathway_info = pd.read_excel(metadata_file_path, sheet_name="pathway")

        # remove empty rows
        self.db_info.dropna(subset=["DB_NAME"], inplace=True)
        self.gene_info.dropna(subset=["GENE_NAME"], inplace=True)
        self.pathway_info.dropna(subset=["PATHWAY_NAME"], inplace=True)

        self.parse_paths()

    def __str__(self):
        return "<GATOR Metadata Class>"

    def remove_temp_files(self):
        for id, db in self.db_info.iterrows():
            os.remove(db['hal_path'])

    def parse_paths(self):
        self.pathways = []
        for id, path in self.pathway_info.iterrows():
            self.pathways.append(pathway.Pathway(path))

    def create_hal_files(self):
        '''
        finds dbs of method kofam and makes a .hal file of all KOs under that db name on the gene sheet.
        Adds hal path to hal_path line in db_info
        '''

        print(f'Making .hal files')

        self.db_info['hal_path'] = ""

        for id, db in self.db_info.iterrows():

            if db['METHOD'] == 'kofam':
                temp_file = uuid.uuid4().hex + ".hal"
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
