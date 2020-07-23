import pandas as pd


class Annotation:
    def __init__(self, pandas_row):
        self.gene = pandas_row['GENE_NAME']
        self.product = pandas_row['GENE_PRODUCT']
        self.note = pandas_row['GENE_NOTE']
        self.reaction = pandas_row['REACTION']
        self.complex = pandas_row['COMPLEX']
        self.db_hits = {}

    def __str__(self):
        return f"<GATOR Annotation Class> {self.gene}"
