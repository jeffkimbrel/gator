import pandas as pd


class Pathway:

    def __init__(self, pandas_raw):
        self.name = pandas_raw['PATHWAY_NAME']
        self.definition = pandas_raw['DEFINITION']
        self.compounds = [x.strip() for x in pandas_raw['COMPOUNDS'].split('->')]
        self.steps = self.parse_definition()

    def __str__(self):
        return f"<GATOR Pathway Class> {self.name}"

    def parse_definition(self):
        steps = []
        step_list = [x.strip() for x in self.definition.split('->')]
        for step in step_list:
            complexes_notparsed = [x.strip() for x in step.split('|')]
            complexes_parsed = []
            for complex in complexes_notparsed:
                genes = set([x.strip() for x in complex.split('+')])
                complexes_parsed.append(list(genes))

            steps.append(complexes_parsed)

        return steps

    def score_pathway(self, annotations, genome="generic"):
        steps = len(self.steps)
        steps_count = 0
        steps_present = []
        genes_present = []
        for i, complexes in enumerate(self.steps):
            step_present = False
            genes_present_step = []
            step = i + 1
            for genes in complexes:
                genes_present_step.append(list(set(annotations) & set(genes)))
                if set(genes).issubset(set(annotations)):
                    step_present = True
            if step_present:
                steps_count += 1
            steps_present.append(step_present)

            genes_present.append(genes_present_step)

        # make fancy compound string
        reaction_string = ""
        for i, step in enumerate(steps_present):
            if step:
                arrow = " ---> "
            else:
                arrow = " -x-> "
            reaction_string += self.compounds[i] + arrow
        reaction_string += self.compounds[-1]

        # score entire pathway present or absent
        if steps == steps_count:
            pathway_present = True
        else:
            pathway_present = False

        # print(f'{self.name}\t{pathway_present}\t{steps}\t{steps_count}\t{reaction_string}\t{genes_present}\t{self.definition}')

        df = pd.Series(data={'GENOME': genome,
                             'PATHWAY': self.name,
                             'PRESENT': pathway_present,
                             'PATHWAY_STEPS': steps,
                             'STEPS_PRESENT': steps_count,
                             'REACTION': reaction_string,
                             'GENES': genes_present,
                             'PATHWAY_DEFINITION': self.definition
                             }
                       )
        return df
