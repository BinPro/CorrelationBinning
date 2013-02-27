import pandas as p
from scipy.stats.mstats import zscore
import corrbin.classification as c

class ScoreCollection(object):
    def __init__(self):
        self.genome = None
        self.group = []
        self.other = []
    def __str__(self):
        string = str(self.genome)
        for score in self.group:
            string += "\n" + str(score)
        for other_group in self.other:
            for score in other_group:
                string += "\n" + str(score)
        return string
        

class Score(object):
    def __init__(self, p_value, contig_genome, compare_genome, contig_id):
        self.p_value = p_value
        self.genome_contig = contig_genome.id
        self.species_contig = contig_genome.species
        self.genus_contig = contig_genome.genus
        self.family_contig = contig_genome.family
        self.genome_compare = compare_genome.id
        self.species_compare = compare_genome.species
        self.genus_compare = compare_genome.genus
        self.family_compare = compare_genome.family
        self.contig_id = contig_id

    def __str__(self):
        return "%f\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%i" % (self.p_value,self.family_contig, self.genus_contig, self.species_contig, self.genome_contig, self.family_compare, self.genus_compare, self.species_compare, self.genome_compare, self.contig_id)


class ExperimentData(object):
    def __init__(self):
        self.df = None
        self.classification = None
        
    def load_data_frame(self,input_file):
        self.df = p.io.parsers.read_table(input_file, sep='\t')
    
    def standardize(self):
        self.df['p_value_standardized'] = p.Series(zscore(self.df.p_value), index=self.df.index)

    def classify(self):
        df_class = c.classify_bool(self)
        self.classification = df_class

