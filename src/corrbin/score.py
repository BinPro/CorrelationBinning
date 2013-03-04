import pandas as p
import corrbin.classification as c
from scipy.stats.mstats import zscore
import numpy as np
import corrbin.exceptions as exc

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
    """Making sense out of the experiment data"""

    def __init__(self, roc_axis_funs):
        self.roc_axis_funs = roc_axis_funs
        self.df = None
        self.classification = {}
        self.roc_data = {}

    def load_data_frame(self,input_file):
        self.df = p.io.parsers.read_table(input_file, sep='\t')
    
    def standardize(self):
        no_inf_df = self.df.replace(-np.inf,np.nan)
        self.df['p_value_standardized'] = \
            p.Series(no_inf_df.p_value, index=self.df.index)
        for contig_id in self.df.contig_id.unique():
            cond = self.df.contig_id == contig_id
            self.df.p_value_standardized[cond] =\
                p.Series(zscore(no_inf_df.p_value[cond]), index=self.df.index[cond])

    def classify(self, level):
        try:
            df_class = \
                c.classify_bool(self,level)
            self.classification[level] = df_class
        except exc.LevelError as e:
            print "Wrong value of level: ", e.level

    def calculate_roc(self):
        if len(self.classification)>0:
            for level in self.classification.keys():
                self.roc_data[level] = \
                    c.calculate_roc(self.classification[level],\
                                        self.roc_axis_funs)
        else:
            pass


class RocAxisFuns(object):
    FUNS = ["true_positive_rate","false_positive_rate",\
                "precision", "false_discovery_rate", \
                "included_contigs_ratio", "accuracy"]
    def __init__(self, x_fun_name, y_fun_name):
        if (x_fun_name in self.FUNS) and (y_fun_name in self.FUNS):
            self.x_fun = eval("c." + x_fun_name)
            self.y_fun = eval("c." + y_fun_name)
        else:
            self.x_fun = None
            self.y_fun = None
            print "Please specify a valid function name."

            
