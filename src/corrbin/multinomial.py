from copy import copy

from probin.dna import DNA
from probin.model.composition import multinomial as mn

from corrbin.score import ScoreCollection, Score
from corrbin.contig_generation import SampleSetting,sample_contig


def par_subtraction(signature, contig):
    signature.subtract(contig.signature)
    return signature_to_parameters(signature)

def signature_to_parameters(signature):
    pa = {}
    n = sum(signature.values())
    for i in xrange(DNA.kmer_hash_count):
        pa[i] = signature[i]/float(n)
    return pa

def score_contig(contig, par, genome, genome_index, test):
    score_obj = ScoreCollection()
    score_gen = Score(mn.log_probability(contig.signature, par), genome, genome, contig.contig_id)
    score_obj.genome = score_gen
    group_genomes = test.group.all_genomes_but_index(genome_index)
    for gen in group_genomes:
        score_obj.group.append(Score(mn.log_probability(contig.signature, gen.par()),genome, gen, contig.contig_id))
    for group in test.rest_groups:
        outside_group = []
        for gen in group.genomes:
            score = Score(mn.log_probability(contig.signature,gen.par()), genome, gen, contig.contig_id)
            outside_group.append(score)
        score_obj.other.append(outside_group)
    return score_obj

# Quite ugly way of adding a new function to DNA class
# Will be added as soon as any other function from this
# module is loaded...
def par(self):
    try:
        return self.parameters
    except:
        self.parameters = signature_to_parameters(self.signature)
        return self.parameters

DNA.par = par

class Experiment(object):
    def __init__(self,x_st, group, rest_groups, id_gen):
        self.n = len(group)
        self.x_st = x_st
        self.count_per_g = self.count_per_g(x_st.no_contigs)
        self.group = group
        self.rest_groups = rest_groups
        self.id_gen = id_gen

    def count_per_g(self, no_contigs):
        if self.x_st.prio == "genomes":
            return [round(no_contigs/float(self.n))]*self.n
        elif self.x_st.prio == "groups":
            base = no_contigs/self.n
            count_l = [base]*self.n
            list_index = 0
            while sum(count_l)<no_contigs:
                try:
                    count_l[list_index] += 1
                    list_index += 1
                except:
                    list_index = 0
            return count_l

    def execute(self):
        all_scores = []
        for genome_index in range(self.n):
            score_objs = self.test_genome(int(genome_index))
            all_scores.append(score_objs)
        return all_scores

    def test_genome(self, genome_index):
        score_objs = []
        i = 0
        genome = self.group.genomes[genome_index]
        while i < self.count_per_g[genome_index]:
            c = sample_contig(genome, self.x_st, self.id_gen.id())
            c.calculate_signature()
            new_par = par_subtraction(copy(genome.signature),c) 

            score_obj = score_contig(c, new_par, genome, genome_index, self)
            score_objs.append(score_obj)
            i+=1
        return score_objs
